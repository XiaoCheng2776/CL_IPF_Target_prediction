#!/usr/bin/env python3
"""
06_geneformer.py — Geneformer fine-tuning + in silico perturbation (mouse)
Run from repo root:  conda run -n gf python scripts/06_geneformer.py

Uses Geneformer V1-10M (6-layer, 256-hidden, 25k vocab).
Mouse genes are mapped to human 1:1 orthologs (ENSG IDs) for tokenisation —
required because Geneformer V1 was pre-trained on human scRNA-seq.

GPU:   NVIDIA RTX A2000 12 GB — fine-tuning runs at batch_size=4 with fp16.
       Reduce forward_batch_size if OOM during perturbation inference.

Input:  results/markers/epi_obj.rds       (from 03_phenotype.R)
Output: results/geneformer/gf_lung/       (fine-tuned model)
        results/geneformer/oe/receptors_ranked_oe.csv
        results/geneformer/del/receptors_ranked_del.csv
"""

import os
import logging
import pickle
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd
import requests
import scanpy as sc
import scipy.io
import torch
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

# =============================================================================
# Paths
# =============================================================================
ROOT       = Path(__file__).resolve().parents[1]
DATA_DIR   = ROOT / "data" / "mouse"
TOKEN_DIR  = ROOT / "results" / "geneformer" / "tokenized"
MODEL_DIR  = ROOT / "results" / "geneformer" / "gf_lung"
OE_DIR     = ROOT / "results" / "geneformer" / "oe"
DEL_DIR    = ROOT / "results" / "geneformer" / "del"
for d in [DATA_DIR, TOKEN_DIR, MODEL_DIR, OE_DIR, DEL_DIR]:
    d.mkdir(parents=True, exist_ok=True)

GMAP_PATH     = DATA_DIR / "gene_symbol_to_ensembl_mouse.csv"  # mouse sym → human ENSG
RAW_H5AD      = DATA_DIR / "epi_raw.h5ad"                      # raw counts AnnData
GF_H5AD       = DATA_DIR / "epi_gf.h5ad"                       # ortholog-remapped h5ad
TOK_DATASET   = TOKEN_DIR / "lung_epi.dataset"
PRETRAIN_ID   = "ctheodoris/Geneformer"
PRETRAIN_SUB  = "Geneformer-V1-10M"

# V1 token dictionary (gc30M = 30M-cell pretrain vocabulary, 25 426 human genes)
GF_PKG    = Path(__file__).resolve().parents[1] / ".."   # resolved below
import geneformer as _gf
GF_PKG    = Path(_gf.__file__).parent
TOK_DICT_PATH = GF_PKG / "gene_dictionaries_30m" / "token_dictionary_gc30M.pkl"

# Candidate receptors for in silico perturbation (mouse gene symbols)
CANDIDATE_RECEPTORS = [
    "Fgfr2", "Fgfr1",           # FGFR2b axis — positive control
    "Egfr", "Erbb2", "Erbb3", "Met", "Igf1r",
    "Axl",                       # pro-fibrotic: expect negative regen_gain
    "Bmpr1a", "Bmpr2", "Acvrl1",
    "Tgfbr1", "Tgfbr2",          # TGFb — pro-fibrotic internal negative control
    "Pdgfra", "Pdgfrb",
    "Notch1", "Notch2",
    "Fzd1", "Fzd2",
]

# =============================================================================
# 0.  Download mouse → human 1:1 ortholog gene map (Ensembl BioMart)
# =============================================================================
if not GMAP_PATH.exists():
    log.info("Downloading mouse→human ortholog map from Ensembl BioMart …")
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<!DOCTYPE Query>'
        '<Query virtualSchemaName="default" formatter="TSV" header="1"'
        ' uniqueRows="0" count="" datasetConfigVersion="0.6">'
        '<Dataset name="mmusculus_gene_ensembl" interface="default">'
        '<Filter name="with_hsapiens_homolog" excluded="0"/>'
        '<Attribute name="external_gene_name"/>'
        '<Attribute name="hsapiens_homolog_ensembl_gene"/>'
        '<Attribute name="hsapiens_homolog_orthology_type"/>'
        '</Dataset></Query>'
    )
    r = requests.get(
        "https://www.ensembl.org/biomart/martservice",
        params={"query": xml},
        timeout=300,
    )
    r.raise_for_status()
    df_map = pd.read_csv(pd.io.common.StringIO(r.text), sep="\t")
    df_map.columns = ["symbol", "human_ensembl", "orthology_type"]
    df_map = df_map[
        (df_map["orthology_type"] == "ortholog_one2one")
        & df_map["human_ensembl"].notna()
        & df_map["human_ensembl"].str.startswith("ENSG")
    ].drop_duplicates("symbol")[["symbol", "human_ensembl"]].rename(
        columns={"human_ensembl": "ensembl_id"}
    )
    df_map.to_csv(GMAP_PATH, index=False)
    log.info(f"  Saved {len(df_map)} 1:1 ortholog pairs → {GMAP_PATH}")

gmap    = pd.read_csv(GMAP_PATH)
sym2ens = dict(zip(gmap["symbol"], gmap["ensembl_id"]))
ens2sym = {v: k for k, v in sym2ens.items()}
log.info(f"Gene map: {len(sym2ens)} mouse → human 1:1 orthologs")

# =============================================================================
# 1.  Export epi_obj from R → sparse-matrix files → AnnData h5ad
#     (R Seurat → Python AnnData without SeuratDisk/anndataR dependencies)
# =============================================================================
if not RAW_H5AD.exists():
    mtx_path  = DATA_DIR / "epi_counts.mtx"
    bar_path  = DATA_DIR / "epi_barcodes.csv"
    feat_path = DATA_DIR / "epi_features.csv"
    meta_path = DATA_DIR / "epi_meta.csv"

    if not mtx_path.exists():
        log.info("Exporting epi_obj RNA counts from R …")
        r_code = r"""
suppressPackageStartupMessages(library(Seurat)); library(Matrix)
obj <- readRDS("results/markers/epi_obj.rds")
DefaultAssay(obj) <- "RNA"; obj <- JoinLayers(obj)
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
writeMM(counts, "data/mouse/epi_counts.mtx")
write.csv(data.frame(barcode = colnames(counts)), "data/mouse/epi_barcodes.csv", row.names = FALSE)
write.csv(data.frame(gene   = rownames(counts)), "data/mouse/epi_features.csv", row.names = FALSE)
write.csv(obj@meta.data, "data/mouse/epi_meta.csv")
message("R export: ", ncol(counts), " cells  ", nrow(counts), " genes")
"""
        rscript = shutil.which("Rscript") or "/home/shiqi/miniforge3/envs/ipf/bin/Rscript"
        proc = subprocess.run(
            [rscript, "--vanilla", "-"],
            input=r_code, capture_output=True, text=True,
        )
        if proc.returncode != 0:
            log.error(proc.stderr[-2000:])
            raise RuntimeError("R export failed")
        log.info(proc.stderr.strip().split("\n")[-1])

    log.info("Building AnnData from R export …")
    counts   = scipy.io.mmread(str(mtx_path)).T.tocsr()   # cells × genes
    barcodes = pd.read_csv(bar_path)["barcode"].values
    features = pd.read_csv(feat_path)["gene"].values
    meta     = pd.read_csv(meta_path, index_col=0)
    # Align obs rows to barcode order
    meta = meta.reindex(barcodes)
    adata = sc.AnnData(
        X   = counts,
        obs = meta,
        var = pd.DataFrame(index=features),
    )
    adata.obs_names = barcodes
    adata.var_names = features
    log.info(f"  Raw AnnData: {adata.n_obs} × {adata.n_vars}")
    adata.write_h5ad(RAW_H5AD)

adata = sc.read_h5ad(RAW_H5AD)
log.info(f"Loaded raw AnnData: {adata.n_obs} × {adata.n_vars}")

# =============================================================================
# 2.  Apply fate labels and filter to completing + arrested cells
# =============================================================================
assert "fate_label" in adata.obs.columns, "fate_label missing — rerun 03_phenotype.R"
adata = adata[adata.obs["fate_label"].isin(["completing", "arrested"])].copy()
adata.obs["regen_label"] = (
    adata.obs["fate_label"].map({"completing": 1, "arrested": 0}).astype(int)
)
log.info(
    f"Cells after fate filter: {adata.n_obs}  "
    f"(completing={int((adata.obs.regen_label == 1).sum())}, "
    f"arrested={int((adata.obs.regen_label == 0).sum())})"
)

# =============================================================================
# 3.  Map mouse symbols → human ENSG IDs; keep only genes in Geneformer vocab
# =============================================================================
with open(TOK_DICT_PATH, "rb") as f:
    token_dict = pickle.load(f)
valid_ensg = set(token_dict.keys())

adata.var["ensembl_id"] = adata.var_names.map(sym2ens)
before = adata.n_vars
adata  = adata[
    :, adata.var["ensembl_id"].notna() & adata.var["ensembl_id"].isin(valid_ensg)
].copy()
adata.var_names = adata.var["ensembl_id"].values
log.info(f"Genes with human ENSG in Geneformer V1 vocab: {adata.n_vars}/{before}")

# Candidate receptors available in vocab
cand_sym_to_ens = {
    s: sym2ens[s]
    for s in CANDIDATE_RECEPTORS
    if s in sym2ens and sym2ens[s] in valid_ensg
}
log.info(
    f"Candidate receptors in vocab: {len(cand_sym_to_ens)}/{len(CANDIDATE_RECEPTORS)}  "
    f"({sorted(cand_sym_to_ens)})"
)

# Geneformer needs n_counts per cell; use raw count totals
adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).flatten()

if not GF_H5AD.exists():
    adata.write_h5ad(GF_H5AD)

# =============================================================================
# 4.  Tokenise (TranscriptomeTokenizer V1)
#     Ranks genes by expression × global median weight, no CLS/EOS tokens.
# =============================================================================
if not TOK_DATASET.exists():
    log.info("Tokenising …")
    from geneformer import TranscriptomeTokenizer
    tk = TranscriptomeTokenizer(
        # "label" is the column name required by Classifier for cell classification
        custom_attr_name_dict={"regen_label": "label", "cell_type": "cell_type"},
        nproc=4,
        model_version="V1",
    )
    tok_tmp = tempfile.mkdtemp()
    shutil.copy(str(GF_H5AD), os.path.join(tok_tmp, "epi_gf.h5ad"))
    tk.tokenize_data(tok_tmp, str(TOKEN_DIR), "lung_epi", file_format="h5ad")
    shutil.rmtree(tok_tmp)
    log.info(f"Tokenised → {TOK_DATASET}")
else:
    log.info("Tokenised dataset found; skipping.")

# =============================================================================
# 5.  Fine-tune: completing (label=1) vs arrested (label=0)
#     Uses Geneformer's Classifier wrapper (correct collator + pad logic for V1).
# =============================================================================
from datasets import load_from_disk

ds_all = load_from_disk(str(TOK_DATASET))
log.info(f"Tokenised dataset: {len(ds_all)} cells, columns: {ds_all.column_names}")
ds     = ds_all.train_test_split(test_size=0.2, seed=42)

SAVED_MODEL = MODEL_DIR / "pytorch_model.bin"
if not SAVED_MODEL.exists():
    from geneformer.classifier import Classifier

    clf = Classifier(
        classifier      = "cell",
        cell_state_dict = {
            "state_key": "label",   # column in tokenized dataset (0=arrested, 1=completing)
            "states"   : "all",     # use all unique label values
        },
        model_version= "V1",
        training_args= {
            "num_train_epochs"             : 5,
            "per_device_train_batch_size"  : 4,    # 12 GB GPU
            "per_device_eval_batch_size"   : 8,
            "fp16"                         : True,
            "save_strategy"                : "epoch",
            "load_best_model_at_end"       : True,
            "metric_for_best_model"        : "macro_f1",
            "logging_steps"                : 50,
        },
        nproc = 4,
        ngpu  = 1,
    )
    log.info("Starting fine-tuning (V1-10M → completing vs arrested) …")
    # Download V1-10M from HuggingFace on first run
    from huggingface_hub import snapshot_download
    pretrain_local = ROOT / "data" / "mouse" / "geneformer_v1_10m"
    if not pretrain_local.exists():
        log.info("  Downloading Geneformer-V1-10M weights …")
        snapshot_download(
            repo_id      = PRETRAIN_ID,
            allow_patterns=[f"{PRETRAIN_SUB}/*"],
            local_dir    = pretrain_local,
        )
    pretrain_dir = pretrain_local / PRETRAIN_SUB
    clf.train_classifier(
        model_directory  = str(pretrain_dir),
        num_classes      = 2,
        train_data       = ds["train"],
        eval_data        = ds["test"],
        output_directory = str(MODEL_DIR),
    )
    log.info("Fine-tuning complete.")
else:
    log.info(f"Fine-tuned model found at {MODEL_DIR}; skipping training.")

# =============================================================================
# 6.  In silico perturbation — direct classifier probability shift
#
#     OE  = prepend gene token (highest-rank position in input sequence)
#     DEL = remove gene token from input sequence
#
#     regen_gain  (OE)  = mean[ P(completing | perturbed) - P(completing | original) ]
#     regen_loss  (DEL) = mean[ P(completing | deleted)   - P(completing | original) ]
#
#     Positive control: Fgfr2 OE → regen_gain > 0 (rescue completing fate)
#                       Fgfr2 DEL → regen_loss < 0 (Fgf10-KO equivalent)
# =============================================================================
from transformers import BertForSequenceClassification

DEVICE  = torch.device("cuda" if torch.cuda.is_available() else "cpu")
MAX_LEN = 2048
PAD_ID  = 0
INFER_BATCH = 32   # reduce if OOM during perturbation

log.info(f"Loading fine-tuned model for perturbation  (device={DEVICE}) …")
model_ft = BertForSequenceClassification.from_pretrained(str(MODEL_DIR))
model_ft = model_ft.to(DEVICE).eval()


def batch_predict_prob1(seqs: list, batch_size: int = INFER_BATCH) -> np.ndarray:
    """Return P(label=1 / completing) for each sequence in seqs."""
    probs = []
    for i in range(0, len(seqs), batch_size):
        chunk  = seqs[i : i + batch_size]
        max_l  = min(max(len(x) for x in chunk), MAX_LEN)
        padded = [list(x)[:max_l] + [PAD_ID] * (max_l - min(len(x), max_l))
                  for x in chunk]
        attn   = [[1] * min(len(x), max_l) + [0] * (max_l - min(len(x), max_l))
                  for x in chunk]
        ids_t  = torch.tensor(padded, dtype=torch.long, device=DEVICE)
        att_t  = torch.tensor(attn,   dtype=torch.long, device=DEVICE)
        with torch.no_grad():
            logits = model_ft(input_ids=ids_t, attention_mask=att_t).logits
        probs.append(torch.softmax(logits, dim=-1)[:, 1].cpu().numpy())
    return np.concatenate(probs)


# Original probabilities
orig_ids = [list(x) for x in ds_all["input_ids"]]
log.info(f"Computing baseline probabilities for {len(orig_ids)} cells …")
p_orig   = batch_predict_prob1(orig_ids)
log.info(f"  Baseline mean P(completing) = {p_orig.mean():.4f}")

results_oe  = []
results_del = []

for symbol, ensg in cand_sym_to_ens.items():
    gene_tok = token_dict.get(ensg)
    if gene_tok is None:
        log.warning(f"{symbol} ({ensg}) not in token dict; skipping")
        continue

    oe_ids  = []
    del_ids = []
    n_expressing = 0

    for ids in orig_ids:
        ids_list = list(ids)
        present  = gene_tok in ids_list
        if present:
            n_expressing += 1
            ids_list.remove(gene_tok)
        # OE: prepend to rank-1, truncate if needed
        oe_ids.append([gene_tok] + ids_list[: MAX_LEN - 1])
        # DEL: gene already removed (or absent); truncate
        del_ids.append(ids_list[:MAX_LEN])

    p_oe  = batch_predict_prob1(oe_ids)
    p_del = batch_predict_prob1(del_ids)

    gain = float(np.mean(p_oe  - p_orig))
    loss = float(np.mean(p_del - p_orig))
    results_oe.append(
        {"symbol": symbol, "ensembl_id": ensg,
         "regen_gain": gain, "n_expressing": n_expressing}
    )
    results_del.append(
        {"symbol": symbol, "ensembl_id": ensg,
         "regen_loss": loss, "n_expressing": n_expressing}
    )
    log.info(f"  {symbol:10s}  OE regen_gain={gain:+.4f}  DEL regen_loss={loss:+.4f}")

df_oe  = pd.DataFrame(results_oe).sort_values("regen_gain", ascending=False)
df_del = pd.DataFrame(results_del).sort_values("regen_loss", ascending=True)

# =============================================================================
# 7.  Rank + positive-control gate
# =============================================================================
print("\n=== POSITIVE CONTROL: Fgfr2 (in silico Fgf10-KO equivalent) ===")
fgfr2_oe  = df_oe[df_oe["symbol"]  == "Fgfr2"]
fgfr2_del = df_del[df_del["symbol"] == "Fgfr2"]
if len(fgfr2_oe):
    print(f"  OE  regen_gain : {fgfr2_oe['regen_gain'].values[0]:+.4f}   (expect > 0)")
if len(fgfr2_del):
    print(f"  DEL regen_loss : {fgfr2_del['regen_loss'].values[0]:+.4f}  (expect < 0  ≡ Fgf10-KO)")

oe_ok  = len(fgfr2_oe)  and fgfr2_oe["regen_gain"].values[0]  > 0
del_ok = len(fgfr2_del) and fgfr2_del["regen_loss"].values[0] < 0

if oe_ok and del_ok:
    print("Positive control PASSED.")
else:
    print(
        "WARNING: positive control not clean — OE ok=%s, DEL ok=%s. "
        "Classifier may need more epochs or class-balancing." % (oe_ok, del_ok)
    )

print("\nTop candidates by OE regen_gain (overexpression → completing fate):")
print(df_oe.to_string(index=False))

print("\nTop candidates by DEL regen_loss (deletion → arrested fate):")
print(df_del.to_string(index=False))

df_oe.to_csv(OE_DIR  / "receptors_ranked_oe.csv",  index=False)
df_del.to_csv(DEL_DIR / "receptors_ranked_del.csv", index=False)
log.info("06_geneformer.py complete.")
