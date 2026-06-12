#!/usr/bin/env bash
# =============================================================================
# 01_download.sh — fetch scRNA-seq datasets for the pulmonary fibrosis pipeline
# Run from the repo root:  bash scripts/01_download.sh
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$ROOT/data"

mkdir -p "$DATA"/human/{adams_gse136831,habermann_gse135893} \
         "$DATA"/mouse/{strunz_gse141259,kobayashi_gse141634,choi_gse145031} \
         "$ROOT"/results/{qc,markers,nichenet,cellchat,geneformer,figures}

# GEO supplementary RAW tarballs follow this URL pattern:
#   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE<NNN>nnn/GSE<FULL>/suppl/GSE<FULL>_RAW.tar
geo_raw () {  # $1 = full accession (e.g. GSE141259)
  local acc="$1" stem
  stem="${acc:0:$(( ${#acc} - 3 ))}nnn"   # GSE141259 -> GSE141nnn
  echo "https://ftp.ncbi.nlm.nih.gov/geo/series/${stem}/${acc}/suppl/${acc}_RAW.tar"
}

fetch_tar () {  # $1 = accession, $2 = target dir
  local acc="$1" dir="$2" url
  url="$(geo_raw "$acc")"
  echo ">> $acc  ->  $dir"
  if [ -f "$dir/${acc}_RAW.tar" ]; then
    echo "   already downloaded, skipping"
  else
    wget -c -P "$dir" "$url"
  fi
  echo "   extracting…"
  tar -xvf "$dir/${acc}_RAW.tar" -C "$dir" >/dev/null
  # many GEO files are gzipped after extraction
  find "$dir" -name '*.gz' -exec gunzip -kf {} + 2>/dev/null || true
}

# -----------------------------------------------------------------------------
# MOUSE  (primary engine of the pipeline — time course + transitional states)
# -----------------------------------------------------------------------------
echo "=== MOUSE datasets ==="
fetch_tar GSE141259 "$DATA/mouse/strunz_gse141259"     # Strunz bleomycin time course
fetch_tar GSE141634 "$DATA/mouse/kobayashi_gse141634"  # Kobayashi PATS
fetch_tar GSE145031 "$DATA/mouse/choi_gse145031"       # Choi DATPs

# -----------------------------------------------------------------------------
# HUMAN  (disease-state anchor + translation)
# Note: Adams/Habermann RAW tars are large. Processed objects are often easier:
#   - IPF Cell Atlas portal:   https://www.ipfcellatlas.com   (Adams = Kaminski/Yale)
#   - Habermann processed:     Broad Single Cell Portal SCP890 / GEO GSE135893 suppl
# Try RAW here; if you only need processed objects, comment these out and use portals.
# -----------------------------------------------------------------------------
echo "=== HUMAN datasets ==="
fetch_tar GSE136831 "$DATA/human/adams_gse136831"      # Adams (Sci Adv 2020)
fetch_tar GSE135893 "$DATA/human/habermann_gse135893"  # Habermann (Sci Adv 2020)

# -----------------------------------------------------------------------------
# HEALTHY HUMAN REFERENCE (HLCA, Sikkema 2023) — download manually from CELLxGENE
# -----------------------------------------------------------------------------
cat <<'EOF'

[manual] Human Lung Cell Atlas (healthy AT2/AT1 reference):
  Download the .h5ad from CELLxGENE and place it in data/human/hlca/
    https://cellxgene.cziscience.com/   (search "Human Lung Cell Atlas")

EOF

echo "Done. Contents of data/:"
find "$DATA" -maxdepth 2 -type d | sort
