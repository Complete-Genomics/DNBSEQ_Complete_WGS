## Build the AlphaGenome AVI lookup table used by modules/avi.nf (run once per database).
## usage: bash avi_prep.sh <alphagenome_variant_impact_score_snvs.tsv.gz> <out.tsv.gz> [threads=8]
## input columns:  #CHROM POS REF ALT raw_score PHRED
## output columns: #CHROM POS REF ALT AVI(=PHRED), bgzip + tabix
## default location expected by the pipeline: ${DB}/hg38/alphagenome/alphagenome_avi_phred.tsv.gz (or pass --avi_table)
## tabix needs the table sorted by CHROM, POS; if indexing fails, sort first (the full table is ~83G).

set -euo pipefail

in=$1
out=$2
threads=${3:-8}

mkdir -p "$(dirname "$out")"
{ printf '#CHROM\tPOS\tREF\tALT\tAVI\n'
  gzip -dc "$in" | grep -v '^#' | cut -f1-4,6
} | bgzip -@ "$threads" -c > "$out"
tabix -f -s1 -b2 -e2 "$out"
