#!/bin/bash
# Compute the CMRG report fields from an ALREADY FINISHED run (no re-alignment / calling).
#
# usage:  bash cmrg_from_results.sh <outdir> <sample> [align_tool=lariat] [var_tool=dv]
#   outdir = the --outdir of the finished run (contains <sample>/align and <sample>/phase)
# env:    DB   CWGS_db dir (default: v1.0.6 on storeData)      REF  reference fasta (needs .fai)
#
# needs bedtools + python3 with pandas and intervaltree; inside the container:
#   singularity exec -B /storeData:/storeData -B /mnt/hustor-04/zebra/ycai:/mnt/hustor-04/zebra/ycai \
#       $SIF_DIR/CWGS.sif bash cmrg_from_results.sh <outdir> <sample>
# result: <outdir>/<sample>/cmrg_manual/<sample>.cmrg.report  (+ the hist/mean beds, per gene)
set -euo pipefail

OUT=${1:?usage: $0 <outdir> <sample> [align_tool] [var_tool]}
ID=${2:?sample id}
ALIGN=${3:-lariat}
VAR=${4:-dv}
HERE=$(cd "$(dirname "$0")" && pwd)
DB=${DB:-/storeData/cg_teams/cg_research/ycai/completeWGS/pipeline/v1.0.6/CWGS_db}
REF=${REF:-${DB}/hg38/panGenome/GCA_000001405.15_GRCh38_no_alt_analysis_set_corrected.fasta}

# same python env as the report0 process, when present
if [ -f /usr/local/miniconda3/bin/activate ] && [ -d /usr/local/miniconda3/envs/six ]; then
    set +u; source /usr/local/miniconda3/bin/activate /usr/local/miniconda3/envs/six; set -u
fi

BAM=${OUT}/${ID}/align/${ID}.${ALIGN}.merge.bam
PVCF=${OUT}/${ID}/phase/${ID}.${ALIGN}.${VAR}.phased.vcf.gz
HB=${OUT}/${ID}/phase/${ID}.${ALIGN}.${VAR}.hapblock
GENES=${DB}/hg38/GRCh38_CMRG_benchmark_gene_coordinates.bed
EXON=${HERE}/cmrg273_exon.bed

missing=0
for f in "$BAM" "${BAM}.bai" "$PVCF" "$HB" "${REF}.fai" "$GENES" "$EXON"; do
    [ -s "$f" ] || { echo "MISSING: $f" >&2; missing=1; }
done
[ $missing -eq 0 ] || { echo "ls ${OUT}/${ID}/align ${OUT}/${ID}/phase to see the real file names; override with ALIGN/VAR args" >&2; exit 1; }

W=${OUT}/${ID}/cmrg_manual
mkdir -p "$W"; cd "$W"

# coverage fraction + mean depth per CMRG gene (same commands as the coverage / coverageMean processes)
bedtools coverage -sorted -g "${REF}.fai" -a "$GENES" -b "$BAM"        > "${ID}.cmrg.hist.bed"
bedtools coverage -sorted -g "${REF}.fai" -a "$GENES" -b "$BAM" -mean  > "${ID}.cmrg.mean.bed"

# coding variants in CMRG exons + phase-block coverage of genes (same inputs as report0)
ln -sf "$GENES" bed
ln -sf "$HB" hapblock
bedtools intersect -a "$PVCF" -b "$EXON" -wb > cmrg_exon.vcf

python3 "${HERE}/cmrg_from_results.py" "$ID" | tee "${ID}.cmrg.report"
echo "written: ${W}/${ID}.cmrg.report" >&2
