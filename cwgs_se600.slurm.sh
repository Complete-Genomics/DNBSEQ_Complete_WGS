#!/usr/bin/env bash
#SBATCH --job-name=CWGS_SE600
#SBATCH --partition=dev
#SBATCH --output=cwgs_se600_%j.out
#SBATCH --error=cwgs_se600_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=320G
#SBATCH --time=100:00:00

# Run cWGS SE600 plus PCR-free FASTQ on one Slurm node.  The SE600 input must
# be listed in the stlfr21 column and already have #bc1_bc2_bc3/[12] in its FASTQ
# header; its sequence must contain gDNA only (the 42 bp barcode is removed).
#
# Submit from the directory containing sample.list:
#   sbatch /mnt/rw/hustor-04/zebra/ycai/projects/completeWGS/pipeline/se600/DNBSEQ_Complete_WGS/cwgs_se600.slurm.sh
#
# sample.list must contain sample, stlfr21, pcrfree1, and pcrfree2 columns.
# REF is GRCh38 because both SE600 and PCR-free use vg giraffe.

set -euo pipefail

src_path=/mnt/rw/hustor-04/zebra/ycai/projects/completeWGS/pipeline/se600
cwgs_path=/mnt/rw/zdspro-01/cg_teams/cg_research/ycai

modules=${src_path}/DNBSEQ_Complete_WGS/modules
scripts=${src_path}/DNBSEQ_Complete_WGS/scripts
ref=${cwgs_path}/completeWGS/pipeline/v1.0.6/CWGS_db/hg38/panGenome/GCA_000001405.15_GRCh38_no_alt_analysis_set_corrected.fasta
ref_len=2933974420

sif_dir=${cwgs_path}/completeWGS/pipeline/sif
db=${cwgs_path}/completeWGS/pipeline/v1.0.6/CWGS_db

# nextflow_env=${cwgs_path}/tools/miniconda3/envs/nextflow-25.04.6
export PATH=/home-03/ycai/.mamba/bin:$PATH

CWGS=${src_path}/DNBSEQ_Complete_WGS/CWGS
sample_list=${PWD}/sample.list
outdir=${PWD}/result
run_dir=${PWD}/CWGS_run
binds=/mnt/rw:/mnt/rw,/mnt:/mnt,/mnt/rw/qnap-archive:/mnt/rw/qnap-archive
cpu3=${SLURM_CPUS_PER_TASK:-24}

for path in "$sample_list" "$sif_dir" "$db" "$modules" "$scripts" "$CWGS"; do
    [[ -e "$path" ]] || { echo "Missing required path: $path" >&2; exit 2; }
done

[[ "$ref_len" =~ ^[0-9]+$ ]] || { echo "ref_len must be an integer: $ref_len" >&2; exit 2; }
command -v nextflow >/dev/null || {
    echo "nextflow not found in PATH: $PATH" >&2
    exit 2
}
mkdir -p "$outdir" "$run_dir"

java -version
nextflow -version

"$CWGS" run "$sample_list" \
    -sifs "$sif_dir" \
    -B "$binds" \
    -db "$db" \
    -exec local \
    -module "$modules" \
    -script "$scripts" \
    -run "$run_dir" \
    -out "$outdir" \
    --use_megabolt false \
    --frombam false \
    --skipBarcodeSplit true \
    --pfmapq 3 \
    --sampleBam true \
    --stLFR_sampling_cov 30 \
    --PF_sampling_cov 40 \
    --pfAligner vg \
    --ref "$ref" \
    --ref_len "$ref_len" \
    --cpu2 4 \
    --cpu3 "$cpu3" \
    2>&1 | tee "$outdir/cwgs_se600_${SLURM_JOB_ID:-local}.log"
