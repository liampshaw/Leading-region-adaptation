#!/bin/bash
# Script to run DefenseFinder (with anti-defense genes) on computing cluster
# Provided for information only: identifying information required to run (username, account) has been replaced with XXXX
# Liam Shaw, April 2026

#SBATCH --job-name=defensefinder
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:15:00 
#SBATCH --mem=200M
#SBATCH --account=XXXX
#SBATCH --array=1-1751%50  # -1751%50 # -1751%50

set -euo pipefail
# To reduce launch failure 
sleep $((RANDOM % 15))


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" prodigal_samples.txt)
#name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rerun_2.txt)

hostname
# Activate conda environment
source activate defensefinder-2.0.1

OUT_TMP=/tmp/defense_finder_tmp_${SLURM_ARRAY_TASK_ID}
OUT_FINAL=/user/home/XXXX/trieste/prodigal-meta-defense-finder-outputs
mkdir -p $OUT_TMP
# copy to tmp folder? cp "$name" 
# Check that file exists
if [ -z "prodigal/$name" ]; then
    echo "No input file for array index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [ ! -f "prodigal/$name" ]; then
    echo "File prodigal/$name does not exist"
    exit 1
fi


cp prodigal/"$name" $OUT_TMP
cd $OUT_TMP
defense-finder run -a "$name" --workers ${SLURM_CPUS_PER_TASK} -o $OUT_TMP
#rsync -av --remove-source-files "${OUT_TMP}/*.tsv" "${OUT_FINAL}/"
# Copy the annotations and the defensefinder systems
cp ${OUT_TMP}/* "${OUT_FINAL}/"
rm -r "${OUT_TMP}"
