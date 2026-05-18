#!/usr/bin/bash


#!/usr/bin/bash

# shellcheck disable=SC2090
# shellcheck disable=SC2086
# shellcheck disable=SC2089

set -euo pipefail

# Specify step
step_id="00"
step_n="rsem_prepare_reference"

# Specify project
REF="$HOME/myRDS/PRJ-MRefD/databases"

# Specific parameters for the function
threads=8
genomeDir="$REF/RSEM-idx/v1.3.3/hg38_primary/hg38_primary"
FASTA="$REF/Gencode/v48/GRCh38.primary_assembly.genome.fa"
GTF="$REF/Gencode/v48/gencode.v48.primary_assembly.annotation.gtf"
LOG="$genomeDir/${step_id}_${step_n}.log"

# mkdir -p "$TMP_DIR"
echo "Resetting genomeDir if it exists ..."
rm -rf "${genomeDir:?}"
mkdir -p "$genomeDir"

########################################
# Function: BWAMEM2_MEM per SAMPLE + LANE
########################################

echo "====== $(date): Processing ======" > "$LOG"

echo "Checking inputs exist ..." >> "$LOG"

if [[ ! -f "$FASTA" ]]; then
    echo "[ERROR] Missing FASTQ file." >> "$LOG"
    exit 2
fi

if [[ ! -f "$GTF" ]]; then
    echo "[ERROR] Missing GTF file." >> "$LOG"
    exit 2
fi

echo "[DONE] All inputs are valid!" >> "$LOG"

echo "Running $step_n ..." >> "$LOG"
CMD="/usr/bin/time -v \
  conda run -n rsem-1.3.3 \
  rsem-prepare-reference \
  -p $threads \
  --gtf $GTF \
  $FASTA \
  $genomeDir"

echo $CMD >> "$LOG"
eval $CMD >> "$LOG" 2>&1
status=$? # capture exit status
if [[ $status -ne 0 ]]; then
    echo "[ERROR] $step_n failed (exit code $status)" >> "$LOG"
    exit 3
fi
echo "[DONE] $step_n completed!" >> "$LOG"

echo "====== $(date): All done! ======" >> "$LOG"