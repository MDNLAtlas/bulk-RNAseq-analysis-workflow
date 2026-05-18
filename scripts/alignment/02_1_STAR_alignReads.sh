#!/usr/bin/bash

# shellcheck disable=SC2090
# shellcheck disable=SC2086
# shellcheck disable=SC2089

set -euo pipefail

# Specify step
step_id="02_1"
step_n="STAR_alignReads"

# Specify project
n_jobs=2
PRJ_n="PRJ-UNITI/DATA_ANALYSIS_UNITI/BulkRNA_KO_PACT"
PRJ_DIR="$HOME/myRDS/$PRJ_n"
TMP_DIR="$HOME/Documents/Tam_Jenny/$PRJ_n"
REF="$HOME/myRDS/PRJ-MRefD/databases"

# Specific parameters for the function
threads=8
genomeDir="$REF/STAR-aligner-idx/v2.7.10b/hg38_primary"
RAW="$PRJ_DIR/raw"
ALIGNED="$PRJ_DIR/STAR-aligned"
TRIMMED="$PRJ_DIR/trimmed"

mkdir -p "$ALIGNED"
mkdir -p "$TMP_DIR"

########################################
# Function: FASTQC per SAMPLE + LANE
########################################

STAR_ALIGNREADS () {
  SAMPLE="$1"
  LANE="$2"
  LIB="$3"

  raw="$RAW/$SAMPLE/$LANE"
  trimmed="$TRIMMED/$SAMPLE/$LANE"
  aligned="$ALIGNED/$SAMPLE/$LANE"
  tmp_dir="$TMP_DIR/STAR-aligned/$SAMPLE/$LANE"
  LOG="$aligned/${step_id}_${step_n}.log"

  mkdir -p $aligned
  mkdir -p $trimmed
  mkdir -p $TMP_DIR/STAR-aligned/$SAMPLE # not create folder

  {
    echo "====== $(date): Processing $SAMPLE - $LANE ======"

    echo "Checking inputs exist ..."
    if [[ ! -d "$trimmed" ]]; then
        echo "[ERROR] Missing input directory"
        exit 1
    fi  

    echo "Resetting tmp dir if it exists ..."
    rm -rf "${tmp_dir:?}"

    echo "Resetting trimmed dir if it exists ..."
    rm -rf "${trimmed:?}"
    mkdir -p "$trimmed"

    if [[ "$LIB" == "PAIRED" ]]; then

      R1="$trimmed/${LANE}_R1.fastq.gz"
      src_R1="$raw/${LANE}_R1.fastq.gz"
      ln -s "$src_R1" "$R1"

      R2="$trimmed/${LANE}_R2.fastq.gz"
      src_R2="$raw/${LANE}_R2.fastq.gz"
      ln -s "$src_R2" "$R2"

      if [[ ! -f "$R1" || ! -f "$R2" ]]; then
          echo "[ERROR] Missing FASTQ files"
          exit 2
      fi
      echo "[DONE] All inputs are valid!"

      echo "Running $step_n ..."
      CMD="/usr/bin/time -v STAR \
        --runMode alignReads \
        --genomeDir $genomeDir \
        --outTmpDir $tmp_dir \
        --genomeLoad NoSharedMemory \
        --readFilesIn $R1 $R2 \
        --runThreadN $threads \
        --readFilesCommand zcat \
        --outFilterType BySJout \
        --alignSJoverhangMin 8 \
        --outFilterMultimapNmax 20 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMattributes NH HI AS NM MD \
        --quantMode TranscriptomeSAM \
        --outFileNamePrefix $aligned/${LANE}. \
        --outSAMtype BAM Unsorted"
      echo $CMD
      eval $CMD
      status=$?
      if [[ $status -ne 0 ]]; then
          echo "[ERROR] $step_n failed (exit code $status)"
          exit 3
      fi

      echo "Unlink tmp fastq from source raw folder ..."
      unlink "$R1"
      unlink "$R2"

    elif [[ "$LIB" == "SINGLE" ]]; then

      R1="$trimmed/${LANE}.fastq.gz"
      src_R1="$raw/${LANE}.fastq.gz"
      ln -s "$src_R1" "$R1"

      if [[ ! -f "$R1" ]]; then
          echo "[ERROR] Missing FASTQ file."
          exit 2
      fi
      echo "[DONE] All inputs are valid!"

      echo "Running $step_n ..."
      CMD="/usr/bin/time -v STAR \
        --runMode alignReads \
        --genomeDir $genomeDir \
        --outTmpDir $tmp_dir \
        --genomeLoad NoSharedMemory \
        --readFilesIn $R1 \
        --runThreadN $threads \
        --readFilesCommand zcat \
        --outFilterType BySJout \
        --alignSJoverhangMin 8 \
        --outFilterMultimapNmax 20 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMattributes NH HI AS NM MD \
        --quantMode TranscriptomeSAM \
        --outFileNamePrefix $aligned/${LANE}. \
        --outSAMtype BAM Unsorted"
      echo $CMD
      eval $CMD
      status=$?
      if [[ $status -ne 0 ]]; then
          echo "[ERROR] $step_n failed (exit code $status)"
          exit 3
      fi

      echo "Unlink tmp fastq from source raw folder ..."
      unlink "$R1"

    fi

    echo "[DONE] $step_n completed!"

    echo "====== $(date): All done! ======"

  } 2>&1 | tee "$LOG"

}

########################################
# export function
########################################
export RAW
export TRIMMED
export ALIGNED
export TMP_DIR
export genomeDir
export threads
export step_id
export step_n
export -f STAR_ALIGNREADS

########################################
# Define target samples/lanes for parallel processing
########################################
target="PRJNA1064892"
target_file="$PRJ_DIR/config/samples_lanes_${target}.txt"
parallel_log="$TMP_DIR/log/parallel/$target/${step_id}_${step_n}.log"

mkdir -p "$(dirname "$parallel_log")"

########################################
# Run parallel
########################################

parallel --joblog "$parallel_log" \
  --resume \
  -j "$n_jobs" \
  --colsep '\t' \
  "conda run -n STAR-2.7.10b STAR_ALIGNREADS {1} {2} {3}" \
  :::: "$target_file"

########################################
# Summary log file
########################################
LOG="$PRJ_DIR/log/completion_status/$target"
mkdir -p "$LOG"

cat "$target_file"| cut -f1,2| while read -r SAMPLE LANE; do
  echo "====== $SAMPLE - $LANE ======"
  grep -B40 "Exit status:" "$ALIGNED/$SAMPLE/$LANE/${step_id}_${step_n}.log" >&2
  grep "Exit status:" "$ALIGNED/$SAMPLE/$LANE/${step_id}_${step_n}.log"
done > "$LOG/${step_id}_${step_n}.log"
