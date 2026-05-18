#!/bin/bash

# shellcheck disable=SC2086
# shellcheck disable=SC2089
# shellcheck disable=SC2090

#############################################################################################################################
### This script performs read trimming using trim-galore for paired-end sequencing data.
### It takes the raw fastq files as input and outputs the trimmed fastq files in the specified output directory.
### The script also generates a log file for each sample, which contains the commands executed and the timestamps for each step of the trimming process.
### The script uses GNU parallel to process multiple samples in parallel, which can significantly reduce the overall runtime of the trimming step.
### The script assumes that the raw fastq files are organized in a specific directory structure, where each sample has its own subdirectory under the raw data directory, 
### and the fastq files are named in a specific format
#############################################################################################################################

# Specify step
step_id="01_2"
step_n="trim_galore"

# Specify project
PRJ_n="PRJ-UNITI/DATA_ANALYSIS_UNITI/BulkRNA_KO_PACT"
PRJ_DIR="$HOME/myRDS/$PRJ_n"
TMP_DIR="$HOME/Documents/Tam_Jenny/$PRJ_n"
RAW="$PRJ_DIR/raw"

# Specific parameters for the function
TRIMMED="$PRJ_DIR/trimmed"
n_cores=4
threads=8
min_len=50
min_qual=20

# Parameters for parallel
n_jobs=2

mkdir -p "$PRJ_DIR"
mkdir -p "$RAW" 
mkdir -p "$TRIMMED"
mkdir -p "$TMP_DIR"

########################################
# Function: TRIMGALORE per SAMPLE + LANE
########################################

TRIMGALORE () {
  SAMPLE="$1"
  LANE="$2"
  LIB="$3"

  raw="$RAW/$SAMPLE/$LANE"
  trimmed="$TRIMMED/$SAMPLE/$LANE"
  tmp_dir="$TMP_DIR/trimmed/$SAMPLE/$LANE"
  LOG="$trimmed/${step_id}_${step_n}.log"

  mkdir -p "$tmp_dir" "$OUT_DIR"

  {
    echo "====== $(date): Processing $SAMPLE - $LANE ======"

    echo "Checking inputs exist ..."
    if [[ ! -d "$raw" ]]; then
        echo "[ERROR] Missing input directory"
        exit 1
    fi

    if [[ "$LIB" == "PAIRED" ]]; then

      R1="$raw/${LANE}_R1.fastq.gz"
      R2="$raw/${LANE}_R2.fastq.gz"
      if [[ ! -f "$R1" || ! -f "$R2" ]]; then
          echo "[ERROR] Missing FASTQ files"
          exit 2
      fi
      echo "Resetting tmp dir if it exists ..."
      rm -rf "${tmp_dir:?}"
      mkdir -p "$tmp_dir"
      echo "[DONE] All inputs are valid!"

      echo "Running $step_n ..."
      CMD="/usr/bin/time -v trim_galore --cores $n_cores \
          -q $min_qual \
          --length $min_len \
          --paired \
          --fastqc_args \"--threads $threads --dir $tmp_dir\" \
          --output_dir $trimmed \
          $R1 $R2"
      echo $CMD
      eval $CMD
      status=$?
      if [[ $status -ne 0 ]]; then
          echo "[ERROR] $step_n failed (exit code $status)"
          exit 3
      fi
      echo "[DONE] $step_n completed!"

      echo "Standardizing new fastq names ..."
      R1_val="$trimmed/${LANE}_R1_val_1.fq.gz"
      R2_val="$trimmed/${LANE}_R2_val_2.fq.gz"

      R1_renamed="$trimmed/${LANE}_R1.fastq.gz"
      R2_renamed="$trimmed/${LANE}_R2.fastq.gz"

      [[ -f "$R1_val" ]] || { echo "[ERROR] Missing $R1_val"; exit 4; }
      [[ -f "$R2_val" ]] || { echo "[ERROR] Missing $R2_val"; exit 4; }

      CMD="mv $R1_val $R1_renamed; \
          mv $R2_val $R2_renamed"
      echo $CMD
      eval $CMD
      status=$? 
      if [[ $status -ne 0 ]]; then
          echo "[ERROR] Rename failed (exit code $status)"
          exit 4
      fi
      echo "[DONE] Rename completed!"

    elif [[ "$LIB" == "SINGLE" ]]; then

      R1="$raw/${LANE}.fastq.gz"

      if [[ ! -f "$R1" ]]; then
          echo "[ERROR] Missing FASTQ files"
          exit 2
      fi
      echo "Resetting tmp dir if it exists ..."
      rm -rf "${tmp_dir:?}"
      mkdir -p "$tmp_dir"
      echo "[DONE] All inputs are valid!"

      echo "Running $step_n ..."
      CMD="/usr/bin/time -v trim_galore --cores $n_cores \
          -q $min_qual \
          --length $min_len \
          --fastqc_args \"--threads $threads --dir $tmp_dir\" \
          --output_dir $trimmed \
          $R1"
      echo $CMD
      eval $CMD
      status=$?
      if [[ $status -ne 0 ]]; then
          echo "[ERROR] $step_n failed (exit code $status)"
          exit 3
      fi
      echo "[DONE] $step_n completed!"

      echo "Standardizing new fastq names ..."
      R1_val="$trimmed/${LANE}_trimmed.fq.gz"

      R1_renamed="$trimmed/${LANE}.fastq.gz"

      [[ -f "$R1_val" ]] || { echo "[ERROR] Missing $R1_val"; exit 4; }

      CMD="mv $R1_val $R1_renamed"
      echo $CMD
      eval $CMD
      status=$? 
      if [[ $status -ne 0 ]]; then
          echo "[ERROR] Rename failed (exit code $status)"
          exit 4
      fi
      echo "[DONE] Rename completed!"

    fi

    echo "====== $(date): All done! ======"

  } 2>&1 | tee "$LOG"

}

########################################
# export variables and functions globally for parallel
########################################
export RAW
export TRIMMED
export TMP_DIR
export threads
export n_cores
export min_len
export min_qual
export step_id
export step_n
export -f TRIMGALORE

########################################
# Define target samples/lanes for parallel processing
########################################
target="batch_1"
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
  "conda run -n trim_galore-0.6.11 TRIMGALORE {1} {2} {3}" \
  :::: "$target_file"

########################################
# Summary log file
########################################
LOG="$PRJ_DIR/log/completion_status/$target"
mkdir -p "$LOG"

while read -r SAMPLE LANE; do
  echo "====== $SAMPLE - $LANE ======"
  grep "Exit status:" "$TRIMMED/$SAMPLE/$LANE/${step_id}_${step_n}.log"
done < "$target_file" > "$LOG/${step_id}_${step_n}.log"