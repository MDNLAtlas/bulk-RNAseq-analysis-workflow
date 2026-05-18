#!/bin/bash

#############################################################################################################################
### This script performs quality control using FastQC for paired-end sequencing data.
### It takes the raw fastq files as input and outputs the FastQC reports in the same directory as the raw fastq files.
### The script also generates a log file for each sample, which contains the commands executed and the timestamps for each step of the quality control process.
### The script uses GNU parallel to process multiple samples in parallel, which can significantly reduce the overall runtime of the quality control step.
### The script assumes that the raw fastq files are organized in a specific directory structure, where each sample has its own subdirectory under the raw data directory, 
### and the fastq files are named in a specific format
#############################################################################################################################

# shellcheck disable=SC2090
# shellcheck disable=SC2086
# shellcheck disable=SC2089

set -euo pipefail

# Specify step
step_id="01_1"
step_n="fastqc"

# Specify project
n_jobs=2
PRJ_n="PRJ-UNITI/DATA_ANALYSIS_UNITI/BulkRNA_KO_PACT"
PRJ_DIR="$HOME/myRDS/$PRJ_n"
TMP_DIR="$HOME/Documents/Tam_Jenny/$PRJ_n"
RAW="$PRJ_DIR/raw"

# Specific parameters for the function
threads=8

mkdir -p "$PRJ_DIR"
mkdir -p "$RAW"
mkdir -p "$TMP_DIR"

########################################
# Function: FASTQC per SAMPLE + LANE
########################################

FASTQC () {
  SAMPLE="$1"
  LANE="$2"
  LIB="$3"

  raw="$RAW/$SAMPLE/$LANE"
  LOG="$DIR/${step_id}_${step_n}.log"

  tmp_dir="$TMP_DIR/raw/$SAMPLE/$LANE"
  mkdir -p "$tmp_dir"

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
      echo "[DONE] All inputs are valid!"

      echo "Running $step_n ..."
      CMD="/usr/bin/time -v fastqc -t $threads --dir $tmp_dir -o $(dirname $R1) $R1 $R2"
      echo $CMD
      eval $CMD
      status=$? # capture exit status
      # Check fastqc result
      if [[ $status -ne 0 ]]; then
          echo "[ERROR] $step_n failed (exit code $status)"
          exit 3
      fi

      echo "Moving html files to TMP_DIR to view ..."
      cp $raw/*.html $tmp_dir

    elif [[ "$LIB" == "SINGLE" ]]; then

      R1="$raw/${LANE}.fastq.gz"

      if [[ ! -f "$R1" ]]; then
          echo "[ERROR] Missing FASTQ files"
          exit 2
      fi
      echo "[DONE] All inputs are valid!"

      echo "Running $step_n ..."
      CMD="/usr/bin/time -v fastqc -t $threads --dir $tmp_dir -o $(dirname $R1) $R1"
      echo $CMD
      eval $CMD
      status=$? # capture exit status
      # Check fastqc result
      if [[ $status -ne 0 ]]; then
          echo "[ERROR] $step_n failed (exit code $status)"
          exit 3
      fi

      echo "Moving html files to TMP_DIR to view ..."
      cp $raw/*.html $tmp_dir

    fi
    echo "[DONE] $step_n completed!"

    echo "====== $(date): All done! ======"

  } 2>&1 | tee "$LOG"

}

########################################
# export function
########################################
export RAW
export TMP_DIR
export threads
export step_id
export step_n
export -f FASTQC

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
  "conda run -n trim_galore-0.6.11 FASTQC {1} {2} {3}" \
  :::: "$target_file"

########################################
# Summary log file
########################################
LOG="$PRJ_DIR/log/completion_status/$target"
mkdir -p "$LOG"

while read -r SAMPLE LANE; do
  echo "====== $SAMPLE - $LANE ======"
  grep "Exit status:" "$RAW/$SAMPLE/$LANE/${step_id}_${step_n}.log"
done < "$target_file" > "$LOG/${step_id}_${step_n}.log"