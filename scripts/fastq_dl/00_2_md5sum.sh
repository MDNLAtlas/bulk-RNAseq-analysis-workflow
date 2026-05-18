#!/bin/bash

### This script computes MD5 checksums for the raw FASTQ files and compares them to the original checksums provided by the sequencing facility. 

# shellcheck disable=SC2086
# shellcheck disable=SC2128
# shellcheck disable=SC2178
# shellcheck disable=SC2329

set -euo pipefail

step_id="00_2"
step_n="md5sum"

n_jobs=2
PRJ_n="PRJ-UNITI/DATA_ANALYSIS_UNITI/BulkRNA_KO_PACT"
PRJ_DIR="$HOME/myRDS/$PRJ_n"
TMP_DIR="$HOME/Documents/Tam_Jenny/$PRJ_n"
RAW="$PRJ_DIR/raw"

mkdir -p "$PRJ_DIR"
mkdir -p "$RAW"
mkdir -p "$TMP_DIR"

########################################
# Function: MD5 check per SAMPLE + LANE
########################################

MD5SUM () {
  SAMPLE="$1"
  LANE="$2"
  LIB="$3"

  DIR="$RAW/$SAMPLE/$LANE"
  LOG="$DIR/${step_id}_${step_n}.log"

  {
    echo "====== $(date): Processing $SAMPLE - $LANE ======"

    if [[ $LIB == "PAIRED" ]]; then

      R1="$DIR/${LANE}_R1.fastq.gz"
      R2="$DIR/${LANE}_R2.fastq.gz"

      # Check files exist
      echo "Checking fastq files exist ..."
      [[ -f "$R1" ]] || { echo "[ERROR] Missing $R1"; exit 2; }
      [[ -f "$R2" ]] || { echo "[ERROR] Missing $R2"; exit 2; }
      echo "[DONE] Checking completed!"

      # Generate md5
      echo "Generating md5 ..."
      CMD=(

        "/usr/bin/time -v md5sum -b $R1 > ${R1}.DC.md5"

        "/usr/bin/time -v md5sum -b $R2 > ${R2}.DC.md5"
        
        )
      for cmd in "${CMD[@]}"; do
        echo $cmd
        eval $cmd
      done
      echo "[DONE] Generating md5 completed!"

      echo "Comparing md5 ..."
      # Compare R1 md5
      if [[ -f "${R1}.md5" ]]; then
        md=$(cut -d' ' -f1 "${R1}.md5")
      else
        md="NA"
      fi
      md_check=$(cut -d' ' -f1 "${R1}.DC.md5")
      echo "R1 Original MD5: $md"
      echo "R1 DC MD5: $md_check"
      [[ "$md" == "$md_check" ]] && echo "R1 OK" || echo "R1 Not OK"

      # Compare R2 md5
      if [[ -f "${R2}.md5" ]]; then
        md=$(cut -d' ' -f1 "${R2}.md5")
      else
        md="NA"
      fi
      md_check=$(cut -d' ' -f1 "${R2}.DC.md5")
      echo "R2 Original MD5: $md"
      echo "R2 DC MD5: $md_check"
      [[ "$md" == "$md_check" ]] && echo "R2 OK" || echo "R2 Not OK"
      echo "[DONE] Comparing md5 completed!"

    elif [[ $LIB == "SINGLE" ]]; then

      R1="$DIR/${LANE}.fastq.gz"

      # Check files exist
      echo "Checking fastq files exist ..."
      [[ -f "$R1" ]] || { echo "[ERROR] Missing $R1"; exit 2; }
      echo "[DONE] Checking completed!"

      # Generate md5
      echo "Generating md5 ..."
      
      CMD="/usr/bin/time -v md5sum -b $R1 > ${R1}.DC.md5"
      echo $CMD
      eval $CMD
      echo "[DONE] Generating md5 completed!"

      echo "Comparing md5 ..."
      # Compare R1 md5
      if [[ -f "${R1}.md5" ]]; then
        md=$(cut -d' ' -f1 "${R1}.md5")
      else
        md="NA"
      fi
      md_check=$(cut -d' ' -f1 "${R1}.DC.md5")
      echo "R1 Original MD5: $md"
      echo "R1 DC MD5: $md_check"
      [[ "$md" == "$md_check" ]] && echo "R1 OK" || echo "R1 Not OK"

      echo "[DONE] Comparing md5 completed!"

    fi

    echo "====== $(date): All done! ======"

  } 2>&1 | tee "$LOG"
}

########################################
# export variables and functions globally for parallel
########################################
export RAW
export -f MD5SUM

########################################
# Define target samples/lanes for parallel processing
########################################
target="PRJNA1064892"
target_file="$PRJ_DIR/config/samples_lanes_${target}.txt"
parallel_log="$TMP_DIR/log/parallel/$target/${step_id}_${step_n}.log"

mkdir -p "$(dirname $parallel_log)"

########################################
# Run parallel
########################################
parallel --joblog "$parallel_log" \
  --resume \
  -j "$n_jobs" \
  --colsep '\t' \
  "MD5SUM {1} {2} {3}" \
  :::: "$target_file"
