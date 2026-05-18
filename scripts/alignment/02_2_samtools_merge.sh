#!/usr/bin/bash

# shellcheck disable=SC2090
# shellcheck disable=SC2086
# shellcheck disable=SC2089
# shellcheck disable=SC2126
# shellcheck disable=SC2030

set -euo pipefail

# Specify step
step_id="02_2"
step_n="samtools_merge"

# Specify project
n_jobs=2
PRJ_n="PRJ-UNITI/DATA_ANALYSIS_UNITI/BulkRNA_KO_PACT"
PRJ_DIR="$HOME/myRDS/$PRJ_n"
TMP_DIR="$HOME/Documents/Tam_Jenny/$PRJ_n"

# Specific parameters for the function
threads=8
ALIGNED="$PRJ_DIR/STAR-aligned"

mkdir -p "$ALIGNED"
mkdir -p "$TMP_DIR"

########################################
# Function: FASTQC per SAMPLE + LANE
########################################

SAMTOOLS_MERGE () {
  SAMPLE="$1"
  n_lanes="$2"

  aligned="$ALIGNED/$SAMPLE"
  LOG="$aligned/${step_id}_${step_n}.log"

  {
    echo "====== $(date): Processing $SAMPLE - $LANE ======"

    echo "Checking inputs exist ..."

    if [[ $n_lanes -gt 1 ]]; then

      bam_list=""

      for LANE in "$aligned"/*; do
        [[ -d "$LANE" ]] || continue
        LANE_NAME=$(basename "$LANE")
        bam="$LANE/${LANE_NAME}.Aligned.toTranscriptome.out.bam"

        if [[ ! -f "$bam" ]]; then
          echo "[ERROR] Missing BAM file $bam."
          exit 2
        fi

        bam_list+=" $bam"
      done

      echo "Resetting tmp dir if it exists ..."
      rm -rf "${tmp_dir:?}"
        echo "[DONE] All inputs are valid!"

      echo "Running $step_n ..."

      merged_bam="$aligned/${SAMPLE}.Aligned.toTranscriptome.out.bam"

      CMD="/usr/bin/time -v \
        samtools merge \
        --threads $threads \
        $merged_bam \
        $bam_list"
      echo $CMD
      eval $CMD
      status=$?
      if [[ $status -ne 0 ]]; then
          echo "[ERROR] $step_n failed (exit code $status)"
          exit 3
      fi

    else

      echo "Sample $SAMPLE has only 1 lane!"

    fi

    echo "[DONE] $step_n completed!"

    echo "====== $(date): All done! ======"

  } 2>&1 | tee "$LOG"

}

########################################
# export function
########################################
export ALIGNED
export threads
export step_id
export step_n
export -f SAMTOOLS_MERGE

########################################
# Define target samples/lanes for parallel processing
########################################
target="PRJNA1064892"
target_file=$(cut -f1 "$PRJ_DIR/config/samples_lanes_${target}.txt"| \
  sort| \
  uniq -c| \
  awk -F' ' '{OFS="\t"}{print $2, $1}')
parallel_log="$TMP_DIR/log/parallel/$target/${step_id}_${step_n}.log"

mkdir -p "$(dirname "$parallel_log")"

########################################
# Run parallel
########################################

parallel --joblog "$parallel_log" \
  --resume \
  -j "$n_jobs" \
  --colsep '\t' \
  "conda run -n bcftools-1.23.1 SAMTOOLS_MERGE {1} {2}" \
  :::: "$target_file"

########################################
# Summary log file
########################################
LOG="$PRJ_DIR/log/completion_status/$target"
mkdir -p "$LOG"

echo "$target_file"| cut -f1| while read -r SAMPLE; do
  echo "====== $SAMPLE - $LANE ======"
  grep "Exit status:" "$RAW/$SAMPLE/${step_id}_${step_n}.log"
done > "$LOG/${step_id}_${step_n}.log"
