#!/usr/bin/bash

# shellcheck disable=SC2090
# shellcheck disable=SC2086
# shellcheck disable=SC2089

set -euo pipefail

# Specify step
step_id="03_1"
step_n="rsem_calculate_expression"

# Specify project
n_jobs=2
PRJ_n="PRJ-UNITI/DATA_ANALYSIS_UNITI/BulkRNA_KO_PACT"
PRJ_DIR="$HOME/myRDS/$PRJ_n"
TMP_DIR="$HOME/Documents/Tam_Jenny/$PRJ_n"
REF="$HOME/myRDS/PRJ-MRefD/databases"

# Specific parameters for the function
threads=8
genomeDir="$REF/RSEM-idx/v1.3.3/hg38_primary/hg38_primary"
ALIGNED="$PRJ_DIR/STAR-aligned"
QUANT="$PRJ_DIR/quantification/RSEM"
# Refer to https://academic.oup.com/bioinformatics/article/38/6/1491/6493231
# rl=100
FRAGMENT_LEN_MEAN=300
FRAGMENT_LEN_SD=50

mkdir -p "$QUANT"
mkdir -p "$TMP_DIR"

########################################
# Function: FASTQC per SAMPLE + LANE
########################################

RSEM_CALCULATE_EXPRESSION () {
  SAMPLE="$1"
  LIB=$(echo "$2"| cut -d";" -f1)
  STRAND=$(echo "$2"| cut -d";" -f2)
  n_lanes="$3"

  aligned="$ALIGNED/$SAMPLE"
  quant="$QUANT/$SAMPLE"
  tmp_dir="$TMP_DIR/quanitification/RSEM/$SAMPLE"
  LOG="$quant/${step_id}_${step_n}.log"

  mkdir -p $quant
  mkdir -p $tmp_dir

  {
    echo "====== $(date): Processing $SAMPLE - $LANE ======"

    echo "Checking inputs exist ..."

    if [[ $n_lanes -gt 1 ]]; then

      bam="$aligned/${SAMPLE}.Aligned.toTranscriptome.out.bam"

    else

      LANE=$(ls $aligned)
      bam="$aligned/$LANE/${LANE}.Aligned.toTranscriptome.out.bam"

    fi

    if [[ ! -f "$bam" ]]; then
        echo "[ERROR] Missing BAM file."
        exit 2
    fi

    echo "Resetting tmp dir if it exists ..."
    rm -rf "${tmp_dir:?}"
    mkdir -p "$tmp_dir"

    echo "[DONE] All inputs are valid!"

    echo "Running $step_n ..."

    if [[ "$LIB" == "PAIRED" ]]; then

      if [[ "$STRAND" == "stranded" ]]; then

        params="--paired-end --forward-prob 0"

      elif [[ "$STRAND" == "non-stranded" ]]; then

        params="--paired-end --forward-prob 0.5"

      fi

    elif [[ "$LIB" == "SINGLE" ]]; then

      params="--fragment-length-mean $FRAGMENT_LEN_MEAN \
        --fragment-length-sd $FRAGMENT_LEN_SD"

    fi

    CMD="/usr/bin/time -v \
      rsem-calculate-expression \
      --bam \
      --no-bam-output \
      -p $threads \
      $params \
      $bam \
      $genomeDir \
      $quant/${SAMPLE} \
      --temporary-folder $tmp_dir"
    echo $CMD
    eval $CMD
    status=$?
    if [[ $status -ne 0 ]]; then
        echo "[ERROR] $step_n failed (exit code $status)"
        exit 3
    fi

    echo "[DONE] $step_n completed!"

    echo "====== $(date): All done! ======"

  } 2>&1 | tee "$LOG"

}

########################################
# export function
########################################
export ALIGNED
export QUANT
export genomeDir
export TMP_DIR
export STRAND
export FRAGMENT_LEN_MEAN
export FRAGMENT_LEN_SD
export threads
export step_id
export step_n
export -f RSEM_CALCULATE_EXPRESSION

########################################
# Define target samples/lanes for parallel processing
########################################
target="PRJNA1064892"
target_file=$(cut -f1,3 "$PRJ_DIR/config/samples_lanes_${target}.txt"| \
  sort| \
  uniq -c| \
  awk '{OFS="\t"} {print $2, $3, $1}')
parallel_log="$TMP_DIR/log/parallel/$target/${step_id}_${step_n}.log"

mkdir -p "$(dirname "$parallel_log")"

########################################
# Run parallel
########################################

parallel --joblog "$parallel_log" \
  --resume \
  -j "$n_jobs" \
  --colsep '\t' \
  "conda run -n rsem-1.3.3 RSEM_CALCULATE_EXPRESSION {1} {2} {3}" \
  ::: "$target_file"

########################################
# Summary log file
########################################
LOG="$PRJ_DIR/log/completion_status/$target"
mkdir -p "$LOG"

echo "$target_file"| cut -f1| while read -r SAMPLE; do
  echo "====== $SAMPLE ======"
  grep -B40 "Exit status:" "$QUANT/$SAMPLE/${step_id}_${step_n}.log" >&2
  grep "Exit status:" "$QUANT/$SAMPLE/${step_id}_${step_n}.log"
done > "$LOG/${step_id}_${step_n}.log"

##################################################################################
# # Specify step
# step_id="03_1"
# step_n="rsem_calculate_expression"

# # Specify project
# n_jobs=2
# PRJ_n="PRJ-UNITI/DATA_ANALYSIS_UNITI/BulkRNA_KO_PACT"
# PRJ_DIR="$HOME/myRDS/$PRJ_n"
# TMP_DIR="$HOME/Documents/Tam_Jenny/$PRJ_n"
# REF="$HOME/myRDS/PRJ-MRefD/databases"

# # Specific parameters for the function
# threads=8
# genomeDir="$REF/RSEM-idx/v1.3.3/hg38_primary/hg38_primary"
# ALIGNED="$PRJ_DIR/STAR-aligned"
# QUANT="$PRJ_DIR/quantification/RSEM"
# # Refer to https://academic.oup.com/bioinformatics/article/38/6/1491/6493231
# # rl=100
# FRAGMENT_LEN_MEAN=300
# FRAGMENT_LEN_SD=50

# mkdir -p "$QUANT"
# mkdir -p "$TMP_DIR"

# ########################################
# # Function: FASTQC per SAMPLE + LANE
# ########################################

# RSEM_CALCULATE_EXPRESSION () {
#   SAMPLE="$1"
#   LANE="$2"
#   LIB=$(echo "$3"| cut -d";" -f1)
#   STRAND=$(echo "$3"| cut -d";" -f2)

#   aligned="$ALIGNED/$SAMPLE/$LANE"
#   quant="$QUANT/$SAMPLE/$LANE"
#   tmp_dir="$TMP_DIR/quanitification/RSEM/$SAMPLE/$LANE"
#   LOG="$quant/${step_id}_${step_n}.log"

#   mkdir -p $quant
#   mkdir -p $tmp_dir

#   {
#     echo "====== $(date): Processing $SAMPLE - $LANE ======"

#     echo "Checking inputs exist ..."

#     bam="$aligned/${SAMPLE}.Aligned.toTranscriptome.out.bam"
#     if [[ ! -f "$bam" ]]; then
#         echo "[ERROR] Missing BAM file."
#         exit 2
#     fi

#     echo "Resetting tmp dir if it exists ..."
#     rm -rf "${tmp_dir:?}"
#     mkdir -p "$tmp_dir"

#     echo "[DONE] All inputs are valid!"

#     echo "Running $step_n ..."

#     if [[ "$LIB" == "PAIRED" ]]; then

#       if [[ "$STRAND" == "stranded" ]]; then

#         params="--paired-end --forward-prob 0"

#       elif [[ "$STRAND" == "non-stranded" ]]; then

#         params="--paired-end --forward-prob 0.5"

#       fi

#     elif [[ "$LIB" == "SINGLE" ]]; then

#       params="--fragment-length-mean $FRAGMENT_LEN_MEAN \
#         --fragment-length-sd $FRAGMENT_LEN_SD"

#     fi

#     CMD="/usr/bin/time -v \
#       rsem-calculate-expression \
#       --bam \
#       --no-bam-output \
#       -p $threads \
#       $params \
#       $bam \
#       $genomeDir \
#       $quant/${SAMPLE} \
#       --temporary-folder $tmp_dir"
#     echo $CMD
#     eval $CMD
#     status=$?
#     if [[ $status -ne 0 ]]; then
#         echo "[ERROR] FastQC failed (exit code $status)"
#         exit 3
#     fi

#     echo "[DONE] $step_n completed!"

#     echo "====== $(date): All done! ======"

#   } 2>&1 | tee "$LOG"

# }

# ########################################
# # export function
# ########################################
# export ALIGNED
# export QUANT
# export genomeDir
# export TMP_DIR
# export STRAND
# export FRAGMENT_LEN_MEAN
# export FRAGMENT_LEN_SD
# export threads
# export step_id
# export step_n
# export -f RSEM_CALCULATE_EXPRESSION

# ########################################
# # Define target samples/lanes for parallel processing
# ########################################
# target="PRJNA1064892"
# target_file="$PRJ_DIR/config/samples_lanes_${target}.txt"
# parallel_log="$TMP_DIR/log/parallel/$target/${step_id}_${step_n}.log"

# mkdir -p "$(dirname "$parallel_log")"

# ########################################
# # Run parallel
# ########################################

# parallel --joblog "$parallel_log" \
#   --resume \
#   -j "$n_jobs" \
#   --colsep '\t' \
#   "conda run -n rsem-1.3.3 RSEM_CALCULATE_EXPRESSION {1} {2} {3}" \
#   :::: "$target_file"

# ########################################
# # Summary log file
# ########################################
# LOG="$PRJ_DIR/log/completion_status/$target"
# mkdir -p "$LOG"

# while read -r SAMPLE LANE; do
#   echo "====== $SAMPLE - $LANE ======"
#   grep "Exit status:" "$RAW/$SAMPLE/$LANE/${step_id}_${step_n}.log"
# done < "$target_file" > "$LOG/${step_id}_${step_n}.log"
