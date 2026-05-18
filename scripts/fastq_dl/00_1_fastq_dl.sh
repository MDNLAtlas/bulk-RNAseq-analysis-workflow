#!/usr/bin/bash

# shellcheck disable=SC2164
# shellcheck disable=SC2086

# Specify step
step_id="00"
step_n="download_fastq"

# Specify project
n_jobs=1
PRJ_n="PRJ-UNITI/DATA_ANALYSIS_UNITI/BulkRNA_KO_PACT"
PRJ_DIR="$HOME/myRDS/$PRJ_n"
TMP_DIR="$HOME/Documents/Tam_Jenny/$PRJ_n"
RAW="$PRJ_DIR/raw"

mkdir -p $PRJ_DIR/config
mkdir -p $RAW
mkdir -p $TMP_DIR

########################################
# Function: DOWNLOAD_FASTQ per SAMPLE + LANE
########################################
DOWNLOAD_FASTQ () {

  SAMPLE="$1"
  LANE="$2"
  lib="$3"
  ftp="$4"
  md5="$5"

  DIR="$RAW/$SAMPLE/$LANE"
  LOG="$RAW/$SAMPLE/$LANE/${step_id}_${step_n}.log"
  mkdir -p "$DIR"

  {
    echo "====== $(date): Processing $SAMPLE - $LANE ======"

    if [[ "$lib" == "SINGLE" ]]; then

      echo "Downloading $ftp ..."
      wget -O "$DIR/${LANE}.fastq.gz" "$ftp"

      echo "Storing $md5 ..."
      echo "$md5" > $DIR/${LANE}.fastq.gz.md5

    elif [[ "$lib" == "PAIRED" ]]; then

      ftp_R1=$(cut -d";" -f1 "$ftp")
      ftp_R2=$(cut -d";" -f2 "$ftp")
      md5_R1=$(cut -d";" -f1 "$md5")
      md5_R2=$(cut -d";" -f2 "$md5")

      echo "Downloading $ftp_R1 ..."
      wget -O "$DIR/${LANE}_R1.fastq.gz" "$ftp_R1"
      echo "Downloading $ftp_R2 ..."
      wget -O "$DIR/${LANE}_R2.fastq.gz" "$ftp_R2"

      echo "Storing $md5_R1 ..."
      echo "$md5_R1" > $DIR/${LANE}_R1.fastq.gz.md5
      echo "Storing $md5_R2 ..."
      echo "$md5_R2" > $DIR/${LANE}_R2.fastq.gz.md5
    fi

    echo "====== $(date): All done! ======"

  } 2>&1 | tee "$LOG"
    
}

########################################
# export variables and functions globally for parallel
########################################
export RAW
export -f DOWNLOAD_FASTQ

########################################
# Define target samples/lanes for parallel processing
########################################
target="PRJNA1064892"
dl_file="$PRJ_DIR/config/filereport_read_run_${target}.tsv"
sampleID_idx=$(head -1 $dl_file| tr "\t" "\n"| nl| grep -wn "sample_title"| cut -d: -f1)
lane_idx=$(head -1 $dl_file| tr "\t" "\n"| nl| grep -wn "run_accession"| cut -d: -f1)
lib_idx=$(head -1 $dl_file| tr "\t" "\n"| nl| grep -wn "library_layout"| cut -d: -f1)
ftp_idx=$(head -1 $dl_file| tr "\t" "\n"| nl| grep -wn "fastq_ftp"| cut -d: -f1)
md5_idx=$(head -1 $dl_file| tr "\t" "\n"| nl| grep -wn "fastq_md5"| cut -d: -f1)

target_file=$(awk -F'\t' \
  -v sampleID="$sampleID_idx" \
  -v lane="$lane_idx" \
  -v lib="$lib_idx" \
  -v ftp="$ftp_idx" \
  -v md5="$md5_idx" \
  'BEGIN{OFS="\t"} NR > 1 {print $sampleID, $lane, $lib, $ftp, $md5}' \
  "$dl_file")

parallel_log="$TMP_DIR/log/parallel/$target/${step_id}_${step_n}.log"

mkdir -p "$(dirname $parallel_log)"

########################################
# Run parallel
########################################
parallel --joblog "$parallel_log" \
  --resume \
  -j "$n_jobs" \
  --colsep '\t' \
  "DOWNLOAD_FASTQ {1} {2} {3} {4} {5}" \
  ::: "$target_file"