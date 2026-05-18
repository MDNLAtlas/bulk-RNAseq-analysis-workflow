#!/usr/bin/bash

REF="$HOME/myRDS/PRJ-MRefD/databases"

# Gunzip GTF primary assembly for genome index
gunzip "$REF/Gencode/v48/gencode.v48.primary_assembly.annotation.gtf.gz"
gunzip "$REF/Gencode/v48/GRCh38.primary_assembly.genome.fa.gz"

# Create gene table for GTF basic genes only
GTF="$REF/Gencode/v48/gencode.v48.basic.annotation.gtf.gz"

header="chr\tstart\tend\tgeneID\tscore\tstrand\tgene_symbol\tgene_type"
(echo -e "$header"; zcat "$GTF" | \
  grep -v "^#"| \
  awk -F'\t' '{OFS="\t"}{split($9,a,"; ")}{if($3 == "gene") print $1, $4, $5, a[1], $6, $7, a[3], a[2]}' | \
  sed 's/gene_id "//g' | \
  sed 's/gene_type "//g' | \
  sed 's/gene_name "//g' | \
  sed 's/"//g') > "${GTF/gtf.gz/genes.tsv}"

header="chr\tstart\tend\ttranscriptID\tscore\tstrand\tgeneID\tgene_symbol\tgene_type\ttranscript_type"
(echo -e "$header"; zcat "$GTF" | \
  grep -v "^#" | \
  awk -F'\t' '{OFS="\t"}{split($9,a,"; ")}{if($3 == "transcript") print $1, $4, $5, a[2], a[1], a[4], a[3], a[5]}' | \
  sed 's/gene_id "//g' | \
  sed 's/transcript_id "//g' | \
  sed 's/gene_name "//g' | \
  sed 's/gene_type "//g' | \
  sed 's/transcript_type "//g' | \
  sed 's/"//g') > "${GTF/gtf.gz/transcripts.tsv}"