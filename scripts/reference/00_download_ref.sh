#!/usr/bin/bash

# shellcheck disable=SC2086

Gencode="$HOME/myRDS/PRJ-MRefD/databases/Gencode"
mkdir -p $Gencode/v48

# Download GTF primary assembly
wget -P "$Gencode/v48" \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz

# Download GTF basic genes only
wget -P "$Gencode/v48" \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.basic.annotation.gtf.gz

# wget -P "$Gencode/v48" \
#   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz

# Download FASTA primary assembly
wget -P "$Gencode/v48" \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz

gunzip $Gencode/v48/GRCh38.primary_assembly.genome.fa.gz

grep "^>" "$Gencode/v48/GRCh38.primary_assembly.genome.fa.gz" | \
  sed 's/>//' > "$Gencode/v48/GRCh38.primary_assembly.contigs.txt"

# Ensembl="$HOME/myRDS/PRJ-MRefD/databases/Ensembl/v115"
# wget -P "$Ensembl/v48" \
#   https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

