#!/bin/bash
# This script downloads paired-end FASTQ files from SRA using accession numbers
# and metadata from a CSV file. It trims each file to ~100k reads for uniformity.

set -e  # Exit immediately if a command exits with a non-zero status

# Path to the metadata CSV file obtained from NCBI SRA Run Selector
METADATA_FILE=~/metagenomics/metadata.csv

cd ~/metagenomics

# Remove header and iterate over each line in the CSV
sed '1d' $METADATA_FILE | while read meta; do
  # Extract accession number and sample name from the CSV fields
  accn=$(echo $meta | cut -d "," -f 1)
  name=$(echo $meta | cut -d "," -f 43)

  echo "Downloading $accn ..."
  fasterq-dump $accn

  # Keep only the first 100k reads (800,000 lines = 100,000 reads)
  echo "Trimming and compressing reads for $name ..."
  head -n 800000 ${accn}_1.fastq | gzip > ${name}_1.fastq.gz
  rm ${accn}_1.fastq

  head -n 800000 ${accn}_2.fastq | gzip > ${name}_2.fastq.gz
  rm ${accn}_2.fastq
done
