#!/bin/bash
# This script runs the Metaflowmics 16S pipeline with Docker and Nextflow.

set -e

# Path definitions (edit if different)
READS_PATH="/home/$USER/metagenomics/*_{1,2}.fastq.gz"
REF_ALN="/home/$USER/databases/silva.seed_v138_1.align"
REF_TAX="/home/$USER/databases/silva.seed_v138_1.tax"
OUTDIR="/home/$USER/metagenomics_results"

# Clone the pipeline if not already cloned
if [ ! -d ~/metaflowmics ]; then
  git clone https://github.com/hawaiidatascience/metaflowmics ~/metaflowmics
fi

cd ~/metaflowmics/metaflowmics

# Apply custom settings to the config file
sed -i 's/min_overlap = .*/min_overlap = 0/' Pipeline-16S/nextflow.config
sed -i 's/silva_release = .*/silva_release = "seed"/' Pipeline-16S/nextflow.config
sed -i 's/clustering_thresholds = .*/clustering_thresholds = "97"/' Pipeline-16S/nextflow.config
sed -i 's/beta_diversity = .*/beta_diversity = "braycurtis-thetayc-sharedsobs-sharedchao"/' Pipeline-16S/nextflow.config

# Run the pipeline
nextflow run Pipeline-16S -profile docker \
  --reads "$READS_PATH" \
  --referenceAln $REF_ALN \
  --referenceTax $REF_TAX \
  --outdir $OUTDIR \
  --cpus 2
