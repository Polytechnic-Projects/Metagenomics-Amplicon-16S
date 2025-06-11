<p align="center">
<a href="https://www.r-project.org/" title="Go to R homepage"><img src="https://img.shields.io/badge/Made%20with-R-FED564?logo=r&logoColor=white" alt="Made with R"></a>
<a href="https://www.ncbi.nlm.nih.gov/sra" title="NCBI SRA">
  <img src="https://img.shields.io/badge/Data%20from-NCBI%20SRA-FED564?logo=databricks&logoColor=white" alt="Data from NCBI SRA">
</a>
</p>

<h3 align="center">📊 Metagenomics-Amplicon-16S - Comparative Analysis of Hydrothermal Vent Microbiomes</h3>
<p align="center">
This project analyzes microbial diversity in hydrothermal vent environments using 16 paired-end amplicon sequencing samples processed with the Metaflowmics pipeline and visualized using R.
</p>

---

## Authors

A group of 4 students from ESTBarreiro with their school id number assigned.

- Ana Fernando, 202101267
- Erica Alaiz, 202300154
- Érica Martins, 202300387
- Melissa Rocha, 202101023


---

## 👋 Index

* [📅 Introduction](#-introduction)
* [⚙️ Features](#-features)
* [💻 Installation](#-installation)
* [🚀 Usage](#-usage)
* [📊 File Structure](#-file-structure)
* [📄 Sample Accession Numbers](#-sample-accession-numbers)
* [📜 License](#-license)
* [🔗 References](#-references)

---

## 📅 Introduction

This repository contains a fully reproducible pipeline for processing 16S rRNA amplicon sequencing data from the BioProject **PRJNA1080611**. The dataset consists of 16 paired-end samples collected from hydrothermal vents and sequenced using two primer sets (V3V4 and V4V5). We performed read preprocessing, ASV inference, taxonomic classification, and diversity analysis.

The goal is to compare the microbial diversity captured by each primer set and assess the community composition across samples.

---

## ⚙️ Features

* 👁️ Metadata-driven project with clean sample annotation
* 📦 Automatic paired-end FASTQ download from SRA
* ⚙️ DADA2 ASV inference and SILVA v138.1 classification via Metaflowmics
* 📈 Alpha (Chao1, Shannon, Simpson) and beta diversity (Bray-Curtis + PERMANOVA)
* 📊 Taxonomic composition (bar & pie charts), indicator species (IndVal), Zetaproteobacteria profile
* 🇷 Fully commented and modular R script

---

## 💻 Installation

### Prerequisites

* Linux-based system
* `Docker`, `Nextflow`, `Java` (`default-jre`)
* `R` (>= 4.5.0)

### Installation Steps

1. Clone this repository:

```sh
git clone https://github.com/Polytechnic-Projects/Metagenomics-Amplicon-16S.git
cd Metagenomics-Amplicon-16S
```

2. Install all dependencies and the SILVA database:

```sh
bash scripts/install_dependencies.sh
```

> Note: The SILVA v138.1 seed database will be downloaded to `~/databases/`

---

## 🚀 Usage

### Step 1: Download and filter FASTQ reads from SRA

```sh
bash scripts/download_and_prepare_reads.sh
```

This step downloads the 16 samples (32 FASTQ files), filters them to \~100k reads, and renames them using metadata from `metadata.csv`.

### Step 2: Run the Metaflowmics pipeline with Docker

```sh
bash scripts/run_pipeline_nextflow.sh
```

This will process the reads using DADA2, assign taxonomy using SILVA, and generate count and taxonomy tables.

### Step 3: Run the R analysis script

```sh
Rscript scripts/analysis.R
```

This performs all statistical analyses (alpha/beta diversity, composition, indicator species) and saves all plots in the `figures/` directory.

---

## 📊 File Structure

```
Metagenomics-Amplicon-16S/
├── figures/               # Output plots (alpha, beta, taxonomy)
├── metadata.tsv           # For mapping and analysis
├── accessions.csv         # Optional CSV mapping of SampleID to SRR IDs
├── results/               # OTU and taxonomy output from Metaflowmics
├── postprocessing/        # Alpha and beta diversity output from Metaflowmics
├── scripts/               # All bash and R scripts
│   ├── install_dependencies.sh
│   ├── download_and_prepare_reads.sh
│   ├── run_pipeline_nextflow.sh
│   └── analysis.R
├── LICENSE
└── README.md
```

---

## 📄 Sample Accession Numbers

| SampleID        | SRA Accession |
|-----------------|----------------|
| Bact674Blue     | SRR28098825    |
| Bact674Green    | SRR28098824    |
| Bact675Black    | SRR28098817    |
| Bact676Black    | SRR28098816    |
| Bact676Green    | SRR28098815    |
| Bact677Black    | SRR28098814    |
| Bact677Blue     | SRR28098813    |
| Bact677Green    | SRR28098812    |
| Uni674Blue      | SRR28098811    |
| Uni674Green     | SRR28098810    |
| Uni675Black     | SRR28098823    |
| Uni676Black     | SRR28098822    |
| Uni676Green     | SRR28098821    |
| Uni677Black     | SRR28098820    |
| Uni677Blue      | SRR28098819    |
| Uni677Green     | SRR28098818    |

---
## 📜 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## 🔗 References

* [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)
* [Metaflowmics](https://github.com/hawaiidatascience/metaflowmics)
* [DADA2](https://benjjneb.github.io/dada2/)
* [SILVA Database](https://www.arb-silva.de/)
* [phyloseq](https://joey711.github.io/phyloseq/)
* [vegan](https://cran.r-project.org/web/packages/vegan/index.html)
* Smith et al. (2024), Microbial communities at hydrothermal vents
