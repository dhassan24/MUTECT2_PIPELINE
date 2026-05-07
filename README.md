# MUTECT2_PIPELINE
Filtering tumor/germline paired variants using MUTECT2 program and detecting somatic SNVs. Docker images and .py filtration pipeline included.
# Somatic SNV Detection Pipeline
### Mutect2 · FACETS · Docker
****

---

## Overview

This repository contains a reproducible computational pipeline for detecting somatic single nucleotide variants (SNVs) from paired tumor/normal BAM files using **GATK4 Mutect2** and **FACETS** copy number analysis. The pipeline is containerized with **Docker** to ensure consistent, reproducible results across computing environments.

The dataset used is the **SMC-HET T2** sample from the ICGC-TCGA DREAM Somatic Mutation Calling Heterogeneity Challenge — a synthetic tumor/normal pair designed to benchmark somatic variant callers under realistic tumor heterogeneity conditions.

---

## Repository Structure

```
MUTECT2_PIPELINE/
├── Docker/
│   └── Dockerfile              # Builds the analysis environment
├── filter_variants.py          # Python filtering script
├── filter_variants.R           # R filtering script
└── README.md                   # This file
```

---

## Docker Image

The Docker image is built on **Ubuntu 18.04** and contains:

| Tool | Purpose |
|------|---------|
| GATK4 (v4.5.0.0) | Somatic SNV calling with Mutect2 |
| samtools (v1.19) | BAM file manipulation and indexing |
| FACETS (R package) | Allele-specific copy number analysis |
| snp-pileup | Generates allele counts from BAMs for FACETS |
| facets-suite | Command-line wrappers for FACETS |
| Python 3.7 + pandas + cyvcf2 | Variant filtering scripts |
| R + Rscript | FACETS R package execution |

### Build the Image

```bash
docker build -t dh-gatk-facets:v1.0 .
```

### Run Interactively

```bash
docker run -it -v /path/to/your/files:/data dh-gatk-facets:v1.0
```

The `-v` flag mounts your local files into `/data` inside the container — no need to copy files into the image.

---

## Filtering Pipeline

### Background

Mutect2 calls somatic variants by comparing tumor and normal BAMs using a Bayesian statistical model. However, raw calls still contain false positives that require additional filtering. FACETS provides allele-specific copy number estimates across the tumor genome, which directly informs how we interpret Mutect2's allele frequency (AF) calls:

- In **amplified regions** (high total copy number), the same variant appears at a lower AF than expected — making it look like noise
- In **LOH regions** (minor copy number = 0), one allele is fully lost — AF expectations break down entirely

By integrating FACETS copy number output with Mutect2 calls, we can flag variants in genomically complex regions that are unreliable regardless of their AF.

### Filters Applied

| Filter | Threshold | Variants Removed |
|--------|-----------|-----------------|
| Mutect2 PASS | Keep FILTER=PASS only | 0 |
| Minimum Allele Frequency | AF < 0.05 | 1,108 |
| Minimum Read Depth | DP < 10 | 8 |
| FACETS High Copy Number | tcn.em > 8 | 0 |
| FACETS LOH Region | lcn.em == 0 | 0 |
| **Final PASS Variants** | | **4,282** |

### Run the Python Script

```bash
docker run -it -v /path/to/files:/data dh-gatk-facets:v1.0 \
  python /data/filter_variants.py \
  --vcf /data/HN002.vcf.gz \
  --facets_seg /data/HG002-T2_hisens.seg \
  --facets_qc /data/HG002-T2.qc.txt \
  --output /data/HN002_filtered.vcf
```

### Run the R Script

```bash
docker run -it -v /path/to/files:/data dh-gatk-facets:v1.0 \
  Rscript /data/filter_variants.R \
  --vcf /data/HN002.vcf.gz \
  --facets_seg /data/HG002-T2_hisens.seg \
  --facets_qc /data/HG002-T2.qc.txt \
  --output /data/HN002_filtered.vcf
```

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_af` | 0.05 | Minimum tumor allele frequency |
| `--min_depth` | 10 | Minimum tumor read depth |
| `--max_tcn` | 8 | Maximum total copy number (FACETS) |

---

## Input Files

| File | Source | Description |
|------|--------|-------------|
| `HN002.vcf.gz` | Mutect2 output | Raw somatic SNV calls |
| `HN002.vcf.stats` | Mutect2 output | Variant calling statistics |
| `HG002-T2_hisens.seg` | FACETS output | High-sensitivity copy number segments |
| `HG002-T2_purity.seg` | FACETS output | Purity-optimized segments |
| `HG002-T2.qc.txt` | FACETS output | Purity and ploidy QC estimates |
| `HG002-T2.txt` | FACETS output | Full FACETS results |
| `HG002-T2.arm_level.txt` | FACETS output | Chromosome arm-level CN calls |
| `HG002-T2.gene_level.txt` | FACETS output | Gene-level CN calls |

---

## Additional Filtering Strategies (With Germline Data)

If additional data is available, false positive filtering can be further improved:

- **Panel of Normals (PoN)** — Build a VCF from 20-40 unrelated normals. Recurrent artifacts across normals are excluded from somatic calls
- **gnomAD population database** — Variants present at >0.1% population frequency are likely germline polymorphisms
- **Matched germline VCF** — Directly subtract patient-specific germline variants; most powerful method for catching rare inherited variants
- **Read orientation model** — GATK's `LearnReadOrientationModel` corrects for OxoG and FFPE artifacts

---

## Dependencies

All dependencies are pre-installed in the Docker image. If running outside Docker:

```bash
# Python
pip install pandas cyvcf2 numpy

# R
Rscript -e "install.packages(c('optparse', 'data.table', 'dplyr', 'vcfR'))"
```

---

## Author

**Danya Hassan**
Keck Graduate Institute
MS Human Genetics and Genomic Data Analytics
