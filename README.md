# COVID-19 Transcriptomic Analysis and Drug Repurposing

## Overview

This project investigates transcriptional changes associated with COVID-19 infection using publicly available RNA-seq data and applies a transcriptomics-based drug repurposing strategy to identify potential therapeutic candidates.

The analysis integrates differential gene expression analysis, functional enrichment analysis, transcriptional signature characterization, and drug–gene perturbation analysis to identify compounds that may potentially reverse disease-associated molecular signatures.

The overall workflow is:

**RNA-seq data → Quality control → Differential expression analysis → Functional enrichment → Disease signature → Drug perturbation analysis → Candidate drug identification**

---

## Research Question

**Which biological pathways and transcriptional changes are associated with COVID-19 infection, and can existing drugs be identified that may potentially reverse these disease-associated molecular signatures?**

The project focuses on comparing COVID-19 samples with healthy controls and subsequently using the resulting molecular signature for drug repurposing.

---

## Objectives

The main objectives of this project are to:

1. Characterize gene-expression changes associated with COVID-19.
2. Identify significantly differentially expressed genes (DEGs).
3. Investigate biological processes and pathways associated with the observed transcriptional changes.
4. Compare COVID-19 transcriptional signatures with perturbational drug-response signatures.
5. Identify existing drugs that may potentially reverse COVID-19-associated gene-expression patterns.
6. Generate reproducible data products and visualizations for downstream interpretation.

---

## Dataset

The analysis uses publicly available RNA-seq data from **GSE163151**.

The repository contains processed count data, sample metadata, normalized expression data, and differential-expression results used throughout the analysis.

### Main comparison

The primary disease comparison is:

**COVID-19 vs Healthy**

An additional comparison between COVID-19 and other viral infections is also incorporated to investigate COVID-19-associated transcriptional signatures beyond general antiviral responses.

---

## Analysis Workflow

### 1. Data preparation

Raw/count-level RNA-seq data and sample metadata are prepared and matched before downstream analysis.

The workflow includes:

* Importing count matrices
* Importing sample metadata
* Matching sample identifiers
* Checking sample consistency
* Removing low-expression genes
* Assessing library size
* Assessing the number of detected genes

Low-expression genes are filtered before differential-expression analysis.

---

### 2. Exploratory data analysis

Several exploratory analyses are performed to assess sample structure and data quality.

These include:

* Library-size analysis
* Number of detected genes
* Principal Component Analysis (PCA)
* Sample correlation analysis
* Correlation heatmaps
* Highly variable gene analysis
* Gene-expression heatmaps

These analyses are used to evaluate sample similarity, group separation, and potential sources of variation.

---

### 3. Differential Expression Analysis

Differential expression analysis is performed using **DESeq2**.

For the main disease comparison:

> **COVID-19 vs Healthy**

Genes are considered significantly differentially expressed using:

* Adjusted p-value (`padj`) < 0.05
* Absolute log2 fold change ≥ 1

Both upregulated and downregulated genes are identified.

The analysis generates:

* Complete DEG table
* Significant DEG table
* MA plot
* Volcano plot
* Top-DEG heatmap
* Normalized expression matrix
* Sample metadata

---

### 4. Functional Enrichment Analysis

Gene Ontology (GO) and pathway-level analyses are performed to characterize the biological functions represented by the differentially expressed genes.

The enrichment analysis focuses on identifying biological processes and pathways associated with the COVID-19 transcriptional response.

The repository contains enrichment results and visualization files, including:

* GO Biological Process results
* GO comparison analyses
* Other-virus-specific GO analysis
* Enrichment visualizations

These analyses help connect individual gene-expression changes with broader biological mechanisms.

---

### 5. Transcriptional Signature Analysis

The significant DEGs are used to define a COVID-19-associated transcriptional signature.

The signature considers both:

* Upregulated genes
* Downregulated genes

This disease-associated signature is subsequently used as the basis for drug perturbation analysis.

The goal is to identify compounds whose transcriptional effects may oppose the observed COVID-19 molecular phenotype.

---

### 6. Drug Repurposing Analysis

A transcriptomic drug repurposing strategy is applied to identify existing compounds that may potentially reverse the COVID-19-associated transcriptional signature.

Drug perturbation information is compared with the disease-associated gene-expression signature.

Candidate compounds are prioritized based on their ability to produce an opposing transcriptional response.

The repository includes drug-repurposing output files containing candidate human drugs and perturbational matching results.

### Conceptual framework

```text
COVID-19 vs Healthy
        │
        ▼
Differentially Expressed Genes
        │
        ├── Upregulated genes
        │
        └── Downregulated genes
        │
        ▼
COVID-19 transcriptional signature
        │
        ▼
Drug perturbation matching
        │
        ▼
Potential signature-reversing compounds
        │
        ▼
Candidate drug repurposing list
```

---

## Repository Structure

```text
Capestone_project/
│
├── GSE163151/
│   └── RNA-seq dataset and related files
│
├── data/
│   └── Input and processed data
│
├── figures/
│   ├── differential_expression/
│   └── enrichment/
│
├── results/
│   ├── differential_expression/
│   ├── signatures/
│   ├── enrichment/
│   └── drug_repurposing/
│
├── COVID19_vs_Healthy_all_DEGs.csv
├── COVID19_vs_Healthy_significant_DEGs.csv
├── COVID19_vs_Healthy_normalized_counts.csv
├── COVID19_vs_Healthy_metadata.csv
│
├── COVID19_drug_repurposing_human_drug_candidates.xlsx
├── ConnectedPerturbations_LIB_5_2026_8_11_15_48_28.xls
│
├── GO.png
├── GO_BP_comparison.png
├── GO_Biological_Process_Comparison.png
├── GO biological process comparison.png
├── OtherVirus_specific_GO_BP.png
├── OtherVirus_specific_GO_BP.csv
├── correlation.png
│
├── Capestone_project.Rmd
├── CMDD_section12_fixed_continuous.Rmd
├── new.Rmd
└── old_capstone.Rmd
```

The current repository contains the main analysis scripts, processed datasets, differential-expression outputs, enrichment results, figures, and drug-repurposing outputs.

---

## Key Output Files

### Differential expression

| File                                       | Description                                                  |
| ------------------------------------------ | ------------------------------------------------------------ |
| `COVID19_vs_Healthy_all_DEGs.csv`          | Complete COVID-19 vs Healthy differential-expression results |
| `COVID19_vs_Healthy_significant_DEGs.csv`  | Significant DEGs based on the predefined thresholds          |
| `COVID19_vs_Healthy_normalized_counts.csv` | DESeq2-normalized expression matrix                          |
| `COVID19_vs_Healthy_metadata.csv`          | Metadata for the COVID-19 vs Healthy comparison              |

The R workflow explicitly exports these datasets after DESeq2 analysis.

### Drug repurposing

| File                                                  | Description                                                            |
| ----------------------------------------------------- | ---------------------------------------------------------------------- |
| `COVID19_drug_repurposing_human_drug_candidates.xlsx` | Candidate human drugs identified through the drug-repurposing workflow |
| `ConnectedPerturbations_LIB_5_2026_8_11_15_48_28.xls` | Perturbational drug-signature matching results                         |

### Functional analysis

| File                                   | Description                                                  |
| -------------------------------------- | ------------------------------------------------------------ |
| `OtherVirus_specific_GO_BP.csv`        | GO Biological Process results for the other-virus comparison |
| `OtherVirus_specific_GO_BP.png`        | Visualization of other-virus-specific GO enrichment          |
| `GO.png`                               | GO enrichment visualization                                  |
| `GO_BP_comparison.png`                 | GO Biological Process comparison                             |
| `GO_Biological_Process_Comparison.png` | Biological-process comparison visualization                  |

---

## Software and Packages

The analysis is primarily implemented in **R** using Bioconductor and CRAN packages.

Key packages include:

* **DESeq2** — differential expression analysis
* **ggplot2** — data visualization
* **dplyr** — data manipulation
* **tidyr** — data processing
* **ComplexHeatmap** — heatmap visualization
* **ggrepel** — non-overlapping labels
* **R Markdown** — reproducible analysis reporting

The main workflow uses DESeq2 for normalization and differential-expression analysis, including variance-stabilizing transformation, PCA, sample correlation analysis, and DEG identification.

---

## Reproducibility

The primary analysis is documented in:

```text
CMDD_section12_fixed_continous.Rmd
```

Additional analysis and drug-repurposing workflow components are documented in:

```text
CMDD_section12_fixed_continuous.Rmd
```

The R Markdown files contain the analysis code required to reproduce the major data-processing, differential-expression, visualization, enrichment, and downstream analyses.

---

## Statistical Criteria

For the main COVID-19 vs Healthy differential-expression analysis, significant genes are defined as:

```text
Adjusted p-value < 0.05
AND
|log2 Fold Change| ≥ 1
```

Genes satisfying:

```text
log2 Fold Change ≥ 1
```

are classified as upregulated, while genes satisfying:

```text
log2 Fold Change ≤ -1
```

are classified as downregulated.

---

## Results Interpretation

The analysis provides a multi-level view of the COVID-19 transcriptional response:

### Gene level

Differentially expressed genes identify individual genes whose expression differs between COVID-19 and healthy samples.

### Pathway level

GO and pathway enrichment analyses identify biological processes associated with these transcriptional changes.

### Signature level

The combined upregulated and downregulated genes define a molecular signature associated with COVID-19.

### Drug level

Drug perturbation matching identifies compounds whose transcriptional effects may oppose the disease-associated signature.

Together, these analyses provide a computational framework for prioritizing existing drugs for further investigation.

---

## Limitations

This project is a **computational drug-repurposing study** and does not demonstrate clinical efficacy.

Important limitations include:

* Transcriptomic reversal does not necessarily imply therapeutic efficacy.
* Drug perturbation signatures may differ across cell types, tissues, doses, and experimental conditions.
* Public RNA-seq datasets may contain biological and technical heterogeneity.
* Candidate drugs require independent validation.
* Computational prioritization should not be interpreted as clinical treatment recommendations.

Experimental validation and clinical evidence are required before any candidate can be considered a validated therapeutic option.

---

## Future Work

Potential future extensions include:

1. Validation of top drug candidates using independent COVID-19 datasets.
2. Cross-dataset replication of the differential-expression signature.
3. Integration of additional drug perturbation databases.
4. Protein–protein interaction network analysis.
5. Identification of hub genes and key regulatory pathways.
6. Molecular docking or structural analysis for prioritized drug–target pairs.
7. Comparison with experimentally validated COVID-19 therapeutics.
8. Development of a reproducible drug-ranking pipeline.
9. External validation of the highest-ranked candidates.

---

## Project Significance

This project demonstrates how transcriptomic data can be integrated with computational drug-repurposing approaches to move from:

**Disease-associated gene expression**

to

**biological mechanism**

and ultimately to

**potential therapeutic candidates**.

The workflow provides a reproducible computational framework for exploring how existing drugs might be repositioned for emerging infectious diseases.

---

## Author

**Jiayi Liang**

Capstone Project — Bioinformatics / Computational Biology

GitHub:
[jiayi-Liang025/Capestone_project](https://github.com/jiayi-Liang025/Capestone_project?utm_source=chatgpt.com)

---

## Disclaimer

This repository is intended for academic and research purposes only.

Drug candidates identified through computational analysis are hypotheses for further investigation and should not be interpreted as clinically validated treatments or medical recommendations.
