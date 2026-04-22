# sncRNA Landscape of Mouse Preimplantation Embryos and Human Blastoids

This repository contains the analysis code for profiling small non-coding RNAs (sncRNAs) — including miRNAs, snoRNAs, snRNAs, tRNAs, rRNAs, and piRNAs — in mouse preimplantation embryos and human blastoids at single-cell resolution.

## Overview

We combine multiple single-cell sequencing modalities to characterize the sncRNA landscape during early embryonic development:

- **CoSeq**: Joint profiling of miRNA and mRNA from the same single cell, enabling direct inference of miRNA–mRNA regulatory relationships
- **Small RNA-seq (smallseq)**: Dedicated small RNA library preparation for comprehensive sncRNA quantification
- **SMARTer smt2-seq**: Full-length transcriptome sequencing for mRNA analysis

The project spans two species — **mouse** (oocyte-to-blastocyst stages) and **human blastoids** — In addition with published human embryo data, allowing cross-species comparison of sncRNA dynamics during early lineage segregation.

## Repository Structure

```
├── bin/
│   ├── mouse_coseq_smt2seq/          # SnakeMake pipeline: mouse smt2-seq quant
│   ├── mouse_embryo_smallseq/         # Shell pipeline: mouse small-seq & Co-seq(sncRNA part)
│   └── human_blastoids_batch/
│       ├── HB_smallseq/               # Human blastoid small-seq pipeline & Co-seq(sncRNA part)
│       └── HB_smt2/                   # Human blastoid smt2-seq pipeline (SnakeMake)
│
└── src/
    ├── pre_processing/                # QC filtering, mature expression quantification
    ├── dimensional_reduction/         # UMAP embedding, sub-clustering (ICM/TE), trajectory inference
    ├── coseq_seq/                     # CoSeq-specific: miRNA–gene & miRNA–TF correlation, permutation tests
    ├── DE_miRNA/                      # Differential expression across lineages
    ├── WGCNA/                         # Weighted gene co-expression network analysis (hdWGCNA)
    ├── miRNA_target/                  # miRNA target gene collection and correlation validation
    ├── cross_species/                 # Mouse vs. human sncRNA comparison
    ├── Inheritance/                   # Parental/maternal inheritance analysis
    ├── human_blastoids_coseq_smallseq/  # Human blastoid small RNA-seq analysis
    └── human_blastoids_coseq_smt2/      # Human blastoid smt2-seq analysis
```

## Data Processing Pipelines

### Quantification Pipelines (`bin/`)

All raw sequencing data are processed through standardized workflows:

1. **STAR alignment** — reads are mapped to the reference genome (mm10 for mouse, GRCh38 for human)
2. **RSEM quantification** — gene-level expression estimates (expected counts, TPM, FPKM)
3. **Small RNA classification** — reads are categorized into miRNA, snoRNA, snRNA, tRNA, rRNA, and piRNA based on genomic annotation
4. **Merge and normalize** — sample-level results are merged into unified expression matrices

The SnakeMake workflows (`*.smk`,`*.bash`) automate these steps for batch processing.

### Analysis Workflow (`src/`)

| Step | Scripts | Description |
|------|---------|-------------|
| **QC** | `pre_processing/QC.meta.R` | Filter cells by mitochondrial ratio, total reads, and feature counts |
| **preprocessing** | `pre_processing/create.mature.expression.R` | Generate mature miRNA/sncRNA expression matrices |
| **Visualization** | `dimensional_reduction/` | UMAP, sub-clustering of ICM/TE lineages |
| **Trajectory** | `dimensional_reduction/UMAP.miRNA.traj.R` | Pseudotime inference using Slingshot |
| **CoSeq** | `coseq_seq/` | miRNA–mRNA and miRNA–TF correlation with permutation testing |
| **DEG** | `DE_miRNA/Lineage_segregation.DEG.R` | Differential expression between lineages |
| **Network** | `WGCNA/` | hdWGCNA for miRNA/gene co-expression module detection |
| **miRNA target** | `miRNA_target/` | Integrate miRNA target predictions with CoSeq correlation |
| **Cross-species** | `cross_species/` | Comparison among mouse preimplantation, human blastoid, human preimplantation  sncRNA profiles |
| **Inheritance** | `Inheritance/` | Analyze parental contribution patterns |


## Requirements
- **R ≥ 4.3** with packages: Seurat, ggplot2, dplyr, tidyr, hdWGCNA, WGCNA, slingshot, pheatmap, glmnet, foreach, doParallel
- **Python ≥ 3.6** with: pandas, numpy, pysam
- **STAR**, **RSEM**, **SnakeMake** (for gene quantification pipelines)


## Publication

This code accompanies a study on sncRNA dynamics during early mammalian development. If you use this code, please cite the associated publication.

## License

This project is available for academic and non-commercial use.
