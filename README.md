# Assessment 2: RNA-Seq analysis

## Author
Student number: 200877385
## Description
## Requirements
### R Version
- R version 4.5.0 (or higher)

### R Packages:
- `BiocManager`
- `tximport` 
- `DESeq2` 
- `biomaRt` 
- `tidyverse`
  - `ggplot2` 
  - `dplyr` 
  - `readr` 
- `pheatmap`
- `RColorBrewer`
- `patchwork`
- `ggVennDiagram`

## Process
All data is downloaded and processed by the code in the 'scripts/assessment_2.R'.
The script should be run from the base directory.

The script will:
1. Create necessary directories
2. Download required data files
3. Perform RNA-Seq analysis
4. Generate all output figures and tables

## Data
The RNA sequencing (RNA-seq) data used is a publicly available data set that can be accessed from the Gene Expression Omnibus (GEO), accession number GSE116583. 
The reference transcriptome used is the GENCODE mouse, release M31.

## Directory stucture and outputs
~~~
├── assessment-2.Rproj
├── data/
│   ├── counts/                   # Salmon outputs - quant.sf for each sample
│   │   ├── SRR7457551/
│   │   ├── SRR7457552/
│   │   ├── ..........
│   │   └── SRR7457562/
│   ├── counts.zip
│   ├── gene_map.csv              # Transcript to gene mapping data
│   └── sample_table.csv          # Sample metadata (Run, ID, Timepoint ect.) - Downloaded from module github repro
├── results/
│   ├── annot_results_24.csv      # Annotated results 24h
│   ├── annot_results_2h.csv      # Annotated results 2h
│   ├── dds.rds                   # DESeq2 dataset
│   ├── degs_24h.csv              # DEGs 24h
│   ├── degs_2h.csv               # DEGs 2h
│   ├── figures/
│   │   ├── degs_24h_heatmap.pdf                  # Euclidean distance heatmap for 250 most differentially expressed genes at 24h
│   │   ├── degs_2h_heatmap.pdf                   # Euclidean distance heatmap for 250 most differentially expressed genes at 2h
│   │   ├── dispersion_plot.pdf                   # Dispersions estimate plot
│   │   ├── figure4_combined_volcano.pdf          # Multipanal volcano plot for 2 and 24 hour groups
│   │   ├── pca_plot.pdf                          # PCA plot
│   │   ├── sample_distance_heatmap.pdf           # Clustering Heatmap of all samples
│   │   ├── venndiagram.pdf                       # Venn Diagram of all DEGs across 2 and 24 hour timepoints
│   │   ├── volcano_plot_24h.pdf                  # Volcano plot of differential gene expression at 24 hours
│   │   └── volcano_plot_2h.pdf                   # Volcano plot of differential gene expression at 2 hours
│   ├── results_24h_table.csv                    # DESeq2 results of Allo24h vs Naive
│   ├── results_2h_table.csv                     # DESeq2 results of Allo2h vs naive
│   └── rld.rds                                  # r-log transformed counts
└── scripts/
    ├── assessment_2.R
    └── directory_structure.R
~~~
