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

## Output

## Directory stucture
├── assessment-2.Rproj
├── data
│   ├── counts
│   │   ├── SRR7457551
│   │   ├── SRR7457552
│   │   ├── SRR7457553
│   │   ├── SRR7457554
│   │   ├── SRR7457555
│   │   ├── SRR7457556
│   │   ├── SRR7457557
│   │   ├── SRR7457558
│   │   ├── SRR7457559
│   │   ├── SRR7457560
│   │   ├── SRR7457561
│   │   └── SRR7457562
│   ├── counts.zip
│   ├── gene_map.csv
│   └── sample_table.csv
├── results
│   ├── annot_results_24.csv
│   ├── annot_results_2h.csv
│   ├── dds.rds
│   ├── degs_24h.csv
│   ├── degs_2h.csv
│   ├── figures
│   │   ├── degs_24h_heatmap.pdf
│   │   ├── degs_2h_heatmap.pdf
│   │   ├── dispersion_plot.pdf
│   │   ├── figure4_combined_volcano.pdf
│   │   ├── pca_plot.pdf
│   │   ├── sample_distance_heatmap.pdf
│   │   ├── venndiagram.pdf
│   │   ├── volcano_plot_24h.pdf
│   │   └── volcano_plot_2h.pdf
│   ├── filtered_results_24.csv
│   ├── results_24h_table.csv
│   ├── results_2h_table.csv
│   └── rld.rds
└── scripts
    ├── assessment_2.R
    └── directory_structure.R