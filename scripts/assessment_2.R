# Load libraries ----
library(tidyverse)
library(tximport)
library(DESeq2)
library(ggplot2) # part of tidyverse
library(pheatmap)
library(RColorBrewer)   
library(biomaRt) 
library(ggVennDiagram)
library(patchwork)
library (ggrepel)
# Colour schemes: ----
group_colours = c("Allo24h" = "#264653", 
                  "Allo2h" = "#e9c46a", 
                  "Naive" = "#e76f51")

heatmap_colours = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                                     "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D",
                                     "#B2182B", "#67001F"))(100)
#create directories if they dont already exist - use vector and for loop ----
print("Creating data, results and scripts directories:")
dir_names = c("data", "results", "scripts")
for (a in dir_names) {
  if(!dir.exists(a)) {
    dir.create(a)
  } else {
    message(a, " directory already exists.")
  }
}
dir.create("results/figures")
#download counts.zip files from github ----
print("Downloading counts data:")
if (!file.exists("data/counts.zip")) {
  download.file('https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/results/counts.zip', destfile = 'data/counts.zip')
  unzip ('data/counts.zip', exdir = 'data/')
} else {
  print("counts.zip already exists.")
}
#create sample_table - 
print("Creating data frame, sample_table, that includes metadata about each sample, and downloading sample_table.csv into the data directory.")
sample_table = read_csv('https://raw.githubusercontent.com/sjcockell/mmb8052/main/practicals/practical_08/data/sample_table.csv')
if (!file.exists("data/sample_table.csv")) {
  write.csv(sample_table, file = "data/sample_table.csv", row.names = FALSE)
} else {
  print("sample_table.csv has already been downloaded.")
}
# create filepath vector ----
#paste0, pasting with no gap the RUN column from sample_table
print("Creating a vector of the file names within the counts directory.")
files = pull(sample_table, Run)
files = paste0('data/counts/', files, '/quant.sf')
names(files) = pull(sample_table, Run)
#import gene map into a dataframe ----
print("Creating data frame called gene_map that maps transcript IDs to gene IDs, and downloading gene_map.csv into the data directory.")
gene_map = read_csv('https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/extdata/gene_map.csv')
if (!file.exists("data/gene_map.csv")) {
  write.csv(gene_map, file = "data/gene_map.csv", row.names = FALSE)
} else {
  print("gene_map.csv has already been downloaded.")
}

#create txi object ----
#files is a vector of filepaths to quants.sf files
#the _quants.sf files contain data outputed by salmon
#creates a lit of 4 matrices: abundance, counts, lengths and countsfromabundance
txi = tximport(files, 
               type='salmon',
               tx2gene=gene_map,
               ignoreTxVersion=TRUE)
#Rename samples ----
renamed_samples = c("24h_1","24h_2","2h_3","2h_4","2h_1","2h_2","N_3","N_4","N_1","N_2","24h_3","24h_4")
sample_table$Sample_Name = renamed_samples
rownames(sample_table) = renamed_samples
colnames(txi$counts) = renamed_samples
colnames(txi$abundance) = renamed_samples
colnames(txi$length) = renamed_samples

#normalisation of data ----
#dds DESeq - Differential Expression analysis of Sequence count data. dds - DESeq Data set
#DESeqDataSetFromTximport is a function that uses the data stored in the txi to create a DESeq data set
#then adding size factor, duspersuions and bionomial wald test witin DDS
dds = DESeqDataSetFromTximport(txi, colData = sample_table_renamed, design = ~ Group)
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)
print("Saving DESeq Data set to results/dds.rds")
if (!file.exists("results/dds.rds")) {
  saveRDS(dds, "results/dds.rds")
} else {
  print("dds.rds has already been saved.")
}

#rlog transformation ----
#rld= regulerized logerthymic transformation - makes the varience more stable
print("Saving the regularized logerythmic transformation of the DEseq data set to results as, rld.rds")
rld = rlog(dds)
if (!file.exists("results/rld.rds")) {
  saveRDS(rld, "results/rld.rds")
} else {
  print("rdk.rds has already been saved.")
}
#PCA PLOT ----
plotPCA(rld, intgroup='Group') +
  ggtitle("PCA of Gene Expression Across Timepoints") +
  xlab("PC1 (60% of variance)") +
  ylab("PC2 (18% of variance)") +
  labs(colour = "Timepoint") +
  geom_point(size =4) +
  scale_colour_manual(values = group_colours,
                      labels = c("Allo 24h", "Allo 2h", "Naive")) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text =element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )
ggsave("results/figures/figure1_pca_plot.pdf")
#Heatmap ----
#create a distance object containing all sample-sample distances - take the matrix of rld, transpose it - so samples are rowa, calcylate distancw using euxlidean method
sample_distance = dist(t(assay(rld)), method='euclidean')
#convert sample_distance to a matrix - rows and columns are samples - each cell is the distance between the 2 samples
sample_distance_matrix = as.matrix(sample_distance)
#create a data frame with the labels for the heatmap- one row per sample, one column - group
heatmap_annotation = data.frame(
  Timepoint = colData(dds)[,c('Group')],
  row.names=rownames(colData(dds))
)
annotation_colors = list(Timepoint = group_colours)
diag(sample_distance_matrix) = NA
pheatmap(sample_distance_matrix,
         color =heatmap_colours,
         clustering_distance_rows=sample_distance,
         clustering_distance_cols=sample_distance,
         annotation_col = heatmap_annotation,
         annotation_row = heatmap_annotation,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         fontsize = 10,
         border_color = "black",
         cellwidth = 20,
         cellheight = 20,
         angle_col = 45,
         na_col = "white",
         main = "Euclidean Distance Heatmap of rlog-transformed samples",
         filename = "results/figures/figure2_sample_distance_heatmap.pdf")
#Dispersion estimate plot ----
pdf("results/figures/figure3_dispersion_plot.pdf", width = 8, height = 6)
plotDispEsts(dds, 
             main = "Dispersion Estimates",
             xlab = "Mean of normalized counts",
             ylab = "Dispersion")
dev.off()
#24hour vs niave comparasion ----
# log2FoldChange > 0 → higher expression in Allo24h, log2FoldChange < 0 → higher expression in Naive
results_24h_table = results(dds, contrast= c('Group', 'Allo24h', 'Naive'))
summary(results_24h_table)
write.csv(results_24h_table, file = "results/results_24h_table.csv", row.names = FALSE)
# create a tibble - row names in results_table were gene ids, create a column for gene ids
results_24h_tibble = as_tibble(results_24h_table, rownames='ensembl_gene_id')
# remove rows with missing data
filtered_results_24h = filter(results_24h_tibble, complete.cases(results_24h_tibble))
#create a -log10padj column - makes the padj larger and positive - easier to see differences
filtered_results_24h = mutate(filtered_results_24h, logPVal = -log10(padj))
#Add signifiance true/false column
filtered_results_24h = mutate(filtered_results_24h,
                              Significance = ifelse(log2FoldChange > 1 & padj < 0.05, "Up-regulated",
                                                    ifelse(log2FoldChange < -1 & padj < 0.05, "Down-regulated",
                                                           "Not significant"))
)
#2h vs naive comparison ----
#log2FoldChange > 0 → higher expression in Allo24h, log2FoldChange < 0 → higher expression in Naive
results_2h_table = results(dds, contrast= c('Group', 'Allo2h', 'Naive'))
summary(results_2h_table)
write.csv(results_2h_table, file = "results/results_2h_table.csv", row.names = FALSE)
# create a tibble - row names in results_table were gene ids, create a column for gene ids
results_2h_tibble = as_tibble(results_2h_table, rownames='ensembl_gene_id')
# remove rows with missing data
filtered_results_2h = filter(results_2h_tibble, complete.cases(results_2h_tibble))
filtered_results_2h = mutate(filtered_results_2h, logPVal = -log10(padj))
#Add signifiance true/false column
filtered_results_2h = mutate(filtered_results_2h,
                             Significance = ifelse(log2FoldChange > 1 & padj < 0.05, "Up-regulated",
                                                   ifelse(log2FoldChange < -1 & padj < 0.05, "Down-regulated",
                                                          "Not significant"))
)
#Annotate results with ensembl gene anotations ----
#use verion 108 and get the mmusculus gene anotations
ensembl108 = useEnsembl(biomart="ensembl", version=108)
ensembl108 = useDataset("mmusculus_gene_ensembl", mart=ensembl108)
#create data frame with the following annotations - filter by ensemb gene id that show up in filtered_results
annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 
                                'start_position', 'end_position', 
                                'strand', 'gene_biotype', 'external_gene_name',
                                'description'),
                   filters = 'ensembl_gene_id', values = filtered_results_24h$ensembl_gene_id,
                   mart = ensembl108)
#annotate 24h filtered results
annot_results_24h = left_join(filtered_results_24h, annotation)
#arrange annoted_results by padj
annot_results_24h = arrange(annot_results_24h, padj)
write.csv(annot_results_24h, file = "results/annot_results_24.csv", row.names = FALSE)
#View(head(annot_results_24h, 10))
#Annotate 2h filtered results
annot_results_2h = left_join(filtered_results_2h, annotation)
annot_results_2h = arrange(annot_results_2h, padj)
write.csv(annot_results_2h, file = "results/annot_results_2h.csv", row.names = FALSE)
#View(head(annot_results_24h, 10))
#Volcano plot colour scheme ----
volcano_colours = c("Up-regulated" = "#9d0208", 
                    "Down-regulated" = "#003049",  
                    "Not significant" = "grey")
#show top 10 significant degs ----
top_24h = bind_rows(degs_24h %>%
                      arrange(padj) %>%
                      slice_head(n = 10))

top_2h = bind_rows(degs_2h %>%
                     arrange(padj) %>%
                     slice_head(n = 10))
#24hr comparison volcano plot ----
#plot log2 fold change against logpVal - plot of how much it changes vs statistical significance
volcano_24h = ggplot(filtered_results_24h, aes(x=log2FoldChange, y=logPVal, colour = Significance)) + 
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = volcano_colours) +
  # Add threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5) +
  #Add Axis lables
  labs(
    title = "Allo 24h vs Naive",
    x = "log2 Fold Change",
    y = "-log10 (adjusted p-value)") +
  #label the sig genes
  geom_label_repel(data = top_24h, aes(label = external_gene_name),
                   fill = "white", colour = "grey5", size = 3, 
                   box.padding = 0.3, point.padding = 0.2,
                   max.overlaps = Inf,
                   segment.color = "grey40", segment.size = 0.4) +
  theme_bw() + #white background
  xlim(-15, 15) +  
  ylim(0, 150)
ggsave("results/figures/volcano_plot_24h.pdf",
       volcano_24h)
#plot 2hr comparison volcano plot ----
volcano_2h = ggplot(filtered_results_2h, aes(x=log2FoldChange, y=logPVal, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = volcano_colours) +
  # Add threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(
    title = "Allo 2h vs Naive",
    x = "log2 Fold Change",
    y = "-log10 (adjusted p-value)") +
  geom_label_repel(data = top_2h, aes(label = external_gene_name),
                   fill = "white", colour = "grey5", size = 3, 
                   box.padding = 0.3, point.padding = 0.2,
                   max.overlaps = Inf,
                   segment.color = "grey40", segment.size = 0.4) +
  theme_bw() +
  xlim(-15, 15) + 
  ylim(0, 150)
ggsave("results/figures/volcano_plot_2h.pdf",
       volcano_2h)
#combined volcano plot figure ----
figure_4 = volcano_2h + volcano_24h +
  plot_annotation(
    title = "Figure 4: Differential Gene Expression for Allo2h and Allo24h vs Naive",
    tag_levels = "A"
  ) & #'&' applys theme to all sub plots
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave("results/figures/figure4_combined_volcano.pdf",
       figure_4,
       width = 12, height = 6)

#Create DEGs tables ----
#degs for 24h
#absolute fold change (pos or neg >1) and adj p value <0.05
degs_24h = filter(annot_results_24h, abs(log2FoldChange) > 1 & padj < 0.05)
write.csv(degs_24h, file = "results/degs_24h.csv", row.names = FALSE)
#degs for 2h
#degs = differentially expressed genes, absolute fold change (pos or neg >1) and adj p value <0.05
degs_2h = filter(annot_results_2h, abs(log2FoldChange) > 1 & padj < 0.05)
write.csv(degs_2h, file = "results/degs_2h.csv", row.names = FALSE)
#Create Degs venn diagram ----
#unique removes replicates
degs_24h_list = unique(degs_24h$external_gene_name)
degs_2h_list = unique(degs_2h$external_gene_name)
# Calculate overlaps
shared_genes = intersect(degs_24h_list, degs_2h_list)
unique_24h = setdiff(degs_24h_list, degs_2h_list)
unique_2h = setdiff(degs_2h_list, degs_24h_list)
# Print statistics
cat("DEGs unique to 24h:", length(unique_24h), "\n")
cat("DEGs unique to 2h:", length(unique_2h), "\n")
cat("DEGs shared:", length(shared_genes), "\n")
cat("Total DEGs (24h):", length(degs_24h_list), "\n")
cat("Total DEGs (2h):", length(degs_2h_list),"\n")
#Create venndiagram
ggVennDiagram(list("24h" = degs_24h_list, "2h" = degs_2h_list)) +
  scale_fill_gradient(low = "white", high = "#606c38", name = "Counts") +
  ggtitle("DEGs Overlap at 2h and 24h") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    text = element_text(size = 6),
    plot.margin = unit(c(2, 2, 2, 2), "cm")
  )
ggsave("results/figures/venndiagram.pdf")
# Up and downregulated degs venndiagram ----
# Upregulated ----
# Create objects containing only significant pos or neg fold change
up_degs_2h = filter(degs_2h, log2FoldChange > 0)  
down_degs_2h = filter(degs_2h, log2FoldChange < 0) 
up_degs_24h = filter(degs_24h, log2FoldChange > 0)  
down_degs_24h = filter(degs_24h, log2FoldChange < 0)  
#Create lists
up_degs_2h_list = unique(up_degs_2h$external_gene_name)
down_degs_2h_list = unique(down_degs_2h$external_gene_name)
up_degs_24h_list = unique(up_degs_24h$external_gene_name)
down_degs_24h_list = unique(down_degs_24h$external_gene_name)
#Up regulated stats
up_shared_genes = intersect(up_degs_24h_list, up_degs_2h_list)
up_unique_24h = setdiff(up_degs_24h_list, up_degs_2h_list)
up_unique_2h = setdiff(up_degs_2h_list, up_degs_24h_list)
#Down stats
down_shared_genes = intersect(down_degs_24h_list, down_degs_2h_list)
down_unique_24h = setdiff(down_degs_24h_list, down_degs_2h_list)
down_unique_2h = setdiff(down_degs_2h_list, down_degs_24h_list)
# Print statistics - up reg
cat("Upregulated DEGs unique to 24h:", length(up_unique_24h), "\n")
cat("Upregulated DEGs unique to 2h:", length(up_unique_2h), "\n")
cat("Upregulated DEGs shared:", length(up_shared_genes), "\n")
cat("Total DEGs Upregulated at 24h:", length(up_degs_24h_list), "\n")
cat("Total DEGs upregulated at 2h:", length(up_degs_2h_list),"\n")
# Print statistics - down reg
cat("Downregulated DEGs unique to 24h:", length(down_unique_24h), "\n")
cat("Downregulated DEGs unique to 2h:", length(down_unique_2h), "\n")
cat("Downregulated DEGs shared:", length(down_shared_genes), "\n")
cat("Total DEGs Downregulated at 24h:", length(down_degs_24h_list), "\n")
cat("Total DEGs Downregulated at 2h:", length(down_degs_2h_list),"\n")
#Create upregulated venn diagram
up_venn = ggVennDiagram(list("24h" = up_degs_24h_list, "2h" = up_degs_2h_list)) +
  scale_fill_gradient(low = "white", high = "#9d0208", name = "Counts") +
  ggtitle("Upregulated DEGs Overlap at 2h and 24h") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    text = element_text(size = 6),
    plot.margin = unit(c(2, 2, 2, 2), "cm")
  )
ggsave("results/figures/upregulated_venndiagram.pdf",
       up_venn)
#Create downregulated venn diagram
down_venn = ggVennDiagram(list("24h" = down_degs_24h_list, "2h" = down_degs_2h_list)) +
  scale_fill_gradient(low = "white", high = "#003049", name = "Counts") +
  ggtitle("Downregulated DEGs Overlap at 2h and 24h") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    text = element_text(size = 6),
    plot.margin = unit(c(2, 2, 2, 2), "cm")
  )
ggsave("results/figures/downregulated_venndiagram.pdf",
       down_venn)
#Combine venn diagrams
combined_venn = up_venn + down_venn +
  plot_annotation(
    title = "Figure 5: Overlap of upregulated and downregulated DEGs ",
    tag_levels = "A"
  ) & #'&' applys theme to all sub plots
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave("results/figures/figure5_combined_venn.pdf",
       combined_venn,
       width = 12, height = 6)

#Heatmap of DEGs - 2hrs ----
#get 250 largest fold changes that are also significant
degs_2h = mutate(degs_2h, abslog2FoldChange=abs(log2FoldChange))
degs_2h = arrange(degs_2h, desc(abslog2FoldChange))
degs250_2h = head(degs_2h, 250)
#create heatmap
#creat matrix of top genes x all samples - 
matrix_degs250_2h = assay(rld)[degs250_2h$ensembl_gene_id,]
#create anotations
degs_heatmap_annotation = data.frame(
  group=colData(dds)[,c('Group')],
  row.names=rownames(colData(dds))
)
colnames(degs_heatmap_annotation)[colnames(degs_heatmap_annotation) == "group"] <- "Timepoint"
#create and save heat map of hierarchical clusering of the top 250 degs that arre significant, using rlog transformed counts
degs_2h_heatmap = pheatmap(matrix_degs250_2h,
                           scale='row',
                           clustering_distance_rows='euclidean',
                           clustering_distance_cols='euclidean',
                           clustering_method = 'complete',
                           color =heatmap_colours,
                           annotation_col = degs_heatmap_annotation,
                           annotation_colors = annotation_colors,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           annotation_names_col = FALSE,
                           annotation_names_row = FALSE,
                           fontsize = 10,
                           border_color = "black",
                           na_col = "white",
                           main = "Heatmap of the Top 250 Differentially Expressed Genes (Allo2h vs Naïve)",
                           filename = "results/figures/degs_2h_heatmap.pdf")
#Heatmap of DEGs - 24hrs ----
#get 250 largest fold changes that are also significant
degs_24h = mutate(degs_24h, abslog2FoldChange=abs(log2FoldChange))
degs_24h = arrange(degs_24h, desc(abslog2FoldChange))
degs250_24h = head(degs_24h, 250)
#create matrix of top genes x all samples - 
matrix_degs250_24h = assay(rld)[degs250_24h$ensembl_gene_id,]
#create and save heat map of hierarchical clusering of the top 250 degs that are significant, using rlog transformed counts
degs_24h_heatmap = pheatmap(matrix_degs250_24h,
                            scale='row',
                            clustering_distance_rows='euclidean',
                            clustering_distance_cols='euclidean',
                            clustering_method = 'complete',
                            color =heatmap_colours,
                            annotation_col = degs_heatmap_annotation,
                            annotation_colors = annotation_colors,
                            show_rownames = FALSE,
                            show_colnames = FALSE,
                            annotation_names_col = FALSE,
                            annotation_names_row = FALSE,
                            fontsize = 10,
                            border_color = "black",
                            na_col = "white",
                            main = "Heatmap of the Top 250 Differentially Expressed Genes (Allo24h vs Naïve)",
                            filename = "results/figures/degs_24h_heatmap.pdf")
#Multi-panel figure of DEGs heatmaps ----
#pheatmap creats a list of componants, $gtable i the graphical object that r ues to mke the heatmap
combined_degs_heatmaps =
  wrap_elements(degs_2h_heatmap$gtable) /
  wrap_elements(degs_24h_heatmap$gtable) +
  plot_annotation(tag_levels = "A")
ggsave(
  "results/figures/figure6_degs_heatmaps_combined.pdf",
  combined_degs_heatmaps,
  width = 8,
  height = 10
)