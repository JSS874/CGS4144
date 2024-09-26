#Part 1

#Set Directory
setwd("./Dataset")

#Installs
install.packages("readr")
install.packages("dplyr")
install.packages("ggplot2")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

#Libraries
library(readr)
library(dplyr)
library(ggplot2)
library(DESeq2)

#Read tsv files
counts <- read_tsv("ERP105973/ERP105973.tsv")
metadata <- read_tsv("ERP105973/metadata_ERP105973.tsv")

# Set row names to gene identifiers
rownames(counts) <- counts$Gene
counts <- counts %>% select(-Gene)

# Step 1c: Check the size of the matrix
num_genes <- nrow(counts)
num_samples <- ncol(counts)

cat("Number of genes:", num_genes, "\n")
cat("Number of samples:", num_samples, "\n")

# Step 1c: Log-scale the data
log_counts <- log2(counts + 1)  # Adding 1 to avoid log of zero

# Step 1c: Calculate per-gene median expression ranges
median_expr <- apply(log_counts, 1, median)

# Step 1c: Create a density plot
ggplot(data.frame(median_expr), aes(x = median_expr)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density of Per-Gene Median Expression (Log-scaled)",
       x = "Median Expression (Log2)",
       y = "Density")

# Install biomaRt if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

# Load the biomaRt package
library(biomaRt)

# Connect to the Ensembl database for Sus scrofa
ensembl <- useMart("ensembl", dataset = "sscrofa_gene_ensembl")

# Extract the gene list (Ensembl IDs) from the counts matrix
ensembl_ids <- rownames(counts)

# Query biomaRt to convert Ensembl IDs to Hugo gene names
gene_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                         filters = "ensembl_gene_id",
                         values = ensembl_ids,
                         mart = ensembl)

# Merge the results with your counts data
# Replace Ensembl IDs with Hugo gene names where possible
gene_conversion <- gene_conversion %>%
  filter(hgnc_symbol != "")  # Remove rows where no conversion was found

# Use Hugo gene names for the genes that were successfully converted
counts_hugo <- counts
rownames(counts_hugo)[rownames(counts_hugo) %in% gene_conversion$ensembl_gene_id] <- 
  gene_conversion$hgnc_symbol

# End of Step 1

# Load necessary libraries
install.packages("Rtsne")
install.packages("umap")
library(umap)
library(Rtsne)
library(DESeq2)



summary(counts_hugo)

# Remove genes with zero variance across samples
counts_filtered <- counts_hugo[apply(counts_hugo, 1, var) > 0, ]

# Perform PCA on the filtered data
pca_result <- prcomp(t(counts_filtered), scale. = TRUE)

# Prepare PCA data for plotting
pca_data <- as.data.frame(pca_result$x)
pca_data$Treatment <- metadata$refinebio_treatment  # Add treatment information

# Plot PCA
ggplot(pca_data, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1 (", round(100 * summary(pca_result)$importance[2, 1], 1), "% variance)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result)$importance[2, 2], 1), "% variance)")) +
  ggtitle("PCA Plot") +
  theme_minimal() +
  labs(color="Treatment")

# t-SNE
tsne_result <- Rtsne(t(counts_filtered), dims = 2, perplexity = 20)
tsne_data <- data.frame(tsne_result$Y, Treatment = metadata$refinebio_treatment)

ggplot(tsne_data, aes(X1, X2, color=Treatment)) +
  geom_point(size=3) +
  ggtitle("t-SNE Plot") +
  xlab("t-SNE Dimension 1") +
  ylab("t-SNE Dimension 2") +
  theme_minimal() +
  labs(color="Treatment")

# UMAP
umap_result <- umap(t(counts_filtered))
umap_data <- data.frame(umap_result$layout, Treatment = metadata$refinebio_treatment)

ggplot(umap_data, aes(X1, X2, color=Treatment)) +
  geom_point(size=3) +
  ggtitle("UMAP Plot") +
  xlab("UMAP Dimension 1") +
  ylab("UMAP Dimension 2") +
  theme_minimal() +
  labs(color="Treatment")

# Save the PCA plot
ggsave("PCA_plot.png", width=8, height=6)

# Save the t-SNE plot
ggsave("tSNE_plot.png", width=8, height=6)

# Save the UMAP plot
ggsave("UMAP_plot.png", width=8, height=6)

# End of Step 2
# Start of Step 3

# Install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

# Load required packages
library(limma)

# Prepare the design matrix for the two groups
group <- factor(metadata$refinebio_treatment)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit the linear model
fit <- lmFit(counts_filtered, design)

# Create a contrast for barren vs enriched
contrast_matrix <- makeContrasts(barren_vs_enriched = barren - enriched, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2, trend = TRUE)

# Get the top differentially expressed genes
top_genes <- topTable(fit2, adjust = "fdr", number = 50)

# Save the full results
full_results <- topTable(fit2, adjust = "fdr", number = Inf)

# Create a volcano plot of the results
ggplot(full_results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("red", "grey")) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Volcano Plot of Differentially Expressed Genes") +
  theme_minimal() +
  labs(color = "Adjusted p-value < 0.05")

# Save the plot to a file
ggsave("volcano_plot.png", width=8, height=6)

# Save the top 50 genes to a CSV file
write.csv(top_genes, file = "top_50_genes.csv")

# Save the full results to a CSV file
write.csv(full_results, file = "full_results.csv")

# End of Step 3
# Start of Step 4

# Install ComplexHeatmap if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

# Load required packages
library(ComplexHeatmap)
library(circlize)  # Required for color annotations

# Extract significantly differentially expressed genes (adjusted p-value < 0.05)
signif_genes <- full_results[full_results$adj.P.Val < 0.05, ]

# Subset the count matrix to include only the significant genes
counts_signif <- counts_filtered[rownames(counts_filtered) %in% rownames(signif_genes), ]

# Create a color annotation for the groups (enriched vs barren)
sample_groups <- metadata$refinebio_treatment

# Define colors for the annotation
group_colors <- list(refinebio_treatment = c("barren" = "blue", "enriched" = "red"))

# Create the annotation object
ha <- HeatmapAnnotation(Group = sample_groups, col = group_colors)

# Extract the top 30 significant genes
top_signif_genes <- head(signif_genes[order(signif_genes$adj.P.Val), ], 30)

# Subset the count matrix to include only the top significant genes
counts_top_signif <- counts_filtered[rownames(counts_filtered) %in% rownames(top_signif_genes), ]


# Generate the heatmap for significantly differentially expressed genes
Heatmap(as.matrix(counts_signif),
        name = "Expression",
        top_annotation = ha,  # Add sample group annotations
        show_row_names = FALSE,  # Hide row (gene) names for clarity
        show_column_names = TRUE,  # Show sample names
        clustering_method_rows = "complete",  # Hierarchical clustering for rows
        clustering_method_columns = "complete",  # Hierarchical clustering for columns
        row_title = "Significantly DE Genes",
        column_title = "Samples",
        col = colorRamp2(c(min(counts_signif), 0, max(counts_signif)), c("blue", "white", "red"))  # Color scale
)

# Increase plot size and display top 30 genes with row labels
pdf("heatmap_top_significant_genes.pdf", width = 12, height = 10)  # Increase figure size
Heatmap(as.matrix(counts_top_signif),
        name = "Expression",
        top_annotation = ha,  # Sample group annotation
        show_row_names = TRUE,  # Show gene names (rows)
        show_column_names = TRUE,  # Show sample names
        clustering_method_rows = "complete",  # Hierarchical clustering for rows
        clustering_method_columns = "complete",  # Hierarchical clustering for columns
        row_title = "Top 30 Significantly DE Genes",
        column_title = "Samples",
        col = colorRamp2(c(min(counts_top_signif), 0, max(counts_top_signif)), c("blue", "white", "red"))
)
dev.off()

# End of Step 4
# Start of Step 5

# Install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")

# Load required libraries
library(limma)
library(edgeR)

# Assume 'counts_filtered' is the filtered count matrix, and 'metadata' contains the treatment information
group <- factor(metadata$refinebio_treatment)  # Adjust this based on your metadata column
design <- model.matrix(~0 + group)  # Create the design matrix without intercept
colnames(design) <- levels(group)

# Replace negative values with zero
counts_filtered[counts_filtered < 0] <- 0

# Create a DGEList object (for voom normalization)
dge <- DGEList(counts = counts_filtered)

# Normalize the data using the voom function (log-transformed)
v <- voom(dge, design, plot = TRUE)

# Fit the linear model using limma
fit <- lmFit(v, design)

# Create contrast matrix for comparison (e.g., barren vs enriched)
contrast_matrix <- makeContrasts(barren_vs_enriched = barren - enriched, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)

# Apply empirical Bayes moderation to the fitted model
fit2 <- eBayes(fit2)

# Extract the top differentially expressed genes
top_genes <- topTable(fit2, adjust = "fdr", number = 50)

# Save the top 50 genes to a CSV file
write.csv(top_genes, file = "top_50_genes_limma.csv")

# Create a volcano plot
volcano <- topTable(fit2, adjust = "fdr", number = Inf)
ggplot(volcano, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("red", "grey")) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Volcano Plot of Differentially Expressed Genes") +
  theme_minimal()

# Save the volcano plot
ggsave("volcano_plot_limma.png", width = 8, height = 6)

# Assuming 'top_genes' contains the results with Ensembl IDs
sig_genes <- top_genes
sig_genes$ensembl_gene_id <- rownames(sig_genes)  # Add Ensembl gene IDs

# Install biomaRt if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

# Load biomaRt library
library(biomaRt)

# Connect to Ensembl database for Sus scrofa
ensembl <- useMart("ensembl", dataset = "sscrofa_gene_ensembl")

# Extract the gene list (Ensembl IDs) from your filtered count matrix
ensembl_ids <- rownames(counts_filtered)

# Query biomaRt to map Ensembl IDs to Entrez IDs
final_mapped_df <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                         filters = "ensembl_gene_id",
                         values = ensembl_ids,
                         mart = ensembl)

# Save the mapping to a file for future use
write.table(final_mapped_df, file = "final_mapped_df.txt", sep = "\t", row.names = FALSE, quote = FALSE)



# Load the Ensembl to Entrez ID mapping (ensure the file exists)
final_mapped_df <- read.table("final_mapped_df.txt", header = TRUE)

# Merge significant genes with the mapping file to get Entrez IDs
merged_genes <- merge(as.data.frame(sig_genes), final_mapped_df, by = "ensembl_gene_id", all.x = TRUE)

# Extract the Entrez Gene IDs from the merged data
entrez_gene_list <- na.omit(unique(merged_genes$entrezgene_id))

# Perform GO enrichment analysis
go_results <- enrichGO(gene = entrez_gene_list, 
                       OrgDb = org.Ss.eg.db, 
                       keyType = "ENTREZID", 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.2, 
                       readable = TRUE)

# Save the GO enrichment results
go_results_df <- as.data.frame(go_results)
write.table(go_results_df, file = "GO_enrichment_results_limma.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# End of Step 5
# Start of Step 6

# Load necessary libraries
library(dplyr)

# Load results from each method (replace with actual filenames)
clusterProfiler_results <- read.csv("clusterProfiler_results.csv")
gProfiler_results <- read.csv("gProfiler_results.csv")
wilcoxon_results <- read.csv("wilcoxon_results.csv")

# Ensure all datasets have consistent column names (modify as needed)
colnames(clusterProfiler_results) <- c("Term", "pvalue", "qvalue", "logFC", "Method")
colnames(gProfiler_results) <- c("Term", "pvalue", "qvalue", "logFC", "Method")
colnames(wilcoxon_results) <- c("Term", "pvalue", "qvalue", "logFC", "Method")

# Merge the datasets on the "Term" column, keeping all terms (even if they appear in only one method)
combined_results <- full_join(clusterProfiler_results, gProfiler_results, by = "Term", suffix = c("_cluster", "_gprofiler"))
combined_results <- full_join(combined_results, wilcoxon_results, by = "Term", suffix = c("", "_wilcoxon"))

# Replace missing values (e.g., terms that weren't found in some methods)
combined_results[is.na(combined_results)] <- "NA"

# Add a column indicating how many methods found each term as significantly enriched
combined_results$Num_Methods_Significant <- rowSums(!is.na(combined_results[, c("pvalue_cluster", "pvalue_gprofiler", "pvalue_wilcoxon")]) & combined_results[, c("pvalue_cluster", "pvalue_gprofiler", "pvalue_wilcoxon")] < 0.05)

# Add a column indicating how many methods included each term in the analysis
combined_results$Num_Methods_Included <- rowSums(!is.na(combined_results[, c("pvalue_cluster", "pvalue_gprofiler", "pvalue_wilcoxon")]))

# Save the combined results to a CSV file
write.csv(combined_results, file = "combined_enrichment_results.csv", row.names = FALSE)



