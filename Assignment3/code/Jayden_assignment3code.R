# Install necessary packages if not already installed
if (!require("readr")) install.packages("readr")
if (!require("jsonlite")) install.packages("jsonlite")
if (!require("cluster")) install.packages("cluster")
if (!require("ggplot2")) install.packages("ggplot2")

# Load necessary libraries
library(readr)
library(jsonlite)
library(cluster)
library(ggplot2)

# Load expression data and metadata
expression_data <- read_tsv("expression_data.tsv")
metadata <- fromJSON("aggregated_metadata.json")

# Function to run silhouette analysis for different k values
silhouette_analysis <- function(expression_data, num_genes, max_k = 10) {
  gene_variances <- apply(expression_data[, -1], 1, var)
  top_genes <- order(gene_variances, decreasing = TRUE)[1:num_genes]
  subset_expression_data <- expression_data[top_genes, ]
  
  distance_matrix <- dist(subset_expression_data[, -1], method = "euclidean")
  
  avg_sil_width <- c()  # Store average silhouette width for each k
  
  # Loop through different values of k to calculate silhouette width
  for (k in 2:max_k) {
    clusters <- cutree(hclust(distance_matrix, method = "complete"), k = k)
    sil <- silhouette(clusters, distance_matrix)
    avg_sil_width[k] <- mean(sil[, 3])  # Store average silhouette score for each k
  }
  
  # Plot silhouette scores for each k
  png(filename = paste0("silhouette_analysis_", num_genes, "_genes.png"))
  plot(2:max_k, avg_sil_width[2:max_k], type = "b", xlab = "Number of Clusters (k)", 
       ylab = "Average Silhouette Width", main = paste("Silhouette Analysis for", num_genes, "Genes"))
  dev.off()
  
  # Return the best k based on the highest silhouette score
  best_k <- which.max(avg_sil_width)
  return(best_k)
}

# Function to run gap statistic for different k values
gap_statistic_analysis <- function(expression_data, num_genes, max_k = 10) {
  gene_variances <- apply(expression_data[, -1], 1, var)
  top_genes <- order(gene_variances, decreasing = TRUE)[1:num_genes]
  subset_expression_data <- expression_data[top_genes, ]
  
  # Ensure that the data is a matrix for k-means
  data_matrix <- as.matrix(subset_expression_data[, -1])  # Ensure we don't use the gene names column
  
  # K-means clustering function for clusGap that returns only cluster labels
  kmeans_clustering <- function(X, k) {
    kmeans(X, centers = k, nstart = 25)$cluster
  }
  
  # Compute the gap statistic for k from 1 to max_k using k-means clustering
  set.seed(123)  # Ensure reproducibility
  gap_stat <- clusGap(data_matrix, FUNcluster = kmeans_clustering, K.max = max_k, B = 50)
  
  # Plot the gap statistic
  png(filename = paste0("gap_statistic_", num_genes, "_genes.png"))
  plot(gap_stat, main = paste("Gap Statistic for", num_genes, "Genes (K-means)"))
  dev.off()
  
  # Return the optimal number of clusters (k) based on the gap statistic
  best_k <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "Tibs2001SEmax")
  return(best_k)
}

# Run gap statistic for different gene sets
best_k_100_gap <- gap_statistic_analysis(expression_data, 100)
best_k_1000_gap <- gap_statistic_analysis(expression_data, 1000)

# Print gap statistic results
print(paste("Best k based on Gap Statistic for 100 genes (K-means):", best_k_100_gap))
print(paste("Best k based on Gap Statistic for 1000 genes (K-means):", best_k_1000_gap))


# Run gap statistic for different gene sets
best_k_100_gap <- gap_statistic_analysis(expression_data, 100)
best_k_1000_gap <- gap_statistic_analysis(expression_data, 1000)

# Print gap statistic results
print(paste("Best k based on Gap Statistic for 100 genes:", best_k_100_gap))
print(paste("Best k based on Gap Statistic for 1000 genes:", best_k_1000_gap))

###########################################

# Function to perform hierarchical clustering for a subset of data with a given number of genes
run_clustering <- function(expression_data, num_genes) {
  gene_variances <- apply(expression_data[, -1], 1, var)
  top_genes <- order(gene_variances, decreasing = TRUE)[1:num_genes]
  subset_expression_data <- expression_data[top_genes, ]
  
  # Calculate distance matrix and perform clustering
  distance_matrix <- dist(subset_expression_data[, -1], method = "euclidean")
  hc <- hclust(distance_matrix, method = "complete")
  
  # Cut the dendrogram into a fixed number of clusters (e.g., k = 3)
  clusters <- cutree(hc, k = 3)
  
  # Return the cluster assignments
  return(clusters)
}

# Perform clustering for different numbers of genes
clusters_10 <- run_clustering(expression_data, 10)
clusters_100 <- run_clustering(expression_data, 100)
clusters_1000 <- run_clustering(expression_data, 1000)
clusters_10000 <- run_clustering(expression_data, 10000)

# Now we have clustering results for 10, 100, 1000, and 10000 genes

# Chi-Squared Test Function
perform_chi_squared_test <- function(clusters_1, clusters_2) {
  contingency_table <- table(clusters_1, clusters_2)
  chi_sq_result <- chisq.test(contingency_table)
  return(chi_sq_result)
}

# Identify common samples across all clustering results
common_samples <- Reduce(intersect, list(rownames(expression_data), names(clusters_10), names(clusters_100), names(clusters_1000), names(clusters_10000)))

# Subset each clustering result to include only common samples
clusters_10_common <- clusters_10[common_samples]
clusters_100_common <- clusters_100[common_samples]
clusters_1000_common <- clusters_1000[common_samples]
clusters_10000_common <- clusters_10000[common_samples]

# Perform chi-squared tests using only the common samples
chi_sq_10_100 <- perform_chi_squared_test(clusters_10_common, clusters_100_common)
chi_sq_100_1000 <- perform_chi_squared_test(clusters_100_common, clusters_1000_common)
chi_sq_1000_10000 <- perform_chi_squared_test(clusters_1000_common, clusters_10000_common)

# Print the results
chi_sq_10_100
chi_sq_100_1000
chi_sq_1000_10000


# Create a summary table for chi-squared test results
chi_sq_table <- data.frame(
  Comparison = c("10 vs 100", "100 vs 1000", "1000 vs 10000"),
  Chi_Squared_Value = c(chi_sq_10_100$statistic, chi_sq_100_1000$statistic, chi_sq_1000_10000$statistic),
  p_value = c(chi_sq_10_100$p.value, chi_sq_100_1000$p.value, chi_sq_1000_10000$p.value)
)

# Output the summary table
print(chi_sq_table)

if (!require("ggalluvial")) install.packages("ggalluvial")
library(ggalluvial)

# Use reshape2 to reshape the dataframe for plotting
library(reshape2)
df_alluvial <- data.frame(
  sample = common_samples,  # Use common samples
  `10_genes` = as.factor(clusters_10_common),
  `100_genes` = as.factor(clusters_100_common),
  `1000_genes` = as.factor(clusters_1000_common),
  `10000_genes` = as.factor(clusters_10000_common)
)

# Reshape the data into a long format for the alluvial diagram
df_alluvial_melted <- melt(df_alluvial, id.vars = "sample", variable.name = "Gene_Set", value.name = "Cluster")

# Create the alluvial diagram
alluvial_plot <- ggplot(df_alluvial_melted,
       aes(x = Gene_Set, stratum = Cluster, alluvium = sample, fill = Cluster)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum() +
  theme_minimal() +
  labs(title = "Alluvial Diagram of Cluster Membership Across Gene Sets",
       x = "Gene Set", y = "Samples") +
  theme(legend.position = "bottom")

ggsave("alluvial_diagram_cluster_membership.png", plot = alluvial_plot, width = 10, height = 7)

#######
# CSV File Upload

# Assuming you have your clustering results for 10, 100, 1000, or 10000 genes (e.g., clusters_10_common)

# Create a dataframe with sample IDs and the corresponding cluster
# Limiting to 71 rows for the first 71 samples

# You can choose one of your clustering results (e.g., clusters_10_common)
cluster_results <- clusters_100[1:71]  # Limit to the first 71 samples

cluster_results

# Assuming your sample IDs are stored in the rownames of your data
sample_ids <- names(cluster_results)

sample_ids_prefixed <- paste0("ERR", sample_ids)

# Create a dataframe
sample_cluster_df <- data.frame(
  sampleId = sample_ids_prefixed,
  cluster = cluster_results
)

# Save the dataframe to a CSV file
write.csv(sample_cluster_df, "sample_cluster_assignment.csv", row.names = FALSE)

#######################
# Combined Heat Map

# Load necessary libraries
library(pheatmap)

# File paths (adjust these to your local file paths)
files <- list(
  "grant, a. spectral_results.csv",
  "kmeans_results.csv",
  "pam_clustering_results_k2.csv",
  "propagation.csv",
  "spitulnik,hcluster_results.csv"
)

# Initialize an empty list to store the cluster assignments from each file
clusters_list <- list()

# Loop through each file and load the clustering results
for (file in files) {
  # Read the CSV file
  cluster_data <- read.csv(file, header = TRUE)
  
  # Assuming the file has two columns: sampleId and cluster, we'll extract the cluster information
  clusters_list[[file]] <- cluster_data$cluster  # Store only the cluster column
}

# Combine the clustering results into a matrix
# Each row is a sample, and each column is a clustering method
clusters_matrix <- do.call(cbind, clusters_list)

# Set the sample IDs as rownames (assuming the first column of each file is 'sampleId')
rownames(clusters_matrix) <- cluster_data$sampleId

# Generate the heatmap and save it as a PNG file
png("clustering_heatmap.png", width = 10, height = 8, units = "in", res = 300)

pheatmap(
  clusters_matrix,
  scale = "none",  # No scaling of values since we are dealing with cluster numbers
  clustering_method = "complete",  # Complete linkage hierarchical clustering
  cluster_rows = TRUE,  # Add row dendrogram
  cluster_cols = TRUE,  # Add column dendrogram
  legend = TRUE,  # Add a legend
  main = "Heatmap of Clustering Results",  # Title for the heatmap
  show_rownames = TRUE,  # Show sample IDs as row names
  show_colnames = TRUE,  # Show the clustering methods as column names
  fontsize_row = 6,  # Adjust font size for row names (sample IDs)
  fontsize_col = 10  # Adjust font size for column names (clustering methods)
)

dev.off()  # Close the PNG device and save the file


