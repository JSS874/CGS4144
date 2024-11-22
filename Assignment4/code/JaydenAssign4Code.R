# Load necessary libraries
library(data.table)
library(glmnet)
library(caret)      # For confusion matrix calculation
library(ggplot2)    # For plotting
library(gridExtra)  # For displaying tables as plots
library(grid)       # For saving table grobs as images

# Prepare paths
data_dir <- file.path("data", "ERP105973")
data_file <- file.path(data_dir, "ERP105973.tsv")
metadata_file <- file.path(data_dir, "metadata_ERP105973.tsv")
cluster_results_file <- file.path("data", "spitulnik,hcluster_results.csv")
results_file <- file.path("results", "cluster_prediction_results.csv")

# Load the data
metadata <- fread(metadata_file)
expression_data <- fread(data_file)
cluster_results <- fread(cluster_results_file)

# Prepare the data by selecting top 5,000 most variable genes
gene_variances <- apply(expression_data[, -1, with = FALSE], 1, var)
top_genes <- order(gene_variances, decreasing = TRUE)[1:5000]
expression_data_subset <- expression_data[top_genes, ]

# Transpose expression data for samples to be rows and genes as columns
sample_data <- t(expression_data_subset[, -1, with = FALSE])
colnames(sample_data) <- expression_data_subset$V1
sample_data <- data.table(sample_id = colnames(expression_data_subset)[-1], sample_data)

# Merge metadata and clustering results
sample_data <- merge(sample_data, metadata[, .(refinebio_accession_code)], 
                     by.x = "sample_id", by.y = "refinebio_accession_code")
sample_data <- merge(sample_data, cluster_results, by.x = "sample_id", by.y = "sampleId")
sample_data[, cluster := as.factor(cluster)]

# Run one-vs-all logistic regression for each cluster
clusters <- levels(sample_data$cluster)
predictions <- matrix(NA, nrow = nrow(sample_data), ncol = length(clusters))
colnames(predictions) <- clusters

for (cluster_label in clusters) {
  # Create binary target for current cluster (1 for current cluster, 0 for all others)
  sample_data[, binary_cluster := ifelse(cluster == cluster_label, 1, 0)]
  
  # Run logistic regression using glmnet (ridge regression for stability)
  x <- as.matrix(sample_data[, -c(1, ncol(sample_data), ncol(sample_data)-1), with = FALSE])
  y <- as.factor(sample_data$binary_cluster)
  fit <- glmnet(x, y, family = "binomial", lambda = 0.01)
  
  # Get predictions as probabilities for each sample
  predictions[, cluster_label] <- predict(fit, x, type = "response")[, 1]
}

# Assign each sample to the cluster with the highest probability
sample_data[, predicted_cluster := apply(predictions, 1, function(row) {
  clusters[which.max(row)]
})]

# Save predictions to CSV in specified format
results <- sample_data[, .(sampleId = sample_id, predicted_class = predicted_cluster)]
fwrite(results, results_file)

# Calculate Confusion Matrix and Accuracy
conf_matrix <- table(sample_data$cluster, sample_data$predicted_cluster)
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)

# Save Confusion Matrix to CSV
conf_matrix_df <- as.data.frame(as.table(conf_matrix))
colnames(conf_matrix_df) <- c("Actual", "Predicted", "Freq")
write.csv(conf_matrix_df, "results/confusion_matrix.csv", row.names = FALSE)

# Save Accuracy to CSV
accuracy_df <- data.frame(Metric = "Accuracy", Value = round(accuracy, 2))
write.csv(accuracy_df, "results/accuracy.csv", row.names = FALSE)

# Plot Confusion Matrix
conf_matrix_plot <- ggplot(conf_matrix_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "black") +
  geom_text(aes(label = Freq), color = "white", size = 5) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Confusion Matrix", x = "Predicted Class", y = "Actual Class") +
  theme_minimal()

# Save Confusion Matrix Plot as PNG
ggsave("results/confusion_matrix_plot.png", plot = conf_matrix_plot, width = 6, height = 5)

# Display and save Accuracy Table as PNG
accuracy_table <- tableGrob(accuracy_df)
png("results/accuracy_table.png", width = 400, height = 100)
grid.draw(accuracy_table)
dev.off()

# Display Confusion Matrix Plot and Accuracy Table together (optional, for viewing)
grid.arrange(conf_matrix_plot, accuracy_table, ncol = 1, heights = c(2, 0.5))

#

# Calculate Confusion Matrix and Accuracy
conf_matrix <- table(sample_data$predicted_cluster, sample_data$cluster)
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)

# Save Confusion Matrix in long format (response, truth, Freq)
conf_matrix_long <- as.data.frame(as.table(conf_matrix))
colnames(conf_matrix_long) <- c("response", "truth", "Freq")
write.csv(conf_matrix_long, "results/confusion_matrix_long.csv", row.names = FALSE)

# Save Accuracy to CSV
accuracy_df <- data.frame(Metric = "Accuracy", Value = round(accuracy, 2))
write.csv(accuracy_df, "results/accuracy.csv", row.names = FALSE)

# Plot Confusion Matrix
conf_matrix_plot <- ggplot(conf_matrix_long, aes(x = response, y = truth, fill = Freq)) +
  geom_tile(color = "black") +
  geom_text(aes(label = Freq), color = "white", size = 5) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Confusion Matrix", x = "Predicted Class", y = "Actual Class") +
  theme_minimal()

# Save Confusion Matrix Plot as PNG
ggsave("results/confusion_matrix_plot.png", plot = conf_matrix_plot, width = 6, height = 5)

# Display and save Accuracy Table as PNG
accuracy_table <- tableGrob(accuracy_df)
png("results/accuracy_table.png", width = 400, height = 100)
grid.draw(accuracy_table)
dev.off()

# Display Confusion Matrix Plot and Accuracy Table together (optional, for viewing)
grid.arrange(conf_matrix_plot, accuracy_table, ncol = 1, heights = c(2, 0.5))



### Step 3

# Load necessary libraries
library(data.table)
library(dplyr)

# Load each prediction file
files <- list(
  KNN = fread("data/KNN_prediction_results.csv"),
  RandomForest = fread("data/RandomForest_combined_results.csv"),
  LinearRegression = fread("data/linear_regression_prediction_results.csv"),
  SVM = fread("data/grant, a. svm_predictions.csv"),
  NaiveBayes = fread("data/NaiveBayesGroups.csv")
)

# Rename the `PredictedClass` column in each file with the model name to avoid duplication
for (name in names(files)) {
  setnames(files[[name]], "PredictedClass", name)
}

# Merge all files by "SampleID"
predictions <- Reduce(function(x, y) merge(x, y, by = "SampleID", all = TRUE), files)

# Step 3a: Count predictions for each class label by sample
# Creating columns Class1, Class2, and Class3
class_counts <- predictions %>%
  pivot_longer(cols = -SampleID, names_to = "model", values_to = "PredictedClass") %>%
  group_by(SampleID, PredictedClass) %>%
  summarize(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = PredictedClass, values_from = count, values_fill = 0) %>%
  rename(Class1 = `1`, Class2 = `2`, Class3 = `3`)

# Step 3b: Count predictions for the same cluster by sample
# Calculate the number of models that predicted the same cluster for each sample
predictions$SameCluster <- apply(predictions[, -1, with = FALSE], 1, function(row) {
  max(table(row))  # Count the maximum occurrence of any predicted class
})

# Combine the tables into the final format
final_table <- predictions %>%
  left_join(class_counts, by = "SampleID") %>%
  select(SampleID, KNN, RandomForest, LinearRegression, SVM, NaiveBayes, Class1, Class2, Class3, SameCluster)

# Save the final table to a CSV file
write.csv(final_table, "results/matrix_of_samples_by_model.csv", row.names = FALSE)

#

# Load necessary libraries
library(data.table)
library(dplyr)
library(stats)

# Load matrix data
matrix_data <- fread("results/matrix_of_samples_by_model.csv")

### Part (a): Count Models Predicting Each Class Label
# Assuming columns `Class1`, `Class2`, `Class3` contain counts of models predicting each class
class_counts <- matrix_data %>%
  select(SampleID, Class1, Class2, Class3)

# Save counts for each class label to a CSV file
write.csv(class_counts, "results/class_counts_per_sample.csv", row.names = FALSE)

### Part (b): Count Models Predicting the Same Cluster
# `SameCluster` column contains the count of models that agreed on the same cluster for each sample
same_cluster_counts <- matrix_data %>%
  select(SampleID, SameCluster)

# Save counts of same cluster predictions to a CSV file
write.csv(same_cluster_counts, "results/same_cluster_counts_per_sample.csv", row.names = FALSE)

### Part (c): Correlation between Cluster Stability and Class Label Stability
# Calculate class label stability as the maximum of Class1, Class2, or Class3 for each sample
matrix_data <- matrix_data %>%
  mutate(label_stability = pmax(Class1, Class2, Class3))

# Correlation test between `label_stability` and `SameCluster`
cor_test <- cor.test(matrix_data$label_stability, matrix_data$SameCluster, method = "spearman")

# Multiple testing correction using Benjamini-Hochberg
# Here, only one p-value exists, but if there were multiple tests, you'd use p.adjust()
p_adjusted <- p.adjust(cor_test$p.value, method = "BH")

# Save correlation results to a CSV file
cor_results <- data.frame(
  Correlation_Coefficient = cor_test$estimate,
  P_Value = cor_test$p.value,
  Adjusted_P_Value = p_adjusted,
  Method = cor_test$method
)
write.csv(cor_results, "results/correlation_test_results.csv", row.names = FALSE)



#### Step 4

# Define the number of top variable genes to try
gene_counts <- c(10, 100, 1000, 10000)

# Prepare a data.table to store all results in a single file
final_results <- data.table(sampleId = sample_data$sample_id)  # Initialize with sample IDs

# Loop over each gene count, retrain, and add predictions as new columns
for (n_genes in gene_counts) {
  # Select the top n_genes most variable genes
  top_genes <- order(gene_variances, decreasing = TRUE)[1:n_genes]
  expression_data_subset <- expression_data[top_genes, ]
  
  # Transpose expression data for samples to be rows and genes as columns
  sample_data <- t(expression_data_subset[, -1, with = FALSE])
  colnames(sample_data) <- expression_data_subset$V1
  sample_data <- data.table(sample_id = colnames(expression_data_subset)[-1], sample_data)
  
  # Merge metadata and clustering results
  sample_data <- merge(sample_data, metadata[, .(refinebio_accession_code)], 
                       by.x = "sample_id", by.y = "refinebio_accession_code")
  sample_data <- merge(sample_data, cluster_results, by.x = "sample_id", by.y = "sampleId")
  sample_data[, cluster := as.factor(cluster)]
  
  # Run one-vs-all logistic regression with specified lambda
  clusters <- levels(sample_data$cluster)
  predictions <- matrix(NA, nrow = nrow(sample_data), ncol = length(clusters))
  colnames(predictions) <- clusters
  
  for (cluster_label in clusters) {
    # Create binary target for the current cluster (1 for the current cluster, 0 for others)
    sample_data[, binary_cluster := ifelse(cluster == cluster_label, 1, 0)]
    
    # Run logistic regression with a custom lambda for stability
    x <- as.matrix(sample_data[, -c(1, ncol(sample_data), ncol(sample_data)-1), with = FALSE])
    y <- as.factor(sample_data$binary_cluster)
    fit <- glmnet(x, y, family = "binomial", lambda = 0.01)
    
    # Get predictions as probabilities for each sample
    predictions[, cluster_label] <- predict(fit, x, type = "response")[, 1]
  }
  
  # Assign each sample to the cluster with the highest probability
  sample_data[, predicted_cluster := apply(predictions, 1, function(row) {
    clusters[which.max(row)]
  })]
  
  # Add predictions to final_results as a new column
  final_results[, paste0("top_", n_genes, "_genes") := sample_data$predicted_cluster]
}

# Save the combined results to a single CSV file
fwrite(final_results, "results/combined_cluster_predictions.csv")

#### 4b

# Load necessary libraries
library(data.table)
library(pROC)
library(dplyr)

# Load TrueClass file
true_labels <- fread("data/trueclass.csv")  # assumed to contain columns: SampleID, TrueClass

# Define file paths for each model
file_paths <- list(
  knn = "data/combined_knn.csv",
  lr = "data/combined_lr.csv",
  nbg = "data/combined_nbg.csv",
  rf = "data/combined_rf.csv",
  svm = "data/combined_svm.csv"
)

# Initialize a list to store AUC results for each model and gene count
auc_results <- list()

# Loop through each model file and calculate AUC
for (model_name in names(file_paths)) {
  # Load the data for each model
  data <- fread(file_paths[[model_name]])
  
  # Merge the data with true labels on SampleID
  data <- merge(data, true_labels, by = "SampleID")
  
  # Calculate AUC for each gene count version (10, 100, 1000, 10000)
  gene_counts <- c("top_10_genes", "top_100_genes", "top_1000_genes", "top_10000_genes")
  
  for (gene_count in gene_counts) {
    # Calculate AUC for the current gene count
    roc_obj <- roc(data$TrueClass, as.numeric(data[[gene_count]]))
    auc_val <- auc(roc_obj)
    
    # Store the results with model name and gene count
    auc_results[[paste0(model_name, "_", gene_count)]] <- data.frame(
      Model = model_name,
      Gene_Count = gene_count,
      AUC = auc_val
    )
  }
}

# Combine AUC results into a single data frame
combined_auc_results <- do.call(rbind, auc_results)

# Save the combined AUC results to a CSV
write.csv(combined_auc_results, "results/combined_auc_results.csv", row.names = FALSE)

# Analyze AUC trend across gene counts
auc_trend <- combined_auc_results %>%
  dplyr::group_by(Model, Gene_Count) %>%
  dplyr::summarize(Avg_AUC = mean(AUC, na.rm = TRUE))

# Save the AUC trend to a separate CSV file
write.csv(auc_trend, "results/auc_trend_by_gene_count.csv", row.names = FALSE)

# Print out AUC trends to check if performance increases or decreases with gene count
print(auc_trend)

### Step 5

# Load necessary libraries
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Load the data files
trueclass <- read.csv("data/trueclass.csv", row.names = 1)
matrix_data <- read.csv("results/matrix_of_samples_by_model.csv", row.names = 1)

# Extract the sample groups and convert to enriched/barren labels
sample_groups <- factor(trueclass$TrueClass, levels = c(1, 2), labels = c("enriched", "barren"))

# Create an annotation data frame for the sample groups
annotation <- data.frame(SampleGroups = sample_groups)
rownames(annotation) <- rownames(trueclass)

# Define color palette for heatmap
palette_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

# Set the breaks to range from 0 to 4
breaks <- seq(0, 4, length.out = 101)

# Save the heatmap plot as a .png file
png("plots/heatmap_plot.png", width = 1000, height = 1000)
pheatmap(as.matrix(matrix_data),
         annotation_row = annotation,
         color = palette_colors,
         scale = "none",
         breaks = breaks,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_colors = list(SampleGroups = c(enriched = "salmon", barren = "skyblue")))
dev.off()
