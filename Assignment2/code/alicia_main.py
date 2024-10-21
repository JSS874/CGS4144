import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from gprofiler import GProfiler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap

os.makedirs('results', exist_ok=True)

# Download the expression data and matching metadata from Refine.Bio that you selected in Assignment 1
exp_data = pd.read_csv('ERP105973.tsv', sep='\t', index_col=0)

# What size is your expression matrix? How many genes does it include? How much variation do you see in the data?
# To answer these questions, log-scale the data, calculate per-gene median expression ranges, then make a
# density plot showing those results.
print(f"Expression matrix dimensions: {exp_data.shape}")
print(exp_data.head())

num_expression_data = exp_data.apply(pd.to_numeric, errors='coerce')
num_expression_data = num_expression_data.dropna(axis=1, how='all')
print(f"Numeric expression matrix dimensions: {num_expression_data.shape}")

log_expression_data = np.log2(num_expression_data + 1)

metadata = pd.DataFrame({
    'sample': exp_data.columns,
    'group': ['barren'] * 35 + ['environmentally enriched'] * 35
})

# Generate a PCA plot
pca = PCA(n_components=2)
pca_result = pca.fit_transform(exp_data.T)
pca_dataframe = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
pca_dataframe['sample'] = metadata['sample']
pca_dataframe = pca_dataframe.merge(metadata, on='sample')

plt.figure(figsize=(10, 7))
sns.scatterplot(data=pca_dataframe, x='PC1', y='PC2', hue='group', palette='viridis', s=100)
plt.title('PCA of Gene Expression Data')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title='Group')
plt.savefig('results/PCA_plot.png')
plt.show()

# Generate a t-SNE plot
tsne = TSNE(n_components=2, random_state=42)
tsne_result = tsne.fit_transform(exp_data.T)
tsne_dataframe = pd.DataFrame(data=tsne_result, columns=['t-SNE1', 't-SNE2'])
tsne_dataframe['sample'] = metadata['sample']
tsne_dataframe = tsne_dataframe.merge(metadata, on='sample')

plt.figure(figsize=(10, 7))
sns.scatterplot(data=tsne_dataframe, x='t-SNE1', y='t-SNE2', hue='group', palette='viridis', s=100)
plt.title('t-SNE of Gene Expression Data')
plt.xlabel('t-SNE Component 1')
plt.ylabel('t-SNE Component 2')
plt.legend(title='Group')
plt.savefig('results/tSNE_plot.png')
plt.show()

# Generate a UMAP plot
umap_result = umap.UMAP(n_components=2, random_state=42).fit_transform(exp_data.T)
umap_dataframe = pd.DataFrame(data=umap_result, columns=['UMAP1', 'UMAP2'])
umap_dataframe['sample'] = metadata['sample']
umap_dataframe = umap_dataframe.merge(metadata, on='sample')

plt.figure(figsize=(10, 7))
sns.scatterplot(data=umap_dataframe, x='UMAP1', y='UMAP2', hue='group', palette='viridis', s=100)
plt.title('UMAP of Gene Expression Data')
plt.xlabel('UMAP Component 1')
plt.ylabel('UMAP Component 2')
plt.legend(title='Group')
plt.savefig('results/UMAP_plot.png')
plt.show()

# Perform differential analysis on the samples from your two groups (in diff analysis file)
import pandas as pd

DESeq2_results = pd.read_csv('/Users/aliciagrant/Desktop/BioInfo/results/alicia_diff_expression_results.csv', sep=',')

print("Columns in the deseq_results DataFrame:")
print(DESeq2_results.columns)

if 'padj' in DESeq2_results.columns:
    significant_genes = DESeq2_results[DESeq2_results['padj'] < 0.05]
    top_degs = significant_genes.sort_values(by='padj').head(50)

    DESeq2_results.to_csv('results/full_results.csv', index=False)
    # Create a table of the top 50 differentially expressed genes.
    top_degs.to_csv('results/top_50_degs.csv', index=False)
    print("Full results and top 50 DEGs saved.")
else:
    print("The 'padj' column is not found in the DataFrame. Please check the column names.")


# Create a volcano plot of your Make sure to label the axes and provide a legend.

plt.figure(figsize=(10, 6))
plt.scatter(DESeq2_results['log2FoldChange'], -np.log10(DESeq2_results['padj']),
            c=np.where(DESeq2_results['padj'] < 0.05, 'red', 'grey'), alpha=0.7)

plt.scatter([], [], c='red', label='Significant (padj < 0.05)')
plt.scatter([], [], c='grey', label='Non-significant (padj >= 0.05)')
plt.legend()
plt.title('Volcano Plot')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 Adjusted P-value')
plt.axhline(y=-np.log10(0.05), color='blue', linestyle='--')
plt.axvline(x=0, color='black', linestyle='--')
plt.grid()
plt.savefig('results/volcano_plot.png')
plt.show()

# Extract the list of significantly differentially expressed genes, and generate a heatmap
plt.figure(figsize=(12, 8))
heatmap_data = significant_genes.pivot(index='Gene', columns='log2FoldChange', values='baseMean')
plt.gca().set_facecolor('black')
sns.heatmap(heatmap_data, cmap='coolwarm', annot=True)
plt.title('Heatmap of Significant DEGs')
plt.xlabel('Log2 Fold Change')
plt.ylabel('Gene')
plt.savefig('results/heatmap_significant_degs.png')
plt.show()

diff_exp_genes = pd.read_csv('results/full_results.csv')

gene_list = diff_exp_genes['Gene'].tolist()

# Reactome enrichment analysis
g_profiler = GProfiler(return_dataframe=True)

try:
    results = g_profiler.profile(organism='sscrofa',
                                 query=gene_list,
                                 sources=['REAC'])
    if results.empty:
        print("No results found for the given gene list.")
    else:
        results.to_csv('results/reactome_enrichment_results.csv', index=False)
        print("Reactome enrichment results saved.")
except Exception as e:
    print(f"An error occurred: {e}")