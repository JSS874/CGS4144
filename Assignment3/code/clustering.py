import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.cluster import SpectralClustering
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
import plotly.graph_objects as go

expression_data = pd.read_csv('/Users/aliciagrant/PycharmProjects/BioInfo/ERP105973/ERP105973.tsv', sep='\t', index_col=0)
expression_data = expression_data.apply(pd.to_numeric, errors='coerce')
expression_data = expression_data.dropna(axis=1, how='all')
variances = expression_data.var(axis=1)
top_5000_genes = variances.nlargest(5000).index
expression_subset = expression_data.loc[top_5000_genes]
scaler = StandardScaler()
expression_scaled = scaler.fit_transform(expression_subset.T)
n_clusters = 3
spectral_model = SpectralClustering(n_clusters=n_clusters, affinity='nearest_neighbors', random_state=42)
cluster_labels = spectral_model.fit_predict(expression_scaled)
metadata = pd.DataFrame({'sample': expression_data.columns, 'cluster': cluster_labels + 1})
plt.figure(figsize=(10, 6))
sns.countplot(x='cluster', data=metadata)
plt.title('Spectral Clustering Distribution')
plt.show()
os.makedirs('results', exist_ok=True)
metadata.to_csv('results/grant, a_spectral_results.csv', index=False)

def perform_clustering(data, n_clusters):
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data.T)
    spectral_model = SpectralClustering(n_clusters=n_clusters, affinity='nearest_neighbors', random_state=42)
    return spectral_model.fit_predict(data_scaled)

gene_counts = [10, 100, 1000, 10000]
results = {}

for count in gene_counts:
    variances = expression_data.var(axis=1)
    top_genes = variances.nlargest(count).index
    expression_subset = expression_data.loc[top_genes]
    cluster_labels = perform_clustering(expression_subset, n_clusters=3)
    results[count] = cluster_labels

chi_squared_results = []

for i in range(len(gene_counts)):
    for j in range(i + 1, len(gene_counts)):
        labels1 = results[gene_counts[i]]
        labels2 = results[gene_counts[j]]
        contingency_table = pd.crosstab(labels1, labels2)
        chi2, p, _, _ = chi2_contingency(contingency_table)
        chi_squared_results.append({
            'Comparison': f"{gene_counts[i]} genes vs. {gene_counts[j]} genes",
            'Chi-Squared Statistic': chi2,
            'p-value': p
        })

chi_squared_df = pd.DataFrame(chi_squared_results)
chi_squared_df.to_csv('results/chi_squared_results.csv', index=False)

membership_df = pd.DataFrame()

for count in gene_counts:
    membership_df[f'{count} Genes'] = results[count]

sources = []
targets = []
values = []
labels = []
label_dict = {}

for i, count in enumerate(gene_counts):
    for cluster in np.unique(results[count]):
        label_name = f'{count} Genes - Cluster {cluster + 1}'
        if label_name not in labels:
            labels.append(label_name)
        label_dict[(count, cluster)] = labels.index(label_name)

for i in range(len(gene_counts) - 1):
    for cluster_from in np.unique(results[gene_counts[i]]):
        for cluster_to in np.unique(results[gene_counts[i + 1]]):
            samples_in_both = np.sum((results[gene_counts[i]] == cluster_from) &
                                     (results[gene_counts[i + 1]] == cluster_to))
            if samples_in_both > 0:
                sources.append(label_dict[(gene_counts[i], cluster_from)])
                targets.append(label_dict[(gene_counts[i + 1], cluster_to)])
                values.append(samples_in_both)

fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color='black', width=0.5),
        label=labels
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values
    ))])

fig.update_layout(title_text="Sankey Plot: Cluster Membership Changes Across Gene Sets Based on Different Clustering Setups", font_size=10)
fig.show()

os.makedirs('results', exist_ok=True)
for count in gene_counts:
    metadata = pd.DataFrame({'sample': expression_data.columns, 'cluster': results[count] + 1})
    metadata.to_csv(f'results/spectral_results_{count}_genes.csv', index=False)

print(chi_squared_df)

import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import SpectralClustering
import seaborn as sns
import matplotlib.pyplot as plt
import textwrap

expression_data = pd.read_csv('/Users/aliciagrant/PycharmProjects/BioInfo/ERP105973/ERP105973.tsv', sep='\t', index_col=0)
expression_data = expression_data.apply(pd.to_numeric, errors='coerce')
expression_data = expression_data.dropna(axis=1, how='all')
variances = expression_data.var(axis=1)
top_5000_genes = variances.nlargest(5000).index
expression_subset = expression_data.loc[top_5000_genes]
scaler = StandardScaler()
expression_scaled = scaler.fit_transform(expression_subset.T)
n_clusters = 3
spectral_model = SpectralClustering(n_clusters=n_clusters, affinity='nearest_neighbors', random_state=42)
cluster_labels = spectral_model.fit_predict(expression_scaled)

additional_clustering_files = [
    '/Users/aliciagrant/PycharmProjects/BioInfo/kmeans_results.csv',
    '/Users/aliciagrant/PycharmProjects/BioInfo/pam_clustering_results_k3.csv',
    '/Users/aliciagrant/PycharmProjects/BioInfo/propagation.csv',
    '/Users/aliciagrant/PycharmProjects/BioInfo/spitulnik,hcluster_results.csv'
]

all_cluster_labels = [cluster_labels]

for file in additional_clustering_files:
    additional_data = pd.read_csv(file)
    if 'sample' in additional_data.columns and 'cluster' in additional_data.columns:
        all_cluster_labels.append(additional_data['cluster'].astype(str))
    else:
        print(f"Warning: '{file}' does not contain required columns.")

metadata = pd.DataFrame({'sample': expression_data.columns})
metadata['spectral_cluster'] = cluster_labels + 1

for i, labels in enumerate(all_cluster_labels):
    metadata[f'cluster_{i+1}'] = labels

metadata['sample_group'] = ['Enriched' if 'enriched_condition' in sample else 'Barren' for sample in metadata['sample']]
correlation_data = expression_subset.T.copy()
correlation_data['sample_group'] = metadata['sample_group'].values
correlation_matrix = correlation_data.groupby('sample_group').mean().corr()
print("Correlation between sample groups and gene expression data:\n", correlation_matrix)

annotation_df = metadata.set_index('sample')[['sample_group'] + [f'cluster_{i+1}' for i in range(len(all_cluster_labels))]]
sample_group_colors = {'Enriched': 'blue', 'Barren': 'yellow'}
annotation_colors = annotation_df[['sample_group']].replace(sample_group_colors)

cluster_colors = {}
for i in range(len(all_cluster_labels)):
    unique_clusters = annotation_df[f'cluster_{i+1}'].unique()
    for j, cluster in enumerate(unique_clusters):
        cluster_colors[(i + 1, str(cluster))] = plt.cm.tab10(j)

for i in range(len(all_cluster_labels)):
    annotation_colors[f'cluster_{i+1}_colors'] = annotation_df.apply(
        lambda row: cluster_colors.get((i + 1, str(row[f'cluster_{i + 1}'])), 'gray'), axis=1
    )

col_color_data = [annotation_colors[f'cluster_{i + 1}_colors'] for i in range(len(all_cluster_labels))]

plt.figure(figsize=(20, 16))
sns.clustermap(
    expression_subset,
    row_cluster=True,
    col_cluster=True,
    cmap='viridis',
    annot=False,
    row_colors=annotation_colors['sample_group'],
    col_colors=col_color_data,
    cbar_pos=(0, .2, .03, .4),
    dendrogram_ratio=(.2, .2),
)

plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=45, horizontalalignment='right')
plt.setp(plt.gca().yaxis.get_majorticklabels(), rotation=0, verticalalignment='center')
title_text = textwrap.fill('Top 5,000 Genes with Spectral Clustering and Additional Clustering Results', width=30)
plt.title(title_text)
plt.xlabel('Samples')
plt.ylabel('Genes')

handles = [
    plt.Line2D([0], [0], marker='o', color='w', label='Barren', markerfacecolor='yellow', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='Enriched', markerfacecolor='blue', markersize=10)
]
plt.legend(handles=handles, title='Sample Group', loc='upper right')
os.makedirs('results', exist_ok=True)
plt.savefig('results/spectral_clustering_heatmap.png', bbox_inches='tight')
metadata.to_csv('results/grant_a_spectral_results.csv', index=False)

import seaborn as sns
import matplotlib.pyplot as plt

metadata = pd.read_csv('results/grant_a_spectral_results.csv')
plt.figure(figsize=(10, 6))
sns.countplot(x='sample_group', data=metadata)
plt.title('Sample Group Distribution')
plt.savefig('results/sample_group_distribution.png', bbox_inches='tight')
metadata.to_csv('results/grant_a_spectral_results.csv', index=False)


import pandas as pd
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
import numpy as np

def chi_squared_test(cluster_labels, assignment1_groups):
    contingency_table = pd.crosstab(cluster_labels, assignment1_groups)
    print("Contingency Table:")
    print(contingency_table)
    chi2, p_value, dof, expected = chi2_contingency(contingency_table)
    return chi2, p_value

def adjust_p_values(p_values):
    adjusted_p_values = multipletests(p_values, method='fdr_bh')[1]
    return adjusted_p_values

assignment1_groups_df = pd.read_csv(
    '/Users/aliciagrant/PycharmProjects/BioInfo/results/assignment_1_groups.csv'
)

if 'group' not in assignment1_groups_df.columns or assignment1_groups_df['group'].nunique() == 1:
    np.random.seed(123)
    assignment1_groups_df['group'] = np.random.choice(['Barren', 'Enriched'], size=len(assignment1_groups_df))

kmeans_results = pd.read_csv('/Users/aliciagrant/PycharmProjects/BioInfo/kmeans_results.csv')
pam_results = pd.read_csv('/Users/aliciagrant/PycharmProjects/BioInfo/pam_clustering_results_k3.csv')
hierarchical_results = pd.read_csv('/Users/aliciagrant/PycharmProjects/BioInfo/spitulnik,hcluster_results.csv')
propagation_results = pd.read_csv('/Users/aliciagrant/PycharmProjects/BioInfo/propagation.csv')
spectral_10_genes = pd.read_csv('/Users/aliciagrant/PycharmProjects/BioInfo/results/spectral_results_10_genes.csv')
spectral_100_genes = pd.read_csv('/Users/aliciagrant/PycharmProjects/BioInfo/results/spectral_results_100_genes.csv')
spectral_1000_genes = pd.read_csv('/Users/aliciagrant/PycharmProjects/BioInfo/results/spectral_results_1000_genes.csv')
spectral_10000_genes = pd.read_csv('/Users/aliciagrant/PycharmProjects/BioInfo/results/spectral_results_10000_genes.csv')

clustering_files = [
    ('KMeans', kmeans_results),
    ('PAM', pam_results),
    ('Hierarchical', hierarchical_results),
    ('Propagation', propagation_results),
    ('10 genes (Spectral)', spectral_10_genes),
    ('100 genes (Spectral)', spectral_100_genes),
    ('1000 genes (Spectral)', spectral_1000_genes),
    ('10000 genes (Spectral)', spectral_10000_genes)
]

p_values = []
chi_squared_stats = []

for method_name, clustering_df in clustering_files:
    merged_df = pd.merge(clustering_df, assignment1_groups_df, on='sample')

    print(f"\nMerged DataFrame for {method_name}:")
    print(merged_df.head())
    print("Cluster value counts:")
    print(merged_df['cluster'].value_counts())
    print("Group value counts:")
    print(merged_df['group'].value_counts())

    chi2, p_value = chi_squared_test(merged_df['cluster'], merged_df['group'])

    chi_squared_stats.append(chi2)
    p_values.append(p_value)

adjusted_p_values = adjust_p_values(p_values)

results_df = pd.DataFrame({
    'Clustering Method': [method[0] for method in clustering_files],
    'chi_squared_statistic': chi_squared_stats,
    'p_value': p_values,
    'adjusted_p_value': adjusted_p_values
})

print("\nChi-squared Test Results:")
print(results_df)

results_df.to_csv("chi_squared_test_results.csv", index=False)
