import pandas as pd
import numpy as np
from scipy import stats
import os

expression_data = pd.read_csv('ERP105973.tsv', sep='\t', index_col=0)
print(f"Expression matrix dimensions: {expression_data.shape}")

numeric_expression_data = expression_data.apply(pd.to_numeric, errors='coerce')
numeric_expression_data = numeric_expression_data.dropna(axis=1, how='all')
print(f"Numeric expression matrix dimensions: {numeric_expression_data.shape}")

log_expression_data = np.log2(numeric_expression_data + 1)

metadata = pd.DataFrame({
    'sample': numeric_expression_data.columns,
    'group': ['barren'] * 35 + ['environmentally enriched'] * 35
})
metadata.set_index('sample', inplace=True)

differential_results = []

available_genes = log_expression_data.index.tolist()
print(f"Available genes in log expression data: {available_genes[:10]}")

for gene in available_genes:
    if gene not in log_expression_data.index:
        print(f"Gene {gene} not found in log expression data; skipping...")
        continue

    barren_group = log_expression_data.loc[gene][metadata['group'] == 'barren'].values
    enriched_group = log_expression_data.loc[gene][metadata['group'] == 'environmentally enriched'].values

    if len(barren_group) == 0 or len(enriched_group) == 0:
        print(f"Not enough data for gene {gene}; skipping...")
        continue

    t_stat, p_value = stats.ttest_ind(barren_group, enriched_group, equal_var=False)
    log2_fold_change = np.log2(np.mean(enriched_group) / np.mean(barren_group)) if np.mean(barren_group) != 0 else np.nan

    differential_results.append({
        'Gene': gene,
        'log2FoldChange': log2_fold_change,
        'p_value': p_value,
        'padj': p_value
    })

differential_results_df = pd.DataFrame(differential_results)

result_file = 'results/alicia_diff_expression_results.csv'
differential_results_df.to_csv(result_file, index=False)
print(f"Differential expression results saved to {result_file}.")
