import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns


abundance = pd.read_csv('/home/rachel/kreport/abundancesp.csv', index_col=0)

meta = pd.read_csv('/home/rachel/kreport/meta.csv')


results = []
vaches = meta['vache_id'].unique()

for vache in vaches:
    samples = meta[meta['vache_id'] == vache].sort_values('condition')['sample'].values
    
    if len(samples) != 2:
        print(f" {samples} not found")
        continue
    
    sample_before, sample_after = samples

    if sample_before not in abundance.columns or sample_after not in abundance.columns:
        print(f"Samples {sample_before} or {sample_after} not found")
        continue

    data = abundance[[sample_before, sample_after]].fillna(0)

    ranks_before = data[sample_before].rank()
    ranks_after = data[sample_after].rank()

    rho, pval = spearmanr(ranks_before, ranks_after)

    results.append({
        'vache_id': vache,
        'spearman_rho': rho,
        'p_value': pval
    })


results_df = pd.DataFrame(results)

print("\n=== Résultats Spearman ===")
print(results_df)

results_df.to_csv('spearman_resultssp.csv', index=False)
sns.set(style="whitegrid")

plt.figure(figsize=(10, 6))
barplot = sns.barplot(
    data=results_df,
    x='vache_id',
    y='spearman_rho',
    palette='viridis'
)

for index, row in results_df.iterrows():
    barplot.text(
        index,
        row['spearman_rho'] + 0.02,
        f"p={row['p_value']:.3f}",
        ha='center',
        va='bottom',
        fontsize=9
    )

plt.title('Spearman correlation by cow')
plt.ylabel('Spearman rho')
plt.xlabel('Cow ID')
plt.ylim(0, 1.05)  
plt.tight_layout()

plt.savefig('spearman_barplotgsp.png', dpi=300)
plt.show()


