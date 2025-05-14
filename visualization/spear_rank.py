import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns

abundance = pd.read_csv('/home/rachel/kreport/abundancegen5k.csv', index_col=0)

taxa_interet = [
    "Methanosarcina", "Methanosphaera", "Methanobrevibacter", "Ruminococcus",
    "Fibrobacter",
    "Butyrivibrio",
    "Prevotella", "Clostridium"
]

taxa_mask = abundance.index.to_series().apply(
    lambda x: any(taxon in x for taxon in taxa_interet)
)
abundance_interest = abundance[taxa_mask]

meta = pd.read_csv('/home/rachel/kreport/meta.csv'
condition_map = dict(zip(meta['sample'], meta['condition']))
abundance_interest = abundance_interest.rename(columns=condition_map)

before_samples = [col for col in abundance_interest.columns if 'before' in col]
after_samples = [col for col in abundance_interest.columns if 'after' in col]

before_mean = abundance_interest[before_samples].mean(axis=1)
after_mean = abundance_interest[after_samples].mean(axis=1)

before_rank = before_mean.rank(ascending=False)
after_rank = after_mean.rank(ascending=False)

rho, pval = spearmanr(before_rank, after_rank)
print(f"Spearman correlation for selected taxa: rho = {rho:.3f}, p-value = {pval:.3e}")

rank_change = (before_rank - after_rank).sort_values()

plt.figure(figsize=(10, 6))
sns.barplot(x=rank_change.index, y=rank_change.values, palette='coolwarm')
plt.xticks(rotation=90, fontsize=8)
plt.ylabel('Rank Change (Before - After)')
plt.title(f'Change in Ranking of Selected Taxa\nSpearman rho={rho:.2f}, p={pval:.2e}')
plt.tight_layout()
plt.show()
