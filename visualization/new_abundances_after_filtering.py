import pandas as pd
import glob
import os


kreport_folder = '/home/rachel/kreport/'
meta = pd.read_csv('/home/rachel/kreport/meta.csv')


kreport_files = glob.glob(os.path.join(kreport_folder, '*_kreport'))
abundance_dict = {}

for file in kreport_files:
    sample_name = os.path.basename(file).replace('_kreport', '')

    df = pd.read_csv(file, sep='\t', header=None, names=[
        'percent', 'reads', 'assigned', 'rank_code', 'taxid', 'name'
    ])
    df = df[df['rank_code'].isin(['G'])]
    df = df[df['reads'] > 5000]
    
    total_reads = df['reads'].sum()
    
    if total_reads > 0:
        df['relative_abundance'] = df['reads'] / total_reads
    else:
        
        df['relative_abundance'] = 0
    
    df['taxon'] = df['rank_code'].map({'G': 'Genus'}) + ':' + df['name'].str.strip()

    abundance_dict[sample_name] = df.set_index('taxon')['relative_abundance']


abundance = pd.DataFrame.from_dict(abundance_dict, orient='columns').fillna(0)
abundance.index.name = 'taxon'

print(abundance.head())

abundance.to_csv('/home/rachel/kreport/abundancegen5k.csv')
