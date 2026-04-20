import pandas as pd

# Function to extract gene ID from transcript ID
def extract_gene_id(transcript_id):
    return transcript_id.split('t')[0]

# Function to remove 'ko:' prefix from KEGG identifiers
def clean_kegg_id(kegg_id):
    return kegg_id.replace('ko:', '')

# Read the file into a DataFrame, assuming the first line is the header
df = pd.read_csv('eggnog_kegg_ko.tsv', sep=' ', header=0, names=['query', 'KEGG_ko'])

# Apply the function to extract gene IDs
df['gene_id'] = df['query'].apply(extract_gene_id)

# Clean the KEGG identifiers
df['KEGG_ko'] = df['KEGG_ko'].apply(clean_kegg_id)

# Group by gene_id and aggregate KEGG identifiers
grouped_df = df.groupby('gene_id')['KEGG_ko'].apply(lambda x: ','.join(sorted(set(x)))).reset_index()

# Save the result to a new TSV file
grouped_df.to_csv('eggnog_kegg_genes.tsv', sep='\t', index=False, header=['gene_id', 'KEGG_ko'])