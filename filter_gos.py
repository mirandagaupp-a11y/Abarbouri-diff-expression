import pandas as pd

# Manually define the input and output file names
input_file = "final_gene_gos.tsv"
output_file = "filtered_gene_gos.tsv"

# Get true column names from the first line (remove leading #)
with open(input_file) as f:
    for line in f:
        if line.startswith("#"):
            column_line = line.lstrip("#").strip()
            break

# Split column names and load the data with correct headers
column_names = column_line.split('\t')
df = pd.read_csv(input_file, sep='\t', comment="#", header=None, names=column_names)

# Confirm column names
data_cols = ['gene_ontology_BLASTX', 'gene_ontology_BLASTP', 'EggNM.GOs']

# Filter: keep rows where at least one data column is not "."
filtered_df = df[~(df[data_cols] == ".").all(axis=1)]

# Write output
with open(output_file, 'w') as f:
    f.write("#" + "\t".join(column_names) + "\n")
    filtered_df.to_csv(f, sep='\t', index=False)

