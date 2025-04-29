import sys
import pandas as pd

# Load your data into a DataFrame (assuming 'data.tsv' is your TSV file)
df = pd.read_csv(sys.argv[1], sep='\t')

# Filter for specific SNP classes
classes = ['SNPts', 'SNPtv', 'INS_50', 'DEL_50', 'Indel_50']
df_filtered = df[df['SNPClass'].isin(classes)]

# Calculate mean, min, and max F-score by CmpFootprint and SNPClass
summary = df_filtered.groupby(['CmpFootprint', 'SNPClass']).agg(
    Mean_Fscore=('Fscore', 'mean'),
    Min_Fscore=('Fscore', 'min'),
    Max_Fscore=('Fscore', 'max')
).reset_index()

# Save the summary to a markdown file
summary.to_markdown('summary_results.md', index=False)

# Output to verify
print(summary.to_markdown(index=False))
