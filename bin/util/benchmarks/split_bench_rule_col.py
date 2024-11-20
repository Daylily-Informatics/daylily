import sys
import pandas as pd

# Read command-line arguments
infile = sys.argv[1]
outfile = sys.argv[2]

# Read the table into a pandas DataFrame
df = pd.read_csv(infile, sep="\t")

# Split the 'rule' column into two columns
df[["rule_prefix", "rule_suffix"]] = df["rule"].str.split(".", n=1, expand=True)

# Reorder columns to switch the first and second columns
columns = df.columns.tolist()
columns[0], columns[1] = columns[1], columns[0]  # Swap the first and second columns
df = df[columns]

# Save the modified table to a new file
df.to_csv(outfile, sep="\t", index=False)
