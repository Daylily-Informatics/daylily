import sys
import pandas as pd

# Read command-line arguments
infile = sys.argv[1]
outfile = sys.argv[2]

# Read the table into a pandas DataFrame
df = pd.read_csv(infile, sep="\t")

# Split the 'rule' column into two columns
df[["rule_prefix", "rule_suffix"]] = df["rule"].str.split(".", n=1, expand=True)

# Create a new first column by joining the original first two columns with a '-'
df.insert(0, "combined_rule", df["rule"] + "-" + df["sample"])

# Save the modified table to a new file
df.to_csv(outfile, sep="\t", index=False)
