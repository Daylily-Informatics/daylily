import sys
import pandas as pd

infile=sys.argv[1]
outfile=sys.argv[2]

# read the table into a pandas DataFrame
df = pd.read_csv(infile, sep="\t")

# split the rule column into two columns
df[["rule_prefix", "rule_suffix"]] = df["rule"].str.split(".", n=1, expand=True)

# save the modified table to a new file
df.to_csv(outfile, sep="\t", index=False)
