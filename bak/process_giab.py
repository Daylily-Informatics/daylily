

import sys
import pandas as pd
import re

# Load the TSV file
df = pd.read_csv(sys.argv[1], sep='\t')

# Functions to parse coverage values (handle fractional and hybrid)
def parse_coverage(sample_str, tech_prefix):
    match = re.search(f'{tech_prefix}(\\d+-?\\d*)x', sample_str)
    if match:
        cov_str = match.group(1).replace('-', '.')
        return float(cov_str)
    return 0.0

# Apply parsing functions correctly
df['ILMN_cov'] = df['Sample'].apply(lambda x: parse_coverage(x, 'Ids') or parse_coverage(x, 'I'))
df['ONT_cov'] = df['Sample'].apply(lambda x: parse_coverage(x, 'Ods'))
df['UG_cov'] = df['Sample'].apply(lambda x: parse_coverage(x, 'Uds'))

# Check for any rows where coverage parsing failed completely
failed_parsing = df[(df['ILMN_cov'] + df['ONT_cov'] + df['UG_cov']) <= 0.0]

if not failed_parsing.empty:
    for idx, row in failed_parsing.iterrows():
        raise Exception(f"Coverage parsing failed for Sample: {row['Sample']}")

# Save the modified DataFrame to a new TSV
df.to_csv(sys.argv[2], sep='\t', index=False)
