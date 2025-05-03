import sys
import pandas as pd

# Load data from TSV file (replace 'input.tsv' with your file path)
data = pd.read_csv(sys.argv[1], sep='\t')

# Convert seconds to minutes
data['minutes'] = data['s'] / 60

# Group by 'rule' and calculate required statistics
summary = data.groupby('rule').agg(
    avg_task_cost=('task_cost', 'mean'),
    min_task_cost=('task_cost', 'min'),
    max_task_cost=('task_cost', 'max'),
    avg_minutes=('minutes', 'mean'),
    min_minutes=('minutes', 'min'),
    max_minutes=('minutes', 'max'),
    avg_cpu_efficiency=('cpu_efficiency', 'mean'),
    min_cpu_efficiency=('cpu_efficiency', 'min'),
    max_cpu_efficiency=('cpu_efficiency', 'max')
).reset_index()

# Display the summary
print(summary.to_markdown(index=False, floatfmt=".4f"))
