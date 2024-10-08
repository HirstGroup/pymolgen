import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

"""
# Example data: Y-values for each depth ranging from -20 to 0
np.random.seed(42)  # For reproducibility
data1 = np.random.uniform(-20, 0, 100)  # 100 random values for Depth 1
data2 = np.random.uniform(-20, 0, 100)  # 100 random values for Depth 2
data3 = np.random.uniform(-20, 0, 100)  # 100 random values for Depth 3
data4 = np.random.uniform(-20, 0, 100)  # 100 random values for Depth 4
data5 = np.random.uniform(-20, 0, 100)  # 100 random values for Depth 5

# Depth labels
depths = ['Depth 1', 'Depth 2', 'Depth 3', 'Depth 4', 'Depth 5']

# Create a DataFrame for plotting
df = pd.DataFrame({
    'Depth': ['Depth 1'] * len(data1) + ['Depth 2'] * len(data2) + 
             ['Depth 3'] * len(data3) + ['Depth 4'] * len(data4) + 
             ['Depth 5'] * len(data5),
    'Y Values': np.concatenate([data1, data2, data3, data4, data5])
})
"""

parser = argparse.ArgumentParser(description="Make violin plot")

# Required arguments
parser.add_argument('-i','--input', help='Input txt file from pymolgen',required=True)
parser.add_argument('-o','--output', help='Output File Name',required=True)

args = parser.parse_args()


def count_depth(row):

	return len(row['fragments'].split('-')) - 1


df = pd.read_csv(args.input, sep=':', header=None, names=['fragments', 'bonds', 'bp1', 'bp2', 'mw'])

df['depth'] = df.apply(count_depth, axis=1)

df['log-bp1'] = np.log(df['bp1'])

print(df)

# Create the violin plot
plt.figure(figsize=(10, 6))
sns.violinplot(x='depth', y='log-bp1', data=df)

# Customize the plot
plt.title('Build Probabilities at Different Depths')
plt.xlabel('Depth')
plt.ylabel('log(Build Probability)')

# Show the plot
plt.savefig(args.output)
