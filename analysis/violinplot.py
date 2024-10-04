import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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

# Create the violin plot
plt.figure(figsize=(10, 6))
sns.violinplot(x='Depth', y='Y Values', data=df)

# Customize the plot
plt.title('Violin Plot of Y Values at Different Depths')
plt.xlabel('Depth')
plt.ylabel('Y Values')

# Show the plot
plt.show()
