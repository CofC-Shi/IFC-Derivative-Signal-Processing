import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load datasets
data1 = pd.read_csv('Input\Path1.cvs')
data2 = pd.read_csv('Input\Path2.csv')

# Filter for Channel 2 in both datasets
data1_channel2 = data1[data1['Channel'] == 2]
data2_channel2 = data2[data2['Channel'] == 2]

# Determine the minimum class size for balanced sampling
min_size = min(len(data1_channel2), len(data2_channel2))

# Create the plot
fig, ax = plt.subplots()

# Sample each dataset based on the minimum class size
data1_sampled = data1_channel2.sample(n=min_size, random_state=1)
data2_sampled = data2_channel2.sample(n=min_size, random_state=1)

# Create the plot for panel (a) style
fig, ax = plt.subplots(figsize=(6, 5))  # ax is not a list when creating a single subplot

# Panel (a) style scatter plot
ax.scatter(data1_sampled['Amplitude'], data1_sampled['Time_ms'],
           color='blue', alpha=0.7, label='4μm (Dataset 1)', edgecolor='none')
ax.scatter(data2_sampled['Amplitude'], data2_sampled['Time_ms'],
           color='red', alpha=0.7, label='7μm (Dataset 2)', edgecolor='none')

# Set log scale for both axes in panel (a) style
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Peak Amplitude (V)', fontsize=16)
ax.set_ylabel('Transient Time (ms)', fontsize=16)
ax.legend(fontsize=16)

# Improve layout
plt.tight_layout()

# Display the plot
plt.show()
