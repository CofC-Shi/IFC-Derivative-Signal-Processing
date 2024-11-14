import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

# Load datasets
data1 = pd.read_csv('Path/to/Dataset1')
data2 = pd.read_csv('Path/to/Dataset2')

# Calculate Peak-to-Peak Voltage
data1['PeakToPeakVoltage'] = data1['PosPeak'] - data1['NegPeak']
data2['PeakToPeakVoltage'] = data2['PosPeak'] - data2['NegPeak']

# Filter for Channel 2 in both datasets
data1_channel2 = data1[data1['Channel'] == 2]
data2_channel2 = data2[data2['Channel'] == 2]

# Create the plot
fig, ax = plt.subplots()

# Plot Dataset 1
ax.scatter(data1_channel2['PeakToPeakVoltage'], data1_channel2['Time_ms'], color='blue', alpha=0.7,
           label='4um (Dataset 1)', edgecolor='k')
# draw_ellipse(data1_channel2['PeakToPeakVoltage'], data1_channel2['Time_ms'], ax, 'blue')

# Plot Dataset 2
ax.scatter(data2_channel2['PeakToPeakVoltage'], data2_channel2['Time_ms'], color='red', alpha=0.7,
           label='7um (Dataset 2)')
# draw_ellipse(data2_channel2['PeakToPeakVoltage'], data2_channel2['Time_ms'], ax, 'red')

# Set log scale for both axes
ax.set_xscale('log')
ax.set_yscale('log')

# Set labels and title
ax.set_xlabel('Peak-to-Peak Voltage (V)', fontsize=14)
ax.set_ylabel('Half-width Transient Time (ms)', fontsize=14)
# plt.title('Scatter Plot for Two Datasets (4um and 7um) with Elliptical Distributions')
plt.legend(fontsize=12)
# plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.show()