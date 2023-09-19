import os

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np

import utils

# Data from the table
datasets = ['DS1', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6', 'DS7', 'DS8']
column_labels = ['RAxML', 'IQ-TREE', 'Dodonaphy+', 'Dodonaphy', 'BioNJ']

data = np.array([
    [-6771.2, -6771.2, -6803.9, -6807.4, -6955.2],
    [-25722.2, -25722.2, -25794.3, -25804.3, -26226.8],
    [-30657.5, -30657.6, -30670.1, -31053.9, -33494.1],
    [-12751.1, -12751.1, -12764.3, -12794.7, -13088.2],
    [-7491.8, -7491.8, -7548.6, -7549.3, -7934.8],
    [-6137.6, -6137.6, -6143.6, -6147.4, -6293.2],
    [-33676.1, -33676.2, -33793.2, -34312.0, -36844.0],
    [-7759.1, -7759.2, -7811.1, -7838.1, -8179.8]
])

# Create subplots
fig, ax = plt.subplots()

# Plot trendlines for each column
# for i in range(len(column_labels)):
#     ax.plot(datasets, data[:, i] - data[:, 0], label=column_labels[i])
jitter = 0.01
for i, data_label in enumerate(datasets):
    x_jitter = (i-3.5)*jitter
    ax.scatter(
        np.arange(len(column_labels)) + x_jitter,
        data[i, :] - data[i, 0],
        label=data_label
    )


# Add grid and legend
ax.legend()
ax.grid()
ax.set_yscale("symlog")
ax.set_yscale('symlog', linthreshy=1)

# Set labels and title
ax.set_ylabel('Log Likelihood Difference', fontsize=14)
ax.set_xticks(np.arange(len(column_labels)))
ax.set_xticklabels(column_labels, rotation=0, ha='center')


# ax.ticklabel_format(style='sci', axis='y', scilimits=(-3, 0))
ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax.tick_params(axis='both', which='major', labelsize=12)
plt.tight_layout()

path_save = os.path.join(".", "out", "ml_performance")
utils.savefig_bioinformatics(path_save)
plt.show()
