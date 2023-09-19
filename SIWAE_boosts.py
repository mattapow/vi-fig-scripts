# plot the final log SIWAE value as the number of mixtures changes
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit

import utils

# extract the SIWAE values
path_root = "../../dodo-experiments/analysis/ds1/vi/up_nj/"
siwae = []
std = []
n_mixtures = []
for m in range(1, 11):
    # ELBO estimated from many final samples
    path_vi_log = path_root + f"d3_lr1_i1_b{m}_jc69_boosts/vi.log"
    file_save = "boosts_SIWAE"
    with open(path_vi_log, 'r') as f:
        lines = f.read().splitlines()
        line_siwae = lines[-3]
        this_siwae = -float(line_siwae.split()[-1])
        this_std = 0.0

    siwae.append(this_siwae)
    std.append(this_std)
    n_mixtures.append(m)

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(6, 4)
ax.set_position([0.16, 0.11, 0.83, 0.88])
# ax.annotate("B", (0.05, 0.85), xycoords="figure fraction", fontsize="18")
plt.errorbar(n_mixtures, siwae, yerr=std, fmt='ko-')
plt.xlabel("Number of mixtures")
plt.ylabel("log SIWAE")
plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

path_save = os.path.join("out", file_save)
utils.savefig_bioinformatics(path_save)
plt.show()
