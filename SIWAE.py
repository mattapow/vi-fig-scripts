# plot the final log SIWAE value as the number of mixtures changes
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# extract the SIWAE values
path_root = "../analysis/ds1/vi/up_nj/boosts_jc69_n100/"
siwae = []
std = []
n_mixtures = []
from_log = True
for m in range(1, 9):
    if from_log:
    # option 1: use log file. ELBO estimated from many final samples
        path_vi_log = path_root + f"d3_lr1_i3_b{m}_jc69_exp_boosts/vi.log"
        file_save = "boosts_SIWAE_jc69_exp_boosts.pdf"
        with open(path_vi_log, 'r') as f:
            lines = f.read().splitlines()
            line_siwae = lines[-3]
            this_siwae = -float(line_siwae.split()[-1])
            this_std = 0.0
    else:    
    # option 2: get mean elbo estimates from last n epcohs
        path_elbo = path_root + f"d3_lr1_i3_b{m}_boosts/elbo.txt"
        file_save = "boosts_SIWAE_last_iteration.pdf"
        with open(path_elbo, 'r') as f:
            lines = f.read().splitlines()
            n = 50  # number of lines to consider
            sum_lines = 0.0
            these_siwaes = []
            for i in range(-1, -(n+1), -1):
                these_siwaes.append(float(lines[i]))
            this_siwae = np.mean(these_siwaes)
            this_std = np.std(these_siwaes)

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
plt.savefig(path_save, format='pdf', dpi=600)
plt.show()
