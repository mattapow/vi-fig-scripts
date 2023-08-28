# plot tree length distribution as a function of number of boosts

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import dendropy
import seaborn as sns
import numpy as np


ds = 1
dim = 3
ln_crv = 0
burnin = 1
embed = "up"
ds_path = os.path.join("..", "analysis", f"ds{ds}")
max_mix = 1
trials = [f"d20_lr1_i1_b2_d20_b2_k5_iqtree_GTR" for i in range(1, max_mix+1)]

exp_path = os.path.join(ds_path, "vi", "up_nj")
dodo_paths = [os.path.join(exp_path, f"{trial}", "samples.t") for trial in trials]

# Plot
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(6, 4)
ax.set_position([0.1, 0.11, 0.83, 0.88])
ax.set_xlabel("Tree length")
ax.set_xlim(0.28, 0.62)
ax.annotate("A", (0.05, 0.85), xycoords="figure fraction", fontsize="18")
color_map = cm.get_cmap('viridis', max_mix)
lineStyles = [":", "-", "--", "-.", ":", "-", "--", "-.", ":"]

# Extract MrBayes tree lengths
mb_tree_lengths = []
for i in range(1, 3):
    mb_path = os.path.join(ds_path, "mb", f"gold9.run{i}.t")
    mb_trees = dendropy.TreeList.get(path=mb_path, schema="nexus")
    mb_tree_lengths.append([tree.length() for tree in mb_trees])
sns.kdeplot(
    mb_tree_lengths[0][burnin:],
    ax=ax,
    color="k",
    linestyle="-",
    label=f"MrBayes run {1}",
    linewidth=2,
)
print(len(mb_tree_lengths))
sns.kdeplot(mb_tree_lengths[1][burnin:], ax=ax, color='k', linestyle='-', label=f"MrBayes run {2}")

# Extract Dodonaphy tree lengths
for i, dodo_path in enumerate(dodo_paths):
    print(dodo_path)
    trees = dendropy.TreeList.get(path=dodo_path, schema="nexus")
    dodo_tree_lengths = [tree.length() for tree in trees]
    sns.kdeplot(
        dodo_tree_lengths[burnin:],
        ax=ax,
        color="b",
        linestyle=lineStyles[i],
        label=f"M={i+1}",
    )

plt.legend()
print(np.mean(dodo_tree_lengths[burnin:]))
print(np.mean(mb_tree_lengths[0][burnin:]))
print(np.var(dodo_tree_lengths[burnin:]))
print(np.var(mb_tree_lengths[0][burnin:]))

# save the figure
plt.savefig(f"tree_length.eps", format="eps", dpi=600)
plt.show()
