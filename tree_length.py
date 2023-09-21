# plot tree length distribution as a function of number of boosts

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import dendropy
import seaborn as sns
import numpy as np

import utils

ds = 1
dim = 3
ln_crv = 0
burnin = 1
embed = "up"
ds_path = os.path.join("..", "..", "dodo-experiments", "analysis", f"ds{ds}")
max_mix = 1
trials = [f"boosts_gtr/d3_lr1_i3_b1_boosts" for i in range(1, max_mix+1)]

exp_path = os.path.join(ds_path, "vi", "up_nj")
dodo_paths = [os.path.join(exp_path, f"{trial}", "samples.t") for trial in trials]

# Plot
fig = plt.figure(figsize=(4, 3))
fig.gca().set_position([0.19, 0.15, 0.80, 0.84])

fig.gca().set_xlabel("Tree length")
fig.gca().set_xlim(0.28, 0.62)
# fig.gca().annotate("A", (0.05, 0.85), xycoords="figure fraction", fontsize="18")
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
    ax=fig.gca(),
    color="k",
    linestyle="-",
    label=f"MrBayes run {1}",
    linewidth=2,
)
print(len(mb_tree_lengths))
sns.kdeplot(mb_tree_lengths[1][burnin:], ax=fig.gca(), color='k', linestyle='-', label=f"MrBayes run {2}")

# Extract Dodonaphy tree lengths
for i, dodo_path in enumerate(dodo_paths):
    print(dodo_path)
    trees = dendropy.TreeList.get(path=dodo_path, schema="nexus")
    dodo_tree_lengths = [tree.length() for tree in trees]
    sns.kdeplot(
        dodo_tree_lengths[burnin:],
        ax=fig.gca(),
        color="b",
        linestyle=lineStyles[i],
        label="Dodonaphy",
    )

plt.legend()
print(np.mean(dodo_tree_lengths[burnin:]))
print(np.mean(mb_tree_lengths[0][burnin:]))
print(np.var(dodo_tree_lengths[burnin:]))
print(np.var(mb_tree_lengths[0][burnin:]))

# save the figure
path_save = os.path.join(".", "out", "tree_length")
utils.savefig_bioinformatics(path_save)
