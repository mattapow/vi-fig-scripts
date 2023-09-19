
# tanglegram of iqtree vs dodonaphy tree
import tanglegram as tg
import matplotlib.pyplot as plt
import pandas as pd
from dendropy import Tree as Tree
from dendropy.calculate import treecompare
import numpy as np
import os

import utils

# load the trees
# run_id = "8"
for ds in range(1, 2):
    path1 = f"/Users/151569/Projects/Dodonaphy/dodo-experiments/analysis/ds{ds}/iqtree_GTR/DS.treefile"
    dir2 = f"/Users/151569/Projects/Dodonaphy/dodo-experiments/analysis/ds{ds}/hmap/nj/None/bionj/lr2_tau5_n2000_k2_d3/"

    path2 = os.path.join(dir2, "mape.t")
    # path_save = os.path.join(dir2, "tangelgram"))

    label1 = "IQ-TREE"
    label2 = "Dodonaphy"

    tree1 = Tree.get_from_path(path1, "newick")
    tree2 = Tree.get_from_path(path2, "nexus", taxon_namespace=tree1.taxon_namespace)

    # generate the labels
    iq_labels = tree1.taxon_namespace.labels()
    dodo_labels = tree2.taxon_namespace.labels()

    # generate the distance matrix
    iq_dists = tree1.phylogenetic_distance_matrix()
    dodo_dists = tree2.phylogenetic_distance_matrix()

    # convert distances to numpy
    tree1_dists = np.zeros((len(iq_labels), len(iq_labels)))
    for i, t1 in enumerate(tree1.taxon_namespace[:]):
        for j, t2 in enumerate(tree1.taxon_namespace[:]):
            tree1_dists[i, j] = iq_dists(t1, t2)
    tree2_dists = np.zeros((len(dodo_labels), len(dodo_labels)))
    for i, t1 in enumerate(tree2.taxon_namespace[:]):
        for j, t2 in enumerate(tree2.taxon_namespace[:]):
            tree2_dists[i, j] = dodo_dists(t1, t2)

    # convert into pandas dataframe
    tree1_df = pd.DataFrame(tree1_dists, columns=iq_labels, index=iq_labels)
    tree2_df = pd.DataFrame(tree2_dists, columns=dodo_labels, index=dodo_labels)

    # Plot tanglegram
    fig = tg.plot(tree1_df, tree2_df, sort="step2side")

    # fix up the plot
    fig.set_size_inches(10, 8)
    ax_list = fig.axes
    ax_list[0].set_xlabel(label1)
    ax_list[1].set_xlabel(label2)
    ax_list[0].set_facecolor("w")
    ax_list[1].set_facecolor("w")
    ax_list[0].grid(visible=True, which="major", axis='x', color="whitesmoke")
    ax_list[1].grid(visible=True, which="major", axis='x', color="whitesmoke")
    # ax_list[0].annotate(
    #     f"RF dist = {treecompare.symmetric_difference(tree1, tree2)}\nwRF dist = {treecompare.weighted_robinson_foulds_distance(tree1, tree2)}",
    #     xy=(0.4, 0.94),
    #     xycoords='figure fraction')

    path_save = os.path.join(".", "out", "tanglegram")
    utils.savefig_bioinformatics(path_save)

    print(f"RF dist = {treecompare.symmetric_difference(tree1, tree2)}")
    print(f"wRF dist = {treecompare.weighted_robinson_foulds_distance(tree1, tree2)}")
    plt.show()

