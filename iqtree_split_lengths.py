# compare split lengths between dodonaphy ml and iqtree
import matplotlib.pyplot as plt
import os
import dendropy

ds = 1
dim = 3
ln_crv = 0
burnin = 1000
embed = "up"
ds_path = os.path.join( "analysis", f"ds{ds}")
trials = ['']  #("", "_1", "_2", "_3", "_4", "_5", "_6", "_7", "_8", "_9")

#Plot
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(4, 4)
ax.set_position([0.20, 0.11, 0.79, 0.88])
ax.set_xlabel('Dodonaphy Branch length')
ax.set_ylabel('True Branch length')
# ax.set_xlim(0.28, 0.52)
ax.annotate("C", (.05, .85), xycoords="figure fraction", fontsize="18")

# dictionary of split lengths
both_split_lens_int = {}
both_split_lens_leaf = {}

# iqtree
path_iqtree = os.path.join(ds_path, "iqtree", 'DS.treefile')
iqtree = dendropy.Tree.get(
    path=path_iqtree,
    schema="nexus")

for edge in iqtree.split_bitmask_edge_map.values():
    bit_string = edge.bipartition.split_as_bitstring()
    if edge.length is not None:
        if edge.is_leaf():
            if bit_string in both_split_lens_leaf:
                both_split_lens_leaf[bit_string][0].append(edge.length)
            else:
                both_split_lens_leaf[bit_string] = [[edge.length], []]
        else:
            if bit_string in both_split_lens_int:
                both_split_lens_int[bit_string][0].append(edge.length)
            else:
                both_split_lens_int[bit_string] = [[edge.length], []]

#Dodonaphy
dodo_only_lens_int = {}
dodo_only_lens_leaf = {}
dodo_path = os.path.join(ds_path, "hmap", "nj", "gammadir", "lr2_tau5_iqtree_GTR", "mape.t")
dodo_tree = dendropy.Tree.get(
    path=dodo_path,
    schema="nexus")
for edge in dodo_tree.split_bitmask_edge_map.values():
    bit_string = edge.bipartition.split_as_bitstring()
    if edge.length is not None:
        if edge.is_leaf():
            if bit_string in both_split_lens_leaf:
                both_split_lens_leaf[bit_string][1].append(edge.length)
            else:
                dodo_only_lens_leaf[bit_string] = edge.length
        else:
            if bit_string in both_split_lens_int:
                both_split_lens_int[bit_string][1].append(edge.length)
            else:
                dodo_only_lens_int[bit_string] = edge.length

# plot mean lengths
dodo_means = []
iq_means = []
full = False
stats = (both_split_lens_int, both_split_lens_leaf)
colors = ('r', 'b')
markers = ('d', 'o')
for split_dict, color, marker in zip(stats, colors, markers):
    for splits, lengths in split_dict.items():
        if len(lengths[1]) > 0:
            if full:
                for iq_len in lengths[0]:
                    for dodo_len in lengths[1]:
                        plt.plot(dodo_len, iq_len, linestyle='none', color=color, marker=marker)
                        print(f"{dodo_len} vs {iq_len}")
                        # dodo_means.append(dodo_len)
                        # mb_means.append(mb_len)
            else:
                dodo_len = sum(lengths[1]) / len(lengths[1])
                iq_len = sum(lengths[0]) / len(lengths[0])
                alpha=(len(lengths[0])/25000)**.5
                plt.plot(dodo_len, iq_len, linestyle='none', color=color, marker=marker, alpha=alpha)
                dodo_means.append(dodo_len)
                iq_means.append(iq_len)
            # print(f'{x}, {y}')


x = (0, max(dodo_means))
plt.plot(x, x, color='k')

# fp = os.path.join(meta.path_save, "split_lengths.eps")
# fp = os.path.join(ds_path, "mcmc", f"{embed}_nj", f"d{dim}_k{ln_crv}{trial}", "results", "split_lengths_colour.svg")
# fp = os.path.join(ds_path, "vi", f"{embed}_nj", "con", f"d5_lr2_i3_b10", "split_lengths.svg")
# fp = os.path.join(ds_path, "hmap", f"nj", f"lr2_tau4_raw_lrHalved", "results", "split_lengths.eps")

# plt.savefig(fp, format="svg",transparent=True)
plt.show()

# plt.hist2d(lengths[1], lengths[0], bins=50)
# plt.show()
