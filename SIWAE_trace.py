# plot the trace of the SIWAE value as the number of mixtures changes
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


fig, ax = plt.subplots(1, 1)
fig.set_size_inches(6, 4)
ax.set_position([0.16, 0.11, 0.83, 0.88])

max_mix = 10
color_map = cm.get_cmap('viridis', max_mix)
lineStyles = ["-", "--", "-.", ":"] * 3 
# ax.set_xlim(0., 1e3)
# ax.annotate("B", (0.05, 0.85), xycoords="figure fraction", fontsize="18")

path_root = "../analysis/ds1/vi/up_nj/"
for m in range(1, max_mix+1):
    path_vi = path_root + f"d3_lr1_i3_b{m}_boosts/elbo.txt"
    siwae = np.genfromtxt(path_vi, dtype=float, skip_header=0, unpack=True)
    # siwae = siwae_raw[:1001]
    plt.plot(siwae, color=color_map(m-1), label=f"M={m}", linestyle=lineStyles[m-1])

plt.xlabel("Iteration")
plt.ylabel("log SIWAE")
plt.legend(title="Number of mixtures", ncol=2)
# plt.savefig("boosts_trace.eps", format='eps', dpi=600)
plt.show()
