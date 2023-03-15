import matplotlib.pyplot as plt
import xarray as xr
from cycler import cycler
from utils import add_annotation, setup_rcParams, remove_xticks, remove_yticks

setup_rcParams()

RESOLUTIONS = [32, 64, 128, 256]
default_cycler = cycler(color=["#ff6666", "#990000", "#6666ff", "#000099", "#000000"])
plt.rc("axes", prop_cycle=default_cycler)
ds = xr.open_dataset("humidity_pdfs.nc")

fig, axes = plt.subplots(2, 2, figsize=(8, 6), sharex=True, sharey=True)

for resolution in RESOLUTIONS:
    axes[0, 0].stairs(ds["hist_epic_" + str(resolution)], ds["bin_edges"])
    axes[0, 1].stairs(ds["hist_mpic_" + str(resolution)], ds["bin_edges"])
    axes[1, 0].stairs(ds["hist_monc_" + str(resolution)], ds["bin_edges"])
axes[0, 0].legend(["$" + str(i) + "^3$" for i in RESOLUTIONS], loc="upper left")
axes[1, 1].stairs(ds["hist_epic_" + str(RESOLUTIONS[-1])], ds["bin_edges"])
axes[1, 1].stairs(ds["hist_mpic_" + str(RESOLUTIONS[-1])], ds["bin_edges"])
axes[1, 1].stairs(ds["hist_monc_" + str(RESOLUTIONS[-1])], ds["bin_edges"])
for ax in axes.flat:
    ax.set_ylim(0, 2.0)
    ax.set_xlim(0, 0.08)
    ax.grid(linestyle="dashed")
axes[0, 0].set_ylabel("probability density")
axes[1, 0].set_ylabel("probability density")
axes[1, 1].set_xlabel("$q_l$ (-)")
axes[1, 0].set_xlabel("$q_l$ (-)")
remove_xticks(axes[0, 1])
remove_xticks(axes[0, 0])
remove_yticks(axes[0, 1])
remove_yticks(axes[1, 1])
axes[1, 1].legend(["EPIC", "MPIC", "MONC"], loc="upper left")
add_annotation(axes[0, 0], "EPIC", [0.96, 0.9], ha="right")
add_annotation(axes[0, 1], "MPIC", [0.96, 0.9], ha="right")
add_annotation(axes[1, 0], "MONC", [0.96, 0.9], ha="right")
add_annotation(
    axes[1, 1],
    "All models at \n $" + str(RESOLUTIONS[-1]) + "^3$ grid points",
    [0.96, 0.9],
    ha="right",
)
fig.tight_layout()
plt.savefig("compare_hl_pdf.png", bbox_inches="tight")
plt.savefig("compare_hl_pdf.jpeg", bbox_inches="tight")
plt.savefig("compare_hl_pdf.pdf", bbox_inches="tight")
