import matplotlib.pyplot as plt
import xarray as xr
from cycler import cycler
from utils import add_annotation

plt.rcParams.update(
    {
        "figure.figsize": (8, 10),
        "figure.dpi": 200,
        "font.family": "serif",
        "font.size": 11,
        "text.usetex": True,
        "text.latex.preamble": "\n".join(
            [
                r"\usepackage{amsmath}",
                r"\usepackage[utf8]{inputenc}",
                r"\usepackage[T1]{fontenc}",
                r"\usepackage{siunitx}",
            ]
        ),
    }
)
RESOLUTIONS = [32, 64, 128, 256]
default_cycler = cycler(color=["#ff6666", "#990000", "#6666ff", "#000099", "#000000"])
plt.rc("axes", prop_cycle=default_cycler)
ds = xr.open_dataset("vorticity_pdfs.nc")

fig, axis = plt.subplots(4, 1)
titles = [
    "EPIC",
    "MPIC",
    "MONC with subsampling",
    "all models at $" + str(RESOLUTIONS[-1]) + "^3$",
]

for resolution in RESOLUTIONS:
    axis[0].stairs(ds["hist_epic_" + str(resolution)], ds["bin_edges"])
    axis[1].stairs(ds["hist_mpic_" + str(resolution)], ds["bin_edges"])
    axis[2].stairs(ds["hist_monc_" + str(resolution)], ds["bin_edges"])
axis[0].legend(["$" + str(i) + "^3$" for i in RESOLUTIONS])
axis[3].stairs(ds["hist_epic_" + str(RESOLUTIONS[-1])], ds["bin_edges"])
axis[3].stairs(ds["hist_mpic_" + str(RESOLUTIONS[-1])], ds["bin_edges"])
axis[3].stairs(ds["hist_monc_" + str(RESOLUTIONS[-1])], ds["bin_edges"])
for axnr, ax in enumerate(axis):
    ax.set_ylim(1e-8, 100000.0)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlim(1e-6, 400)
    ax.set_ylabel("probability density")
    ax.grid(linestyle="dashed")
axis[3].set_xlabel("$|\omega|$ (-)")
axis[3].legend(["EPIC", "MPIC", "MONC"])
add_annotation(axis[0], "EPIC", [0.5, 0.9], ha="center")
add_annotation(axis[1], "MPIC", [0.5, 0.9], ha="center")
add_annotation(axis[2], "MONC", [0.5, 0.9], ha="center")
add_annotation(axis[3], "All models at $" + str(RESOLUTIONS[-1]) + "^3$ grid points", [0.5, 0.9], ha="center")
fig.tight_layout()
plt.savefig("compare_o_pdf.png", bbox_inches="tight")
plt.savefig("compare_o_pdf.jpeg", bbox_inches="tight")
plt.savefig("compare_o_pdf.pdf", bbox_inches="tight")
