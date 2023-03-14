import matplotlib.pyplot as plt
import xarray as xr
from cycler import cycler

plt.rcParams.update(
    {
        "figure.figsize": (8, 10),
        "figure.dpi": 200,
        "font.family": "serif",
        "font.size": 11,
        "text.usetex": False,
    }
)

#    'text.latex.preamble': "\n".join([
#        r"\usepackage{amsmath}",
#        r"\usepackage[utf8]{inputenc}",
#        r"\usepackage[T1]{fontenc}",
#        r"\usepackage{siunitx}",
#        ])

RESOLUTIONS = [32, 64, 128, 256]
default_cycler = cycler(color=["#ff6666", "#990000", "#6666ff", "#000099", "#000000"])
plt.rc("axes", prop_cycle=default_cycler)
ds = xr.open_dataset("humidity_pdfs.nc")

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
    ax.set_ylim(0, 2.0)
    ax.set_xlim(0, 0.08)
    ax.set_ylabel("probability density")
    ax.set_title(titles[axnr])
    ax.grid(linestyle="dashed")
axis[3].set_xlabel("$q_l$ [kg/kg]")
axis[3].legend(["EPIC", "MPIC", "MONC"])
fig.tight_layout()
plt.savefig("compare_hl_pdf.png")
plt.savefig("compare_hl_pdf.jpeg")
plt.savefig("compare_hl_pdf.pdf")
