import matplotlib.pyplot as plt
import xarray as xr
from cycler import cycler
import numpy as np
from matplotlib import ticker

plt.rcParams.update(
    {
        "figure.figsize": (4, 8),
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

fig, axis = plt.subplots(3, 1)
titles = ["EPIC ref error", "MPIC ref error", "MONC ref error"]


def calc_err(sim, ref):
    return np.sqrt(np.sum((sim - ref) * (sim - ref)) / len(sim))


epic_ref = ds["hist_epic_" + str(RESOLUTIONS[-1])][1:]
monc_ref = ds["hist_monc_" + str(RESOLUTIONS[-1])][1:]
mpic_ref = ds["hist_mpic_" + str(RESOLUTIONS[-1])][1:]

for resolution in RESOLUTIONS:
    epic = ds["hist_epic_" + str(resolution)][1:]
    mpic = ds["hist_mpic_" + str(resolution)][1:]
    monc = ds["hist_monc_" + str(resolution)][1:]
    axis[0].scatter(
        resolution, calc_err(epic, epic_ref), edgecolors="#ff6666", facecolors="#ff6666"
    )
    axis[0].scatter(
        resolution, calc_err(mpic, epic_ref), edgecolors="#990000", facecolors="none"
    )
    axis[0].scatter(
        resolution, calc_err(monc, epic_ref), edgecolors="#6666ff", facecolors="none"
    )
    axis[1].scatter(
        resolution, calc_err(epic, mpic_ref), edgecolors="#ff6666", facecolors="#ff6666"
    )
    axis[1].scatter(
        resolution, calc_err(mpic, mpic_ref), edgecolors="#990000", facecolors="none"
    )
    axis[1].scatter(
        resolution, calc_err(monc, mpic_ref), edgecolors="#6666ff", facecolors="none"
    )
    axis[2].scatter(
        resolution, calc_err(epic, monc_ref), edgecolors="#ff6666", facecolors="#ff6666"
    )
    axis[2].scatter(
        resolution, calc_err(mpic, monc_ref), edgecolors="#990000", facecolors="none"
    )
    axis[2].scatter(
        resolution, calc_err(monc, monc_ref), edgecolors="#6666ff", facecolors="none"
    )
for axnr, ax in enumerate(axis):
    ax.set_xscale("log")
    ax.set_xticks([32, 64, 128, 256])
    ax.minorticks_off()
    ax.set_title(titles[axnr])
    ax.grid(linestyle="dashed")
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:d}"))
    ax.set_ylim(0, 0.5)
axis[2].set_xlabel("grid points")
axis[0].legend(["EPIC", "MPIC", "MONC"])
fig.tight_layout()
plt.savefig("err_hl_pdf.png")
plt.savefig("err_hl_pdf.jpeg")
plt.savefig("err_hl_pdf.pdf")
