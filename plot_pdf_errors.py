import matplotlib.pyplot as plt
import xarray as xr
from cycler import cycler
import numpy as np
from matplotlib import ticker
from utils import add_annotation, setup_rcParams

setup_rcParams()

RESOLUTIONS = [32, 64, 128, 256]
default_cycler = cycler(color=["#ff6666", "#990000", "#6666ff", "#000099", "#000000"])
plt.rc("axes", prop_cycle=default_cycler)
ds = xr.open_dataset("humidity_pdfs.nc")

fig, axis = plt.subplots(1, 3, figsize=(8, 4))


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
    ax.grid(linestyle="dashed")
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("${x:d}^3$"))
    ax.set_ylim(0, 0.55)
    ax.set_xlabel("grid points")
axis[0].set_ylabel("RMS error")
axis[0].legend(["EPIC", "MPIC", "MONC"],loc=(0.4,0.65))
add_annotation(
    axis[0], "Reference: EPIC, $" + str(RESOLUTIONS[-1]) + "^3$", [0.5, 0.9], ha="center"
)
add_annotation(
    axis[1], "Reference: MPIC, $" + str(RESOLUTIONS[-1]) + "^3$", [0.5, 0.9], ha="center"
)
add_annotation(
    axis[2], "Reference: MONC, $" + str(RESOLUTIONS[-1]) + "^3$", [0.5, 0.9], ha="center"
)
fig.tight_layout()
plt.savefig("err_hl_pdf.png", bbox_inches="tight")
plt.savefig("err_hl_pdf.jpeg", bbox_inches="tight")
plt.savefig("err_hl_pdf.pdf", bbox_inches="tight")
