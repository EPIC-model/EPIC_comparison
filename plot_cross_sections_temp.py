import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import colorcet as cc
from utils import add_annotation, setup_rcParams, remove_xticks, remove_yticks
import matplotlib as mpl
from numba import njit

setup_rcParams()
fig = plt.figure(figsize=(8, 5))


@njit
def replace_nans_by_zero(aa):
    for ii in range(np.shape(aa)[0]):
        for jj in range(np.shape(aa)[1]):
            if aa[ii, jj] != aa[ii, jj]:
                aa[ii, jj] = 0.0
    return aa


grid = ImageGrid(
    fig,
    111,
    nrows_ncols=(2, 3),
    aspect=True,
    axes_pad=(0.1, 0.1),
    direction="row",
    share_all=True,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="4%",
    cbar_pad=0.1,
)

files = [
    "epic_32_sharp_t6.nc",
    "epic_32_lim3r_t6.nc",
    "epic_32_lim5r_t6.nc",
    "mpic_32_sharp_t6.nc",
    "mpic_32_lim3r_t6.nc",
    "mpic_32_lim5r_t6.nc",
]

my_cmap = mpl.colors.LinearSegmentedColormap.from_list(
    "", ["white", *plt.cm.Blues(np.arange(255))]
)

for findex, ax in enumerate(grid):
    # Iterating over the grid returns the Axes.
    if findex < len(files):
        ds = xr.open_dataset(files[findex])
        ds["humidity"].values = replace_nans_by_zero(ds["humidity"].values)
        im = ax.imshow(
            ds["humidity"],
            vmin=0.0,
            vmax=0.09,
            cmap=my_cmap,
            origin="lower",
            interpolation="lanczos",
            extent=[ds["X"][0], ds["X"][-1], ds["Z"][0], ds["Z"][-1]],
        )
        ax.set_xlim(2.900, 5.100)
        ax.set_ylim(2.900, 5.100)
        ax.set_aspect(1)
        if findex / 3 >= 1:
            ax.set_xlabel("x (-)")
            ax.set_xticks([3, 4, 5])
        else:
            remove_xticks(ax)
        if findex % 3 == 0:
            ax.set_ylabel("z (-)")
            ax.set_yticks([3, 4, 5])
        else:
            remove_yticks(ax)
    else:
        ax.axis("off")
add_annotation(grid[0], "Parcels", [0.5, 1.2], ha="center")
add_annotation(grid[1], "Supergaussian (3x)", [0.5, 1.2], ha="center")
add_annotation(grid[2], "Supergaussian (5x)", [0.5, 1.2], ha="center")
add_annotation(grid[0], "EPIC", [-0.5, 0.5], va="center")
add_annotation(grid[3], "MPIC", [-0.5, 0.5], va="center")
cb = grid.cbar_axes[0].colorbar(im)
grid.cbar_axes[0].set_ylabel("humidity (-)")
grid.cbar_axes[0].yaxis.set_label_position("right")
plt.savefig("cross_comp_hum.png", bbox_inches="tight")
plt.savefig("cross_comp_hum.jpeg", bbox_inches="tight")
plt.savefig("cross_comp_hum.pdf", bbox_inches="tight")
