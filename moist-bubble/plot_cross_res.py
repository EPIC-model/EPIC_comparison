import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import colorcet as cc
import matplotlib as mpl
from utils import add_annotation, setup_rcParams, remove_xticks, remove_yticks

my_cmap = mpl.colors.LinearSegmentedColormap.from_list(
    "", ["white", *plt.cm.Blues(np.arange(255))]
)

setup_rcParams()
fig = plt.figure(figsize=(8, 10))

grid = ImageGrid(
    fig,
    111,
    nrows_ncols=(4, 3),
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
    "epic_32_lim3r_t6.nc",
    "mpic_32_lim3r_t6.nc",
    "monc_32_t6.nc",
    "epic_64_lim3r_t6.nc",
    "mpic_64_lim3r_t6.nc",
    "monc_64_t6.nc",
    "epic_128_lim3r_t6.nc",
    "mpic_128_lim3r_t6.nc",
    "monc_128_t6.nc",
    "epic_256_lim3r_t6.nc",
    "mpic_256_lim3r_t6.nc",
    "monc_256_t6.nc",
]

for findex, ax in enumerate(grid):
    # Iterating over the grid returns the Axes.
    if findex < len(files):
        ds = xr.open_dataset(files[findex])
        im = ax.imshow(
            ds["humidity"],
            vmin=0.0,
            vmax=0.08,
            cmap=my_cmap,
            origin="lower",
            interpolation="bilinear",
            extent=[ds["X"][0], ds["X"][-1], ds["Z"][0], ds["Z"][-1]],
        )
        ax.set_xlim(2.900, 5.100)
        ax.set_ylim(2.900, 5.100)
        ax.set_xticks([3, 4, 5])
        ax.set_yticks([3, 4, 5])
        ax.set_aspect(1)
        if findex / 3 >= 3:
            ax.set_xlabel("$x$")
        else:
            remove_xticks(ax)
        if findex % 3 == 0:
            ax.set_ylabel("$z$")
        else:
            remove_yticks(ax)
    else:
        ax.axis("off")
cb = grid.cbar_axes[0].colorbar(im)
add_annotation(grid[0], "EPIC", [0.5, 1.1], ha="center")
add_annotation(grid[1], "MPIC", [0.5, 1.1], ha="center")
add_annotation(grid[2], "MONC", [0.5, 1.1], ha="center")
add_annotation(grid[0], "$32^3$", [-0.26, 0.5], va="center", ha="right")
add_annotation(grid[3], "$64^3$", [-0.26, 0.5], va="center", ha="right")
add_annotation(grid[6], "$128^3$", [-0.26, 0.5], va="center", ha="right")
add_annotation(grid[9], "$256^3$", [-0.26, 0.5], va="center", ha="right")
grid.cbar_axes[0].set_ylabel("humidity")
grid.cbar_axes[0].yaxis.set_label_position("right")
plt.savefig("cross_comp_hum_res.png", bbox_inches="tight")
plt.savefig("cross_comp_hum_res.jpeg", bbox_inches="tight")
plt.savefig("cross_comp_hum_res.pdf", bbox_inches="tight")
