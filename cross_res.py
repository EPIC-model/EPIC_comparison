import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import colorcet as cc

plt.rcParams.update(
    {
        "figure.figsize": (8, 5),
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
# })

fig = plt.figure()

grid = ImageGrid(
    fig,
    111,
    nrows_ncols=(2, 4),
    aspect=True,
    axes_pad=(0.25, 0.3),
    direction="row",
    share_all=True,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="4%",
    cbar_pad=0.05,
)

files = [
    "epic_32_lim3r_t6.nc",
    "epic_64_lim3r_t6.nc",
    "epic_128_lim3r_t6.nc",
    "epic_256_lim3r_t6.nc",
    "mpic_32_lim3r_t6.nc",
    "mpic_64_lim3r_t6.nc",
    "mpic_128_lim3r_t6.nc",
    "mpic_256_lim3r_t6.nc",
]
labels = [
    "EPIC 32$^3$",
    "EPIC 64$^3$",
    "EPIC 128$^3$",
    "EPIC 256$^3$",
    "MPIC 32$^3$",
    "MPIC 64$^3$",
    "MPIC 128$^3$",
    "MPIC 256$^3$",
]

for findex, ax in enumerate(grid):
    # Iterating over the grid returns the Axes.
    ds = xr.open_dataset(files[findex])
    im = ax.pcolormesh(
        ds["X"], ds["Z"], ds["humidity"], vmin=0.0, vmax=0.08, cmap=cc.cm["rainbow4"]
    )
    ax.set_xlim(3.000, 5.000)
    ax.set_ylim(3.000, 5.000)
    ax.set_title(labels[findex])
    ax.set_aspect(1)
    if findex / 3 >= 1:
        ax.set_xlabel("x [-]")
    if findex % 3 == 0:
        ax.set_ylabel("y [-]")
cb = grid.cbar_axes[0].colorbar(im)
grid.cbar_axes[0].set_ylabel("humidity [-]")
grid.cbar_axes[0].yaxis.set_label_position("right")
plt.savefig("plotting_res.png")
plt.savefig("plotting_res.jpeg")
