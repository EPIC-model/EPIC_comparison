import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import colorcet as cc
from utils import add_annotation

plt.rcParams.update(
    {
        "figure.figsize": (9, 4),
        "figure.dpi": 960,
        "font.family": "serif",
        "font.size": 11,
        "text.usetex": False,
    }
)

dpi = 960
cm = 1 / 2.54  # centimeters in inches
fig = plt.figure(figsize=(13 * cm, 5 * cm), dpi=dpi)
ax = plt.gca()
my_cmap = mpl.colors.LinearSegmentedColormap.from_list(
    "", ["white", *plt.cm.Blues(np.arange(255))]
)

ffile1 = "epic_256_lim3r_t6.nc"
ffile2 = "epic_256_sharp_t6.nc"

# Iterating over the grid returns the Axes.
ds1 = xr.open_dataset(ffile1)
ds2 = xr.open_dataset(ffile2)


def do_zoom(ax, bounds, xlo, xhi, ylo, yhi, **kwargs):
    cmap = kwargs.pop("cmap", "Blues")
    vmin = kwargs.pop("vmin", None)
    vmax = kwargs.pop("vmax", None)
    axins = ax.inset_axes(bounds, alpha=1.0)
    axins.pcolormesh(
        ds1["X"],
        ds1["Z"],
        ds1["humidity"],
        vmin=0.0,
        vmax=0.08,
        cmap=my_cmap,
        shading="gouraud",
        rasterized=True,
    )
    axins.pcolormesh(
        ds2["X"],
        ds2["Z"],
        ds2["humidity"],
        vmin=0.0,
        vmax=0.08,
        cmap=my_cmap,
        rasterized=True,
    )
    axins.set_xlabel("")
    axins.set_ylabel("")
    axins.set_xlim(xlo, xhi)
    axins.set_ylim(ylo, yhi)
    axins.set_xticks([])
    axins.set_xticklabels([])
    axins.set_yticks([])
    axins.set_yticklabels([])
    axins.set_aspect(1)
    ax.indicate_inset_zoom(axins, edgecolor="black", alpha=0.5)
    return axins


ax.pcolormesh(
    ds1["X"],
    ds1["Z"],
    ds1["humidity"],
    vmin=0.0,
    vmax=0.08,
    cmap=my_cmap,
    shading="gouraud",
    rasterized=True,
)
ax.pcolormesh(
    ds2["X"],
    ds2["Z"],
    ds2["humidity"],
    vmin=0.0,
    vmax=0.08,
    cmap=my_cmap,
    rasterized=True,
)

left = 1.250
right = 5.250
top = 5.000
bottom = 2.000

ax.set_xlim([left, right])
ax.set_ylim([bottom, top])
ax.set_aspect(1)

axins = do_zoom(
    ax=ax,
    bounds=[1.0, 0.47, 0.6, 0.6],
    xlo=3.700,
    xhi=4.485,
    ylo=3.400,
    yhi=4.185,
)

axins2 = do_zoom(
    ax=axins,
    bounds=[1.2, 0.06, 0.9, 0.9],
    xlo=3.900,
    xhi=4.2925,
    ylo=3.600,
    yhi=3.9925,
)

axins3 = do_zoom(
    ax=axins2,
    bounds=[0.1, -1.04, 0.9, 0.9],
    xlo=4.000,
    xhi=4.19625,
    ylo=3.700,
    yhi=3.89625,
)

axins4 = do_zoom(
    ax=axins3,
    bounds=[-1.4, 0.06, 0.9, 0.9],
    xlo=4.050,
    xhi=4.148125,
    ylo=3.750,
    yhi=3.848125,
)

ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)

ax.tick_params(axis="y", which="both", bottom=False, top=False, labelbottom=False)

ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

plt.tight_layout()
for axis in [ax]:
    aa = axis.get_position()
    axis.set_position([aa.x0 - 0.12, aa.y0 - 0.06, aa.x1 - 0.12, aa.y1 - 0.06])
for axis in [axins, axins2, axins3, axins4]:
    aa = axis.get_position()
    axis.set_position([aa.x0 - 0.06, aa.y0 - 0.06, aa.x1 - 0.06, aa.y1 - 0.06])
plt.savefig("cross_zoom.png")
plt.savefig("cross_zoom.jpeg")
plt.savefig("cross_zoom.pdf")
plt.close()
