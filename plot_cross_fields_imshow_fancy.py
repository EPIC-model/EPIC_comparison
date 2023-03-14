import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import colorcet as cc
from nc_reader import nc_reader
import matplotlib.colors as cls

# ~ from PIL import Image
from numba import njit

plt.rcParams.update(
    {
        "figure.figsize": (9, 4),
        "figure.dpi": 960 * 4,
        "font.family": "serif",
        "font.size": 11,
        "text.usetex": False,
    }
)

dpi = 960 * 4
cm = 1 / 2.54  # centimeters in inches
fig = plt.figure(figsize=(13 * cm, 5 * cm), dpi=dpi)
ax = plt.gca()
my_cmap = mpl.colors.LinearSegmentedColormap.from_list(
    "", ["white", *plt.cm.Blues(np.arange(255))]
)

ffile1 = "epic_zooms_0_lim3r_t6.nc"
sfile1 = "epic_zooms_0_sharp_t6.nc"


@njit
def replace_nans(aa, bb):
    for ii in range(np.shape(aa)[0]):
        for jj in range(np.shape(aa)[1]):
            if aa[ii, jj] != aa[ii, jj]:
                aa[ii, jj] = bb[ii, jj]
    return aa


def combine_ds(ds1, ds2):
    da2 = ds2["humidity"]
    da1 = ds1["humidity"].interp_like(da2)
    da2.values = replace_nans(da2.values, da1.values)
    return da2


# Iterating over the grid returns the Axes.
ds1 = xr.open_dataset(ffile1)
ds2 = xr.open_dataset(sfile1)


def do_zoom(ax, level, bounds, xlo, xhi, ylo, yhi, **kwargs):
    cmap = kwargs.pop("cmap", "Blues")
    ffile1 = "epic_zooms_" + str(level) + "_lim3r_t6.nc"
    sfile1 = "epic_zooms_" + str(level) + "_sharp_t6.nc"

    ds1 = xr.open_dataset(ffile1)
    ds2 = xr.open_dataset(sfile1)
    da = combine_ds(ds1, ds2)

    axins = ax.inset_axes(bounds, alpha=1.0)
    img = axins.imshow(
        da,
        vmin=0.0,
        vmax=0.0012,
        cmap=my_cmap,
        origin="lower",
        interpolation="lanczos",
        #        interpolation_stage='rgba',
        extent=[xlo, xhi, ylo, yhi],
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


left = 1250
right = 5250
top = 5000
bottom = 2000

ax.set_xlim([left, right])
ax.set_ylim([bottom, top])
ax.set_position([0.05, 0.10, 0.5, 0.8])

da = combine_ds(ds1, ds2)
ax.imshow(
    da,
    vmin=0.0,
    vmax=0.0012,
    cmap=my_cmap,
    origin="lower",
    interpolation="lanczos",
    #    interpolation_stage='rgba',
    extent=[left, right, bottom, top],
)

ax.set_aspect(1)

axins = do_zoom(
    ax=ax,
    level=1,
    bounds=[1.0, 0.47, 0.6, 0.6],
    xlo=3700,
    xhi=4485,
    ylo=3400,
    yhi=4185,
)

axins2 = do_zoom(
    ax=axins,
    level=2,
    bounds=[1.2, 0.06, 0.9, 0.9],
    xlo=3900,
    xhi=4292.5,
    ylo=3600,
    yhi=3992.5,
)

axins3 = do_zoom(
    ax=axins2,
    level=3,
    bounds=[0.1, -1.04, 0.9, 0.9],
    xlo=4000,
    xhi=4196.25,
    ylo=3700,
    yhi=3896.25,
)

axins4 = do_zoom(
    ax=axins3,
    level=4,
    bounds=[-1.4, 0.06, 0.9, 0.9],
    xlo=4050,
    xhi=4148.125,
    ylo=3750,
    yhi=3848.125,
)

ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
ax.tick_params(axis="y", which="both", bottom=False, top=False, labelbottom=False)

ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

plt.savefig("graphical_abstract.png", dpi=960)
plt.savefig("graphical_abstract.pdf", dpi=960)
plt.close()
