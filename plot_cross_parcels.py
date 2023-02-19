import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import colorcet as cc
from nc_reader import nc_reader
import matplotlib.colors as cls

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

ffile1 = "epic_unsc_256_lim3r_t6.nc"
efile1 = "intersected_ellipses_step_12_from_moist_0000000012_parcels.nc"

# Iterating over the grid returns the Axes.
ds1 = xr.open_dataset(ffile1)


def do_zoom(ax, bounds, xlo, xhi, ylo, yhi, **kwargs):
    cmap = kwargs.pop("cmap", "Blues")
    axins = ax.inset_axes(bounds, alpha=1.0)
    axins.pcolormesh(
        ds1["X"],
        ds1["Z"],
        ds1["humidity"],
        vmin=0.0,
        vmax=0.0012,
        cmap=my_cmap,
        shading="gouraud",
        rasterized=True,
        snap=True,
    )
    encr = nc_reader()
    encr.open(efile1)
    dx = 6.28 / 256
    x_pos = encr.get_dataset(0, "x_position")
    z_pos = encr.get_dataset(0, "z_position")
    ind = np.argwhere(
        (x_pos >= xlo - dx)
        & (x_pos <= xhi + dx)
        & (z_pos >= ylo - dx)
        & (z_pos <= yhi + dx)
    )
    ind = ind.squeeze()
    hum = encr.get_dataset(0, "humidity", indices=ind)
    sort = np.argsort(hum)
    ell = encr.get_ellipses(step=0, indices=ind[sort])
    axins.add_collection(ell)
    ell.set_offset_transform(axins.transData)
    ell.set_clip_box(axins.bbox)
    ell.set_alpha(1.0)
    norm = cls.Normalize(vmin=0.0, vmax=0.0012)
    ell.set_facecolor(my_cmap(norm(hum[sort])))
    encr.close()
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
    vmax=0.0012,
    cmap=my_cmap,
    shading="gouraud",
    rasterized=True,
    snap=True,
)

left = 1250
right = 5250
top = 5000
bottom = 2000

encr = nc_reader()
encr.open(efile1)
dx = 6280 / 256
x_pos = encr.get_dataset(0, "x_position")
z_pos = encr.get_dataset(0, "z_position")
ind = np.argwhere(
    (x_pos >= left - dx)
    & (x_pos <= right + dx)
    & (z_pos >= bottom - dx)
    & (z_pos <= top + dx)
)
ind = ind.squeeze()
hum = encr.get_dataset(0, "humidity", indices=ind)
sort = np.argsort(hum)
ell = encr.get_ellipses(step=0, indices=ind[sort])
ax.add_collection(ell)
ell.set_offset_transform(ax.transData)
ell.set_clip_box(ax.bbox)
ell.set_alpha(1.0)
norm = cls.Normalize(vmin=0.0, vmax=0.0012)
ell.set_facecolor(my_cmap(norm(hum[sort])))
encr.close()


ax.set_xlim([left, right])
ax.set_ylim([bottom, top])
ax.set_position([0.05, 0.10, 0.5, 0.8])
ax.set_aspect(1)

axins = do_zoom(
    ax=ax,
    bounds=[1.0, 0.47, 0.6, 0.6],
    xlo=3700,
    xhi=4485,
    ylo=3400,
    yhi=4185,
)

axins2 = do_zoom(
    ax=axins,
    bounds=[1.2, 0.06, 0.9, 0.9],
    xlo=3900,
    xhi=4292.5,
    ylo=3600,
    yhi=3992.5,
)

axins3 = do_zoom(
    ax=axins2,
    bounds=[0.1, -1.04, 0.9, 0.9],
    xlo=4000,
    xhi=4196.25,
    ylo=3700,
    yhi=3896.25,
)

axins4 = do_zoom(
    ax=axins3,
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

plt.savefig("cross_zoom.png")
plt.savefig("cross_zoom.jpeg")
plt.savefig("cross_zoom.pdf")
plt.close()
