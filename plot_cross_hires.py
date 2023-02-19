import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import colorcet as cc

plt.rcParams.update(
    {
        "figure.figsize": (5, 5),
        "figure.dpi": 300,
        "font.family": "serif",
        "font.size": 11,
        "text.usetex": False,
    }
)
#    "text.usetex": True,
#    'text.latex.preamble': "\n".join([
#        r"\usepackage{amsmath}",
#        r"\usepackage[utf8]{inputenc}",
#        r"\usepackage[T1]{fontenc}",
#        r"\usepackage{siunitx}",
#        ])
# })


ffile1 = "epic_256_lim3r_t6.nc"
ffile2 = "epic_256_sharp_t6.nc"

xlos = [0.0, 3.5, 3.7, 4.1]
xhis = [6.28, 4.8, 4.5, 4.4]
zlos = [0.0, 3.2, 3.4, 3.7]
zhis = [6.28, 4.5, 4.3, 4.0]

# Iterating over the grid returns the Axes.
ds1 = xr.open_dataset(ffile1)
ds2 = xr.open_dataset(ffile2)
for findex in range(len(xlos)):
    fig = plt.figure()
    im = plt.pcolormesh(
        ds1["X"],
        ds1["Z"],
        ds1["humidity"],
        vmin=0.0,
        vmax=0.08,
        cmap=cc.cm["rainbow4"],
        shading="gouraud",
        rasterized=True,
    )
    im2 = plt.pcolormesh(
        ds2["X"],
        ds2["Z"],
        ds2["humidity"],
        vmin=0.0,
        vmax=0.08,
        cmap=cc.cm["rainbow4"],
        rasterized=True,
    )
    plt.xlim(xlos[findex], xhis[findex])
    plt.ylim(zlos[findex], zhis[findex])
    plt.gca().set_aspect(1)
    plt.xlabel("x [-]")
    plt.ylabel("y [-]")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("cross_zoom_" + str(findex) + ".png")
    plt.savefig("cross_zoom_" + str(findex) + ".jpeg")
    plt.savefig("cross_zoom_" + str(findex) + ".pdf")
    plt.close()
