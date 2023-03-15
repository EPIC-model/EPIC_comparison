import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt


def setup_rcParams():
    plt.rcParams.update(
        {
            "figure.dpi": 200,
            "font.family": "serif",
            "font.size": 11,
            "text.usetex": True,
            "text.latex.preamble": "\n".join(
                [
                    r"\usepackage{amsmath}",
                    r"\usepackage[utf8]{inputenc}",
                    r"\usepackage[T1]{fontenc}",
                    r"\usepackage{siunitx}",
                ]
            ),
        }
    )


def remove_xticks(ax):
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)


def remove_yticks(ax):
    ax.tick_params(axis="y", which="both", right=False, left=False)


def add_annotation(ax, label, xy, **kwargs):
    bbox = dict(boxstyle="round", facecolor="wheat", linewidth=0.5)
    ax.annotate(label, xy=xy, xycoords="axes fraction", bbox=bbox, **kwargs)


def check_file(input_file_name):
    file_root, file_ext = os.path.splitext(input_file_name)
    if not (file_ext == ".nc"):
        raise ValueError('Input filename must end on ".nc"')


# Write to netcdf file
def write_field_to_file(
    field, field_x, field_y, field_name, out_file, cross_type, mode
):
    ncfile = nc.Dataset(out_file, mode, format="NETCDF4", zlib=True)
    if cross_type == "xy":
        dimlist = ("Y", "X")
    elif cross_type == "xz":
        dimlist = ("Z", "X")
    elif cross_type == "yz":
        dimlist = ("Z", "Y")
    coords = [field_y, field_x]
    if mode == "w":
        for index, dim in enumerate(dimlist):
            ncfile.createDimension(dim, len(coords[index]))
            nc_dims = {}
            nc_dims[dim] = ncfile.createVariable(
                dim, "f8", (dim,), zlib=True, least_significant_digit=6
            )
            nc_dims[dim].units = "m"
            nc_dims[dim].axis = dim  # Optional
            nc_dims[dim].standard_name = "projection_" + dim.lower() + "_coordinate"
            nc_dims[dim].long_name = dim.lower() + "-coordinate"
            nc_dims[dim][:] = coords[index]
    elif not mode == "a":
        raise ValueError("netCDF access mode not implemented")
    nc_var = ncfile.createVariable(field_name, np.dtype("f4").char, dimlist, zlib=True)
    nc_var.units = "1"
    nc_var[:, :] = np.transpose(field[:, :])
    nc_var.long_name = field_name
    ncfile.close()
