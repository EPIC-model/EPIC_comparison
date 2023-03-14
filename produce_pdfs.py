import numpy as np
import netCDF4 as nc
import xarray as xr
from glob import glob
from numba import jit
import scipy

MAIN_DIR = "/work/e710/e710/shared/epic_comparison"
RESOLUTIONS = [32, 64, 128, 256]
xz_loc = 3.140
MPIC_RES_TSTEP_DICT = {32: "0099", 64: "0193", 128: "0475", 256: "0941"}

# Upscale field, using "constant" interpolation at top and bottom boundary
@jit
def upscale_monc_field_3(monc_field):
    nt = np.shape(monc_field)[0]
    nx = np.shape(monc_field)[1]
    ny = np.shape(monc_field)[2]
    nz = np.shape(monc_field)[3]
    monc_field_scaled = np.zeros((nt, 3 * nx, 3 * ny, 3 * nz))
    f23 = 2.0 / 3.0
    f13 = 1.0 / 3.0
    # perform bilinear interpolation
    for tt in range(nt):
        # Bottom 3 cells
        monc_hfield_scaled = scipy.ndimage.zoom(
            monc_field[tt, :, :, 0], 3, order=1, mode="grid-wrap"
        )
        monc_hfield_scaled_plus = scipy.ndimage.zoom(
            monc_field[tt, :, :, 1], 3, order=1, mode="grid-wrap"
        )
        monc_field_scaled[tt, :, :, 0] = monc_hfield_scaled
        monc_field_scaled[tt, :, :, 1] = monc_hfield_scaled
        monc_field_scaled[tt, :, :, 2] = (
            f23 * monc_hfield_scaled + f13 * monc_hfield_scaled_plus
        )
        # Middle layers
        for kk in range(1, nz - 1):
            monc_hfield_scaled_minus = monc_hfield_scaled
            monc_hfield_scaled = monc_hfield_scaled_plus
            monc_hfield_scaled_plus = scipy.ndimage.zoom(
                monc_field[tt, :, :, kk + 1], 3, order=1, mode="grid-wrap"
            )
            monc_field_scaled[tt, :, :, 3 * kk] = (
                f13 * monc_hfield_scaled_minus + f23 * monc_hfield_scaled
            )
            monc_field_scaled[tt, :, :, 3 * kk + 1] = monc_hfield_scaled
            monc_field_scaled[tt, :, :, 3 * kk + 2] = (
                f23 * monc_hfield_scaled + f13 * monc_hfield_scaled_plus
            )
        # Top layer
        monc_hfield_scaled_minus = monc_hfield_scaled
        monc_hfield_scaled = monc_hfield_scaled_plus
        monc_field_scaled[tt, :, :, 3 * (nz - 1)] = (
            f13 * monc_hfield_scaled_minus + f23 * monc_hfield_scaled
        )
        monc_field_scaled[tt, :, :, 3 * (nz - 1) + 1] = monc_hfield_scaled
        monc_field_scaled[tt, :, :, 3 * (nz - 1) + 2] = monc_hfield_scaled
    return monc_field_scaled


@jit
def upscale_monc_field_5(monc_field):
    nt = np.shape(monc_field)[0]
    nx = np.shape(monc_field)[1]
    ny = np.shape(monc_field)[2]
    nz = np.shape(monc_field)[3]
    monc_field_scaled = np.zeros((nt, 5 * nx, 5 * ny, 5 * nz))
    f15 = 0.2
    f25 = 0.4
    f35 = 0.6
    f45 = 0.8
    # perform bilinear interpolation
    for tt in range(nt):
        # Bottom 3 cells
        monc_hfield_scaled = scipy.ndimage.zoom(
            monc_field[tt, :, :, 0], 5, order=1, mode="grid-wrap"
        )
        monc_hfield_scaled_plus = scipy.ndimage.zoom(
            monc_field[tt, :, :, 1], 5, order=1, mode="grid-wrap"
        )
        monc_field_scaled[tt, :, :, 0] = monc_hfield_scaled
        monc_field_scaled[tt, :, :, 1] = monc_hfield_scaled
        monc_field_scaled[tt, :, :, 2] = monc_hfield_scaled
        monc_field_scaled[tt, :, :, 3] = (
            f45 * monc_hfield_scaled + f15 * monc_hfield_scaled_plus
        )
        monc_field_scaled[tt, :, :, 4] = (
            f35 * monc_hfield_scaled + f25 * monc_hfield_scaled_plus
        )
        # Middle layers
        for kk in range(1, nz - 1):
            monc_hfield_scaled_minus = monc_hfield_scaled
            monc_hfield_scaled = monc_hfield_scaled_plus
            monc_hfield_scaled_plus = scipy.ndimage.zoom(
                monc_field[tt, :, :, kk + 1], 5, order=1, mode="grid-wrap"
            )
            monc_field_scaled[tt, :, :, 5 * kk] = (
                f25 * monc_hfield_scaled_minus + f35 * monc_hfield_scaled
            )
            monc_field_scaled[tt, :, :, 5 * kk + 1] = (
                f15 * monc_hfield_scaled_minus + f45 * monc_hfield_scaled
            )
            monc_field_scaled[tt, :, :, 5 * kk + 2] = monc_hfield_scaled
            monc_field_scaled[tt, :, :, 5 * kk + 3] = (
                f45 * monc_hfield_scaled + f15 * monc_hfield_scaled_plus
            )
            monc_field_scaled[tt, :, :, 5 * kk + 4] = (
                f35 * monc_hfield_scaled + f25 * monc_hfield_scaled_plus
            )
        # Top layer
        monc_hfield_scaled_minus = monc_hfield_scaled
        monc_hfield_scaled = monc_hfield_scaled_plus
        monc_field_scaled[tt, :, :, 5 * (nz - 1)] = (
            f25 * monc_hfield_scaled_minus + f35 * monc_hfield_scaled
        )
        monc_field_scaled[tt, :, :, 5 * (nz - 1) + 1] = (
            f15 * monc_hfield_scaled_minus + f45 * monc_hfield_scaled
        )
        monc_field_scaled[tt, :, :, 5 * (nz - 1) + 2] = monc_hfield_scaled
        monc_field_scaled[tt, :, :, 5 * (nz - 1) + 3] = monc_hfield_scaled
        monc_field_scaled[tt, :, :, 5 * (nz - 1) + 4] = monc_hfield_scaled
    return monc_field_scaled


@jit
def get_monc_enstrophy(uu, vv, ww, dx, dy, dz):
    nt = np.shape(uu)[0]
    nx = np.shape(uu)[1]
    ny = np.shape(uu)[2]
    nz = np.shape(uu)[3]
    oo1 = np.zeros((nt, nx, ny, nz))
    oo2 = np.zeros((nt, nx, ny, nz))
    oo3 = np.zeros((nt, nx, ny, nz))
    oo = np.zeros((nt, nx, ny, nz - 1))
    for tt in range(nt):
        for ii in range(nx):
            for jj in range(ny):
                for kk in range(nz - 1):
                    ip = (ii + 1) % nx
                    jp = (jj + 1) % ny
                    kp = kk + 1
                    dwdy = (ww[tt, ii, jp, kk] - ww[tt, ii, jj, kk]) / dy
                    dvdz = (vv[tt, ii, jj, kp] - vv[tt, ii, jj, kk]) / dz
                    oo1[tt, ii, jj, kk] = (dwdy - dvdz) * (dwdy - dvdz)
                    dudz = (uu[tt, ii, jj, kp] - uu[tt, ii, jj, kk]) / dz
                    dwdx = (ww[tt, ip, jj, kk] - ww[tt, ii, jj, kk]) / dx
                    oo2[tt, ii, jj, kk] = (dudz - dwdx) * (dudz - dwdx)
    for tt in range(nt):
        for ii in range(nx):
            for jj in range(ny):
                for kk in range(nz):
                    ip = (ii + 1) % nx
                    jp = (jj + 1) % ny
                    dvdx = (vv[tt, ip, jj, kk] - vv[tt, ii, jj, kk]) / dx
                    dudy = (uu[tt, ii, jp, kk] - uu[tt, ii, jj, kk]) / dy
                    oo3[tt, ii, jj, kk] = (dvdx - dudy) * (dvdx - dudy)
    for tt in range(nt):
        for ii in range(nx):
            for jj in range(ny):
                for kk in range(nz - 1):
                    ip = (ii + 1) % nx
                    jp = (jj + 1) % ny
                    kp = kk + 1
                    oo[tt, ip, jp, kk] = (
                        0.5 * (oo1[tt, ii, jj, kk] + oo1[tt, ip, jj, kk])
                        + 0.5 * (oo2[tt, ii, jj, kk] + oo2[tt, ii, jp, kk])
                        + 0.5 * (oo3[tt, ii, jj, kk] + oo3[tt, ii, jj, kp])
                    )
    return oo


# Set up bin space for liquid water

nbins = 50
bin_edges = np.linspace(0.0, 0.10, nbins)
bin_centres = 0.5 * (bin_edges[1:] + bin_edges[:-1])

bin_edges_coord = {"bin_edges": ("bin_edges", bin_edges, {"long_name": "Bin edges"})}
bin_centres_coord = {
    "bin_centres": ("bin_centres", bin_centres, {"long_name": "Bin centres"})
}

ds = xr.Dataset(
    coords={
        **bin_edges_coord,
        **bin_centres_coord,
    }
)

for resolution in RESOLUTIONS:
    # Set up file paths
    monc_input_file_name = (
        MAIN_DIR + "/monc_mpic_smooth_" + str(resolution) + "/mpic_diagnostic_3d_6.nc"
    )
    epic_input_file_name = (
        MAIN_DIR + "/epic_rev_" + str(resolution) + "/moist_0000000012_parcels.nc"
    )
    mpic_parcel_glob = (
        MAIN_DIR
        + "/pmpic_smooth_bubble_"
        + str(resolution)
        + "/parcels_?????_0"
        + MPIC_RES_TSTEP_DICT[resolution]
        + ".nc"
    )
    ds_nc = nc.Dataset(epic_input_file_name)
    vol = ds_nc.variables["volume"][0, :]
    h = ds_nc.variables["humidity"][0, :]
    z = ds_nc.variables["z_position"][0, :]
    len_condense = ds_nc.groups["physical_quantities"].scale_height
    q_scale = ds_nc.groups[
        "physical_quantities"
    ].saturation_specific_humidity_at_ground_level
    hl = h / q_scale - np.exp(-z / len_condense)
    hl = hl * (hl > 0.0)
    hist_epic, bins_epic = np.histogram(hl, weights=vol, density=True, bins=bin_edges)
    ds["hist_epic_" + str(resolution)] = (
        ("bin_centres"),
        hist_epic,
    )
    ds["hist_epic_" + str(resolution)].attrs = {
        "long_name": "EPIC histogram for $" + str(resolution) + "^3$ points",
        "units": "-",
    }
    hl_all = np.array([])
    vol_all = np.array([])
    for input_file_name in glob(mpic_parcel_glob):
        ds_nc = nc.Dataset(input_file_name)
        vol = ds_nc.variables["vol"][:]
        h = ds_nc.variables["h"][:]
        z = ds_nc.variables["z"][:]
        hl = h - np.exp(-z)
        hl = hl * (hl > 0.0)
        hl_all = np.hstack((hl_all, hl))
        vol_all = np.hstack((vol_all, vol))
    hist_mpic, bins_mpic = np.histogram(
        hl_all, weights=vol_all, density=True, bins=bin_edges
    )
    ds["hist_mpic_" + str(resolution)] = (
        ("bin_centres"),
        hist_mpic,
    )
    ds["hist_mpic_" + str(resolution)].attrs = {
        "long_name": "MPIC histogram for $" + str(resolution) + "$^3$ points",
        "units": "-",
    }
    ds_nc = nc.Dataset(monc_input_file_name)
    if resolution < 50:
        hh = ds_nc.variables["q_vapour"][:, :, :, 1:]
        hh5 = upscale_monc_field_5(hh)
        zn = ds_nc.variables["zn"][1:]
        zn5 = 0.2 * zn[0] + 0.2 * (zn[1] - zn[0]) * np.arange(np.shape(hh5)[3])
        hl = hh5 - np.exp(-zn5)[None, None, None, :]
        hl = hl * (hl > 0.0)
    elif resolution < 100:
        hh = ds_nc.variables["q_vapour"][:, :, :, 1:]
        hh3 = upscale_monc_field_3(hh)
        zn = ds_nc.variables["zn"][1:]
        zn3 = zn[0] / 3.0 + ((zn[1] - zn[0]) / 3.0) * np.arange(np.shape(hh3)[3])
        hl = hh3 - np.exp(-zn3)[None, None, None, :]
        hl = hl * (hl > 0.0)
    else:
        hl = ds_nc.variables["q_cloud_liquid_mass"][:, :, :, 1:]
    hist_monc, bins_monc = np.histogram(hl, density=True, bins=bin_edges)
    ds["hist_monc_" + str(resolution)] = (
        ("bin_centres"),
        hist_monc,
    )
    ds["hist_monc_" + str(resolution)].attrs = {
        "long_name": "MONC histogram for $" + str(resolution) + "$^3$ points",
        "units": "-",
    }
    hl_unsampled = ds_nc.variables["q_cloud_liquid_mass"][:, :, :, 1:]
    hist_monc_uns, bins_monc_uns = np.histogram(
        hl_unsampled, density=True, bins=bin_edges
    )
    ds["hist_monc_uns_" + str(resolution)] = (
        ("bin_centres"),
        hist_monc_uns,
    )
    ds["hist_monc_uns_" + str(resolution)].attrs = {
        "long_name": "MONC unsampled histogram for $" + str(resolution) + "$^3$ points",
        "units": "-",
    }
ds.to_netcdf("humidity_pdfs.nc")

# Set up bin space for vorticity
nbins = 400
bin_edges = np.array([0.000001 * 1.05**x for x in range(nbins)])
bin_centres = 0.5 * (bin_edges[1:] + bin_edges[:-1])

bin_edges_coord = {"bin_edges": ("bin_edges", bin_edges, {"long_name": "Bin edges"})}
bin_centres_coord = {
    "bin_centres": ("bin_centres", bin_centres, {"long_name": "Bin centres"})
}

ds = xr.Dataset(
    coords={
        **bin_edges_coord,
        **bin_centres_coord,
    }
)


for resolution in RESOLUTIONS:
    # Set up file paths
    monc_input_file_name = (
        MAIN_DIR + "/monc_mpic_smooth_" + str(resolution) + "/mpic_diagnostic_3d_6.nc"
    )
    epic_input_file_name = (
        MAIN_DIR + "/epic_rev_" + str(resolution) + "/moist_0000000012_parcels.nc"
    )
    mpic_parcel_glob = (
        MAIN_DIR
        + "/pmpic_smooth_bubble_"
        + str(resolution)
        + "/parcels_?????_0"
        + MPIC_RES_TSTEP_DICT[resolution]
        + ".nc"
    )
    ds_nc = nc.Dataset(epic_input_file_name)
    epic_t_scale = 100.0
    vol = ds_nc.variables["volume"][0, :]
    p = ds_nc.variables["x_vorticity"][0, :]
    q = ds_nc.variables["y_vorticity"][0, :]
    r = ds_nc.variables["z_vorticity"][0, :]
    oo = np.sqrt(p * p + q * q + r * r)
    oo = oo * epic_t_scale
    hist_epic, bins_epic = np.histogram(oo, weights=vol, density=True, bins=bin_edges)
    ds["hist_epic_" + str(resolution)] = (
        ("bin_centres"),
        hist_epic,
    )
    ds["hist_epic_" + str(resolution)].attrs = {
        "long_name": "EPIC histogram for $" + str(resolution) + "^3$ points",
        "units": "-",
    }
    oo_all = np.array([])
    vol_all = np.array([])
    for input_file_name in glob(mpic_parcel_glob):
        ds_nc = nc.Dataset(input_file_name)
        vol = ds_nc.variables["vol"][:]
        p = ds_nc.variables["p"][:]
        q = ds_nc.variables["q"][:]
        r = ds_nc.variables["r"][:]
        oo_all = np.hstack((oo_all, np.sqrt(p * p + q * q + r * r)))
        vol_all = np.hstack((vol_all, vol))
    hist_mpic, bins_mpic = np.histogram(
        oo_all, weights=vol_all, density=True, bins=bin_edges
    )
    ds["hist_mpic_" + str(resolution)] = (
        ("bin_centres"),
        hist_mpic,
    )
    ds["hist_mpic_" + str(resolution)].attrs = {
        "long_name": "MPIC histogram for $" + str(resolution) + "$^3$ points",
        "units": "-",
    }
    ds_nc = nc.Dataset(monc_input_file_name)
    uu = ds_nc.variables["u"][:, :, :, :]
    vv = ds_nc.variables["v"][:, :, :, :]
    ww = ds_nc.variables["w"][:, :, :, :]
    dx = 6.28 / 32.0
    dy = 6.28 / 32.0
    dz = 6.28 / 32.0
    oo = get_monc_enstrophy(uu, vv, ww, dx, dy, dz)
    ds_oo_t = {
        "t": (
            "t",
            range(np.shape(ds_nc.variables["u"][:, :, :, :])[0]),
            {"long_name": "time"},
        )
    }
    ds_oo_x = {
        "x": (
            "x",
            range(np.shape(ds_nc.variables["u"][:, :, :, :])[1]),
            {"long_name": "x"},
        )
    }
    ds_oo_y = {
        "y": (
            "y",
            range(np.shape(ds_nc.variables["u"][:, :, :, :])[2]),
            {"long_name": "y"},
        )
    }
    ds_oo_z = {"z": ("z", ds_nc.variables["z"][:-1], {"long_name": "z"})}
    ds_oo = xr.Dataset(
        coords={
            **ds_oo_t,
            **ds_oo_x,
            **ds_oo_y,
            **ds_oo_z,
        }
    )
    # ds_oo["oo"] = (("t","x","y","z"),oo)
    hist_monc, bins_monc = np.histogram(np.sqrt(oo), density=True, bins=bin_edges)
    ds["hist_monc_" + str(resolution)] = (
        ("bin_centres"),
        hist_monc,
    )
    ds["hist_monc_" + str(resolution)].attrs = {
        "long_name": "MONC histogram for $" + str(resolution) + "$^3$ points",
        "units": "-",
    }
ds.to_netcdf("vorticity_pdfs.nc")
