from glob import glob
from numba import njit, prange, set_num_threads
import netCDF4 as nc
import numpy as np
from utils import check_file

NTHREADS = 4

set_num_threads(NTHREADS)


def get_mpic_parcel_coordinates(input_file_name):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    x = ds_nc.variables["x"][:]
    y = ds_nc.variables["y"][:]
    z = ds_nc.variables["z"][:]
    ds_nc.close()
    return x, y, z


def get_mpic_parcel_volume(input_file_name):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    vol = ds_nc.variables["vol"][:]
    ds_nc.close()
    return vol


def get_mpic_parcel_buoyancy(input_file_name):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    buoyancy = ds_nc.variables["b"][:]
    ds_nc.close()
    return buoyancy


def get_mpic_parcel_humidity(input_file_name):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    humidity = ds_nc.variables["h"][:]
    ds_nc.close()
    return humidity


def get_mpic_parcel_vorticity(input_file_name):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    x_vorticity = ds_nc.variables["p"][:]
    y_vorticity = ds_nc.variables["q"][:]
    z_vorticity = ds_nc.variables["r"][:]
    ds_nc.close()
    return x_vorticity, y_vorticity, z_vorticity


def get_mpic_projection_coordinates(grids_file, refinement_factor):
    check_file(grids_file)
    ds_nc = nc.Dataset(grids_file)
    # initialisation
    nx = ds_nc.dimensions["x"].size
    ny = ds_nc.dimensions["y"].size
    nz = ds_nc.dimensions["z"].size - 1
    # projection grid
    nxproj = nx * refinement_factor
    nyproj = ny * refinement_factor
    nzproj = nz * refinement_factor + 1
    dx = ds_nc["x"][1] - ds_nc["x"][0]
    dy = ds_nc["y"][1] - ds_nc["y"][0]
    dz = ds_nc["z"][1] - ds_nc["z"][0]
    origin = [
        ds_nc["x"][0] - 0.5 * dx,
        ds_nc["y"][0] - 0.5 * dy,
        ds_nc["z"][0] - 0.5 * dz,
    ]
    extent = [nx * dx, ny * dy, nz * dz]
    xp = np.linspace(origin[0], origin[0] + extent[0], nxproj + 1)[0:nxproj]
    yp = np.linspace(origin[1], origin[1] + extent[1], nyproj + 1)[0:nyproj]
    zp = np.linspace(origin[2], origin[2] + extent[2], nzproj)
    ds_nc.close()
    return xp, yp, zp


def get_mpic_extent(grids_file):
    check_file(grids_file)
    ds_nc = nc.Dataset(grids_file)
    # initialisation
    nx = ds_nc.dimensions["x"].size
    ny = ds_nc.dimensions["y"].size
    nz = ds_nc.dimensions["z"].size - 1
    dx = ds_nc["x"][1] - ds_nc["x"][0]
    dy = ds_nc["y"][1] - ds_nc["y"][0]
    dz = ds_nc["z"][1] - ds_nc["z"][0]
    extent = [nx * dx, ny * dy, nz * dz]
    ds_nc.close()
    return extent


def produce_mpic_cross_section(
    parcel_glob,
    grid_file,
    loc=0.0,
    cross_type="xz",
    kernel="third",
    refinement_factor=8,
    r_limit_fac=None,
):
    if cross_type not in ["xz", "xy", "yz"]:
        raise ValueError("cross-section type not value")
    if kernel not in ["third", "sharp"]:
        raise ValueError("kernel type not valid")
    if r_limit_fac is None:
        if kernel == "third":
            r_limit_fac = 5.0
        elif kernel == "sharp":
            r_limit_fac = 1.0
    xp, yp, zp = get_mpic_projection_coordinates(grid_file, refinement_factor)
    extent = get_mpic_extent(grid_file)
    nxproj = len(xp)
    nyproj = len(yp)
    nzproj = len(zp)
    if cross_type == "xy":
        cross_x_len = nxproj
        cross_y_len = nyproj
    elif cross_type == "xz":
        cross_x_len = nxproj
        cross_y_len = nzproj
    elif cross_type == "yz":
        cross_x_len = nyproj
        cross_y_len = nzproj
    cross_field_sum = np.zeros((cross_x_len, cross_y_len), np.double)
    cross_volg_sum = np.zeros((cross_x_len, cross_y_len), np.double)
    input_files = glob(parcel_glob)
    for input_file_name in input_files:
        x, y, z = get_mpic_parcel_coordinates(input_file_name)
        vol = get_mpic_parcel_volume(input_file_name)
        h = get_mpic_parcel_humidity(input_file_name)
        cross_field_thread, cross_volg_thread = mpic_cross_section_inner(
            cross_type,
            loc,
            x,
            y,
            z,
            vol,
            h,
            xp,
            yp,
            zp,
            extent,
            kernel,
            r_limit_fac,
        )
        cross_field_sum = cross_field_sum + np.sum(cross_field_thread, axis=0)
        cross_volg_sum = cross_volg_sum + np.sum(cross_volg_thread, axis=0)
    return cross_field_sum / cross_volg_sum


@njit
def get_field_volg(dist2, r_max, scalar_now, kernel):
    r_act = np.sqrt(dist2) / r_max
    if kernel == "third":
        r_act3 = r_act * r_act * r_act
        r_fact = np.exp(-r_act3)
    elif kernel == "sharp":
        r_fact = r_act < 1.0  # Alternative: individual ellipsoid visualisation
    return r_fact * scalar_now, r_fact


@njit
def get_abs_dist_periodic(x1, x2, x_extent):
    delx = x1 - x2
    delx = delx - x_extent * round(delx / x_extent)
    return abs(delx)


@njit(parallel=True)
def mpic_cross_section_inner(
    cross_type,
    loc,
    x,
    y,
    z,
    vol,
    scalar,
    xp,
    yp,
    zp,
    extent,
    kernel,
    r_limit_fac,
):
    dx_project = xp[1] - xp[0]
    dy_project = yp[1] - yp[0]
    dz_project = zp[1] - zp[0]
    dxi_project = 1.0 / dx_project
    dyi_project = 1.0 / dy_project
    dzi_project = 1.0 / dz_project
    x_origin = xp[0]
    y_origin = yp[0]
    z_origin = zp[0]
    nxproj = len(xp)
    nyproj = len(yp)
    nzproj = len(zp)
    if cross_type == "xy":
        cross_x_len = nxproj
        cross_y_len = nyproj
    elif cross_type == "xz":
        cross_x_len = nxproj
        cross_y_len = nzproj
    elif cross_type == "yz":
        cross_x_len = nyproj
        cross_y_len = nzproj
    cross_field_thread = np.zeros((NTHREADS, cross_x_len, cross_y_len), np.double)
    cross_volg_thread = np.zeros((NTHREADS, cross_x_len, cross_y_len), np.double)
    ipi = 1.0 / np.pi
    for i_thread in prange(NTHREADS):
        for i in range(len(vol)):
            if not (i % NTHREADS == i_thread):  # Use modulus here
                continue
            x_now = x[i]
            y_now = y[i]
            z_now = z[i]
            vol_now = vol[i]
            r_now = (vol_now * 0.75 * ipi) ** (1.0 / 3.0)
            r_max_safety = r_limit_fac
            r_max = r_max_safety * r_now
            if cross_type == "xy" and (loc < z_now - r_max or loc > z_now + r_max):
                continue
            if (
                cross_type == "xz"
                and get_abs_dist_periodic(loc, y_now, extent[1]) > r_max
            ):
                continue
            if (
                cross_type == "yz"
                and get_abs_dist_periodic(loc, x_now, extent[0]) > r_max
            ):
                continue
            xlower = int(dxi_project * ((x_now - r_max) - x_origin)) + 1
            xupper = int(dxi_project * ((x_now + r_max) - x_origin))
            ylower = int(dyi_project * ((y_now - r_max) - y_origin)) + 1
            yupper = int(dyi_project * ((y_now + r_max) - y_origin))
            zlower = max(int(dzi_project * ((z_now - r_max) - z_origin)), 0)
            zupper = min(int(dzi_project * ((z_now + r_max) - z_origin)), nzproj - 1)
            scalar_now = scalar[i]
            if cross_type == "xy":
                # Include upper values in loops (converted from Fortran)
                for ii in range(xlower, xupper + 1):
                    for jj in range(ylower, yupper + 1):
                        xdist = get_abs_dist_periodic(
                            ii * dx_project + x_origin, x_now, extent[0]
                        )
                        ydist = get_abs_dist_periodic(
                            jj * dy_project + y_origin, y_now, extent[1]
                        )
                        zdist = loc - z_now
                        dist2 = xdist * xdist + ydist * ydist + zdist * zdist
                        if dist2 < r_max * r_max:
                            field_contr, volg_contr = get_field_volg(
                                dist2, r_max, scalar_now, kernel
                            )
                            itarg = ii % nxproj  # Periodic BCs
                            jtarg = jj % nyproj  # Periodic BCs
                            cross_field_thread[i_thread, itarg, jtarg] += field_contr
                            cross_volg_thread[i_thread, itarg, jtarg] += volg_contr
            elif cross_type == "xz":
                # Include upper values in loops (converted from Fortran)
                for ii in range(xlower, xupper + 1):
                    for kk in range(zlower, zupper + 1):
                        xdist = get_abs_dist_periodic(
                            ii * dx_project + x_origin, x_now, extent[0]
                        )
                        ydist = get_abs_dist_periodic(loc, y_now, extent[1])
                        zdist = kk * dz_project + z_origin - z_now
                        dist2 = xdist * xdist + ydist * ydist + zdist * zdist
                        if dist2 < r_max * r_max:
                            field_contr, volg_contr = get_field_volg(
                                dist2, r_max, scalar_now, kernel
                            )
                            itarg = ii % nxproj  # Periodic BCs
                            cross_field_thread[i_thread, itarg, kk] += field_contr
                            cross_volg_thread[i_thread, itarg, kk] += volg_contr
            elif cross_type == "yz":
                # Include upper values in loops (converted from Fortran)
                for jj in range(xlower, xupper + 1):
                    for kk in range(ylower, yupper + 1):
                        xdist = get_abs_dist_periodic(loc, x_now, extent[0])
                        ydist = get_abs_dist_periodic(
                            jj * dy_project + y_origin, y_now, extent[1]
                        )
                        zdist = kk * dz_project + z_origin - z_now
                        dist2 = xdist * xdist + ydist * ydist + zdist * zdist
                        if dist2 < r_max * r_max:
                            field_contr, volg_contr = get_field_volg(
                                dist2, r_max, scalar_now, kernel
                            )
                            jtarg = jj % nyproj  # Periodic BCs
                            cross_field_thread[i_thread, jtarg, kk] += field_contr
                            cross_volg_thread[i_thread, jtarg, kk] += volg_contr
    return cross_field_thread, cross_volg_thread
