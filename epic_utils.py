from numba import njit, prange, set_num_threads
import netCDF4 as nc
import numpy as np
from utils import check_file

MAX_ANISOTROPY_FAC = 5.0
NTHREADS = 8

set_num_threads(NTHREADS)


def get_epic_parcel_coordinates(input_file_name, time_step=0):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    x = ds_nc.variables["x_position"][time_step, :]
    y = ds_nc.variables["y_position"][time_step, :]
    z = ds_nc.variables["z_position"][time_step, :]
    ds_nc.close()
    return x, y, z


def get_epic_parcel_volume(input_file_name, time_step=0):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    vol = ds_nc.variables["volume"][time_step, :]
    ds_nc.close()
    return vol


def get_epic_parcel_shape(input_file_name, time_step=0):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    bb1 = ds_nc.variables["B11"][time_step, :]
    bb2 = ds_nc.variables["B12"][time_step, :]
    bb3 = ds_nc.variables["B13"][time_step, :]
    bb4 = ds_nc.variables["B22"][time_step, :]
    bb5 = ds_nc.variables["B23"][time_step, :]
    ds_nc.close()
    return bb1, bb2, bb3, bb4, bb5


def get_epic_parcel_buoyancy(input_file_name, time_step=0):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    buoyancy = ds_nc.variables["buoyancy"][time_step, :]
    ds_nc.close()
    return buoyancy


def get_epic_parcel_humidity(input_file_name, time_step=0):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    humidity = ds_nc.variables["humidity"][time_step, :]
    ds_nc.close()
    return humidity


def get_epic_parcel_vorticity(input_file_name, time_step=0):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    x_vorticity = ds_nc.variables["x_vorticity"][time_step, :]
    y_vorticity = ds_nc.variables["y_vorticity"][time_step, :]
    z_vorticity = ds_nc.variables["z_vorticity"][time_step, :]
    ds_nc.close()
    return x_vorticity, y_vorticity, z_vorticity


def get_epic_projection_coordinates(input_file_name, refinement_factor):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    # set up spatial arrays
    extent = ds_nc.extent
    origin = ds_nc.origin
    ncells = ds_nc.ncells
    # initialisation
    nx = ncells[0]
    ny = ncells[1]
    nz = ncells[2]
    # projection grid
    nxproj = nx * refinement_factor
    nyproj = ny * refinement_factor
    nzproj = nz * refinement_factor + 1  # Add 1 here
    xp = np.linspace(origin[0], origin[0] + extent[0], nxproj + 1)[0:nxproj]
    yp = np.linspace(origin[1], origin[1] + extent[1], nyproj + 1)[0:nyproj]
    zp = np.linspace(origin[2], origin[2] + extent[2], nzproj)
    ds_nc.close()
    return xp, yp, zp


def get_epic_extent(input_file_name):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    # set up spatial arrays
    extent = ds_nc.extent
    ds_nc.close()
    return extent


def get_idealised_moist_parameters(input_file_name):
    check_file(input_file_name)
    ds_nc = nc.Dataset(input_file_name)
    len_condense = ds_nc.groups["physical_quantities"].scale_height
    q_scale = ds_nc.groups[
        "physical_quantities"
    ].saturation_specific_humidity_at_ground_level
    return len_condense, q_scale


def produce_epic_cross_section(
    input_file_name,
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
    x, y, z = get_epic_parcel_coordinates(input_file_name)
    bb1, bb2, bb3, bb4, bb5 = get_epic_parcel_shape(input_file_name)
    vol = get_epic_parcel_volume(input_file_name)
    # len_condense, q_scale = get_idealised_moist_parameters(input_file_name)
    xp, yp, zp = get_epic_projection_coordinates(input_file_name, refinement_factor)
    h = get_epic_parcel_humidity(input_file_name)
    extent = get_epic_extent(input_file_name)
    cross_field_thread, cross_volg_thread = epic_cross_section_inner(
        cross_type,
        loc,
        x,
        y,
        z,
        bb1,
        bb2,
        bb3,
        bb4,
        bb5,
        vol,
        h,
        xp,
        yp,
        zp,
        extent,
        kernel,
        r_limit_fac,
    )
    return np.sum(cross_field_thread, axis=0) / np.sum(cross_volg_thread, axis=0)


@njit
def get_field_volg(xdist, ydist, zdist, bb, scalar_now, r_limit_fac, kernel):
    bbinv = np.zeros((3, 3), np.double)
    bbinv[:, :] = np.linalg.inv(bb)
    xarr = np.zeros((3, 1), np.double)
    xarr[0, 0] = xdist
    xarr[1, 0] = ydist
    xarr[2, 0] = zdist
    r_act = np.sqrt(
        np.dot(
            np.dot(xarr.transpose(), bbinv),
            xarr,
        )[0, 0]
    )
    if r_act < r_limit_fac:
        if kernel == "third":
            r_act3 = r_act * r_act * r_act
            r_fact = np.exp(-r_act3)
        elif kernel == "sharp":
            r_fact = r_act < 1.0  # Alternative: individual ellipsoid visualisation
        contr_flag = True
    else:
        r_fact = 0.0
        contr_flag = False
    return r_fact * scalar_now, r_fact, contr_flag


@njit
def get_abs_dist_periodic(x1, x2, x_extent):
    delx = x1 - x2
    delx = delx - x_extent * round(delx / x_extent)
    return abs(delx)


@njit(parallel=True)
def epic_cross_section_inner(
    cross_type,
    loc,
    x,
    y,
    z,
    bb1,
    bb2,
    bb3,
    bb4,
    bb5,
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
        bb = np.zeros((3, 3), np.double)
        evals = np.zeros(3, np.double)
        for i in range(len(vol)):
            if not (i % NTHREADS == i_thread):  # Use modulus here
                continue
            x_now = x[i]
            y_now = y[i]
            z_now = z[i]
            vol_now = vol[i]
            r_now = (vol_now * 0.75 * ipi) ** (1.0 / 3.0)
            r_max_safety = r_limit_fac * MAX_ANISOTROPY_FAC
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
            # Now let's do this property
            b1 = bb1[i]
            b2 = bb2[i]
            b3 = bb3[i]
            b4 = bb4[i]
            b5 = bb5[i]
            abc = 0.75 * vol_now * ipi
            b6 = (
                abc * abc - b3 * (b2 * b5 - b3 * b4) + b1 * b5 * b5 - b2 * b3 * b5
            ) / (b1 * b4 - b2 * b2)
            bb[0, 0] = b1
            bb[0, 1] = b2
            bb[0, 2] = b3
            bb[1, 0] = b2
            bb[1, 1] = b4
            bb[1, 2] = b5
            bb[2, 0] = b3
            bb[2, 1] = b5
            bb[2, 2] = b6
            evals[:] = np.linalg.eigvalsh(bb)
            anisotropy_fact = np.sqrt(np.nanmax(evals)) / (
                (evals[0] * evals[1] * evals[2]) ** (1.0 / 6.0)
            )
            r_max_safety = r_limit_fac * anisotropy_fact
            r_max = r_max_safety * r_now
            # Determine which part of the grid the parcel contributes to
            if cross_type == "xy" and (loc < z_now - r_max or loc > z_now + r_max):
                continue
            if cross_type == "xz" and (loc < y_now - r_max or loc > y_now + r_max):
                continue
            if cross_type == "yz" and (loc < x_now - r_max or loc > x_now + r_max):
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
                            field_contr, volg_contr, contr_flag = get_field_volg(
                                xdist, ydist, zdist, bb, scalar_now, r_limit_fac, kernel
                            )
                            if contr_flag:
                                itarg = ii % nxproj  # Periodic BCs
                                jtarg = jj % nyproj  # Periodic BCs
                                cross_field_thread[
                                    i_thread, itarg, jtarg
                                ] += field_contr
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
                            field_contr, volg_contr, contr_flag = get_field_volg(
                                xdist, ydist, zdist, bb, scalar_now, r_limit_fac, kernel
                            )
                            if contr_flag:
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
                            field_contr, volg_contr, contr_flag = get_field_volg(
                                xdist, ydist, zdist, bb, scalar_now, r_limit_fac, kernel
                            )
                            if contr_flag:
                                jtarg = jj % nyproj  # Periodic BCs
                                cross_field_thread[i_thread, jtarg, kk] += field_contr
                                cross_volg_thread[i_thread, jtarg, kk] += volg_contr
    return cross_field_thread, cross_volg_thread
