#!/usr/bin/env python
#
# 3D Beltrami flow:
#   velocity:
#     u(x, y, z) = 1/4 * [sin(z) - 3 * cos(z)] * sin(2x + 2y)
#     v(x, y, z) = 1/4 * [sin(z) + 3 * cos(z)] * sin(2x + 2y)
#     w(x, y, z) = cos(z) * cos(2x + 2y)
#
#   vorticitiy:
#     \xi(x, y, z)   = 3 * u(x, y, z) + \xi_pert
#     \eta(x, y, z)  = 3 * v(x, y, z) + \eta_pert
#     \zeta(x, y, z) = 3 * w(x, y, z) + \zeta_pert
#
# where
#    \xi_pert  = a * cos(2y) * cos(z)
#    \eta_pert = b * cos(2x) * cos(z)
#    \zeta_pert = 0.0
#
# and a = 1/5 and b = 1/10.
#
from tools.nc_fields import nc_fields
import numpy as np
import argparse

try:
    parser = argparse.ArgumentParser(
        description="Create Beltrami flow"
    )

    parser.add_argument(
        "--nx",
        type=int,
        required=False,
        default=32,
        help="number of cells in x"
    )

    parser.add_argument(
        "--ny",
        type=int,
        required=False,
        default=32,
        help="number of cells in y"
    )

    parser.add_argument(
        "--nz",
        type=int,
        required=False,
        default=32,
        help="number of cells in z"
    )

    parser.add_argument(
        "--a",
        type=float,
        required=False,
        default=0.2,
        help="perturbation strength in x"
    )

    parser.add_argument(
        "--b",
        type=float,
        required=False,
        default=0.1,
        help="perturbation strength in y"
    )

    args = parser.parse_args()

    # number of cells
    nx = args.nx
    ny = args.ny
    nz = args.nz

    a = args.a
    b = args.b

    ncf = nc_fields()

    ncf.open('beltrami_' + str(nx) + 'x' + str(ny) + 'x' + str(nz) + '.nc')

    # domain origin
    origin = (-0.5 * np.pi, -0.5 * np.pi, -0.5 * np.pi)

    # domain extent
    extent = (np.pi, np.pi, np.pi)


    # mesh spacings
    dx = extent[0] / nx
    dy = extent[1] / ny
    dz = extent[2] / nz

    xi   = np.zeros((nz+1, ny, nx))
    eta  = np.zeros((nz+1, ny, nx))
    zeta = np.zeros((nz+1, ny, nx))

    # ranges from 0 to nx-1
    for i in range(nx):
        for j in range(ny):
            # ranges from 0 to nz
            for k in range(nz+1):
                x = origin[0] + i * dx
                y = origin[1] + j * dy
                z = origin[2] + k * dz

                # velocity:
                u = 0.25 * (np.sin(z) - 3.0 * np.cos(z)) * np.sin(2.0 * x + 2.0 * y)
                v = 0.25 * (np.sin(z) + 3.0 * np.cos(z)) * np.sin(2.0 * x + 2.0 * y)
                w = np.cos(z) * np.cos(2.0 * x + 2.0 * y)

                # perturbation:
                xi_pert = a * np.cos(2.0 * y) * np.cos(z)
                eta_pert = b * np.cos(2.0 * x) * np.cos(z)
                zeta_pert = 0.0

                # vorticity:
                xi[k, j, i]   = 3.0 * u + xi_pert
                eta[k, j, i]  = 3.0 * v + eta_pert
                zeta[k, j, i] = 3.0 * w + zeta_pert

    # write all provided fields
    ncf.add_field('x_vorticity', xi, unit='1/s')
    ncf.add_field('y_vorticity', eta, unit='1/s')
    ncf.add_field('z_vorticity', zeta, unit='1/s')

    ncf.add_axis('t', [0.0])

    ncf.add_box(origin, extent, [nx, ny, nz])

    ncf.close()

except Exception as err:
    print(err)
