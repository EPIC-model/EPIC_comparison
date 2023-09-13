import numpy as np

def get_exact(t, origin, extent, nx, ny, nz, kk, ll, mm, N, f, what, copy_periodic):
    N2 = N ** 2
    f2 = f ** 2
    sigma2 = (N2 * (kk ** 2 + ll ** 2) + f2 * mm ** 2) / (kk ** 2 + ll ** 2 + mm ** 2)
    sigma = np.sqrt(sigma2)

    # mesh spacings
    dx = extent[0] / nx
    dy = extent[1] / ny
    dz = extent[2] / nz

    nnx = nx
    nny = ny
    if copy_periodic:
        nnx = nx + 1
        nny = ny + 1
        
    xi = np.zeros((nnx, nny, nz+1))
    eta = np.zeros((nnx, nny, nz+1))
    zeta = np.zeros((nnx, nny, nz+1))
    u = np.zeros((nnx, nny, nz+1))
    v = np.zeros((nnx, nny, nz+1))
    w = np.zeros((nnx, nny, nz+1))
    buoy = np.zeros((nnx, nny, nz+1))

    N2f2 = N2 - f2
    s2f2i = 1.0 / (sigma2 - f2)
    sigma = np.sqrt(sigma2)
    sigi = 1.0 / sigma
    N2s2si = (N2 - sigma2) * sigi
    k2l2i = 1.0 / (kk ** 2 + ll ** 2)

    for i in range(nx):
        for j in range(ny):
            for k in range(nz+1):
                x = origin[0] + i * dx
                y = origin[1] + j * dy
                z = origin[2] + k * dz
                
                phi = kk * x + ll * y - sigma * t
                
                sinphi = np.sin(phi)
                cosphi = np.cos(phi)
                cosmz = np.cos(mm * z)
                sinmz = np.sin(mm * z)
                
                xi[i, j, k] = s2f2i * what * cosmz * (f * kk * N2s2si * cosphi - ll * N2f2 * sinphi)
                
                eta[i, j, k] = s2f2i * what * cosmz * (f * ll * N2s2si * cosphi + kk * N2f2 * sinphi)
                
                zeta[i, j, k] = f * mm * what * sigi * np.sin(mm * z) * sinphi

                u[i, j, k] = k2l2i * what * mm * sinmz * (kk * sinphi + sigi * f * ll * cosphi)
                v[i, j, k] = k2l2i * what * mm * sinmz * (ll * sinphi - sigi * f * kk * cosphi)
                w[i, j, k] = what * cosmz * cosphi

                buoy[i, j, k] = N2 * z + sigi * N2 * what * cosmz * sinphi

    if copy_periodic:
        xi[nx, :, :] = xi[0, :, :]
        xi[:, ny, :] = xi[:, 0, :]
        
        eta[nx, :, :] = eta[0, :, :]
        eta[:, ny, :] = eta[:, 0, :]

        zeta[nx, :, :] = zeta[0, :, :]
        zeta[:, ny, :] = zeta[:, 0, :]

        u[nx, :, :] = u[0, :, :]
        u[:, ny, :] = u[:, 0, :]
        
        v[nx, :, :] = v[0, :, :]
        v[:, ny, :] = v[:, 0, :]
        
        w[nx, :, :] = w[0, :, :]
        w[:, ny, :] = w[:, 0, :]

        buoy[nx, :, :] = buoy[0, :, :]
        buoy[:, ny, :] = buoy[:, 0, :]
        
    return xi, eta, zeta, u, v, w, buoy
