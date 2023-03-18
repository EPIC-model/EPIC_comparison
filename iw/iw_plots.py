import numpy as np
from tools.nc_reader import nc_reader
import matplotlib.pyplot as plt
from tools.mpl_style import *
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib as mpl
from matplotlib.legend_handler import HandlerBase
from tools.utils import *
from tools.mpl_beautify import *

#
# IW - evolution of raio of energy and enstrophy to their initial values
#

mpl.rcParams['font.size'] = 13
mpl.rcParams['lines.linewidth'] = 1.5

units['time'] = None

ncr = nc_reader()

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(9, 2*2.8), dpi=200, sharex=True)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
linestyles = ['solid', 'dashed']

for j in range(3):
    axs[1, j].set_xlabel(r'time, $t$')
    remove_xticks(axs[0, j])

    for i in range(2):
        axs[i, j].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        axs[i, j].grid(linestyle='dashed', zorder=-2)
        axs[i, j].set_xticks([0, 2.0*np.pi/np.sqrt(2.0), 4.0*np.pi/np.sqrt(2.0)],
                             [r'$0$', r'$2\pi/\sqrt{2}$', r'$4\pi/\sqrt{2}$'])

labels = [r'EPIC', r'PS3D']
annotate = [r'$64^2\times 16$', r'$128^2\times 32$']
        
for i, r in enumerate(['64x64x16', '128x128x32']):
    
    ncr.open('epic_iw_' + r + '_parcel_stats.nc')
    t = ncr.get_all('t')
    ke = ncr.get_all('ke')
    ape = ncr.get_all('ape')
    en = ncr.get_all('en')
    ncr.close()

    label = None
    if i == 0:
        label = labels[0]

    axs[i, 0].plot(t, ke, color=colors[0])
    axs[i, 1].plot(t, ape, color=colors[0], label=label)
    axs[i, 2].plot(t, en, color=colors[0])
    
    ncr.open('ps3d_iw_' + r + '_field_stats.nc')
    t = ncr.get_all('t')
    ke = ncr.get_all('ke')
    ape = ncr.get_all('ape')
    en = ncr.get_all('en')
    ncr.close()

    label = None
    if i == 0:
        label = labels[1]

    axs[i, 0].plot(t, ke, color=colors[1], linestyle='dashed')
    axs[i, 1].plot(t, ape, color=colors[1], linestyle='dashed', label=label)
    axs[i, 2].plot(t, en, color=colors[1], linestyle='dashed')

    add_annotation(axs[i, 0], annotate[i], xy=(-0.52, 1.1)) 

for i in range(2):
    axs[i, 0].set_ylabel(r'kinetic energy, $\mathcal{K}(t)$')
    axs[i, 1].set_ylabel(r'avail. pot. energy, $\mathcal{A}(t)$')
    axs[i, 2].set_ylabel(r'enstrophy, $\Upsilon(t)$')

axs[0, 1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.4), ncols=2)

plt.tight_layout()

plt.subplots_adjust(wspace=0.5)

plt.savefig('iw_ke_ape_en.pdf', bbox_inches='tight')
plt.close()

#
# IW - difference cross section plot of enstrophy at 3 times
#

# get exact field
def get_exact(t):
    N = 2.0
    f = 1.0
    kk = 0.5
    ll = 0.5
    mm = 1.0
    what = 0.001
    nx = 128
    ny = 128
    nz = 32
    
    N2 = N ** 2
    f2 = f ** 2
    sigma2 = (N2 * (kk ** 2 + ll ** 2) + f2 * mm ** 2) / (kk ** 2 + ll ** 2 + mm ** 2)

    # domain origin
    origin = (-2.0 * np.pi, -2.0 * np.pi, -0.5 * np.pi)

    # domain extent
    extent = (4.0 * np.pi, 4.0 * np.pi, np.pi)


    # mesh spacings
    dx = extent[0] / nx
    dy = extent[1] / ny
    dz = extent[2] / nz

    xi = np.zeros((nx+1, ny+1, nz+1))
    eta = np.zeros((nx+1, ny+1, nz+1))
    zeta = np.zeros((nx+1, ny+1, nz+1))
    buoy = np.zeros((nx+1, ny+1, nz+1))

    N2f2 = N2 - f2
    s2f2i = 1.0 / (sigma2 - f2)
    sigma = np.sqrt(sigma2)
    sigi = 1.0 / sigma
    N2s2si = (N2 - sigma2) * sigi

    # ranges from 0 to nx-1
    for i in range(nx+1):
        for j in range(ny+1):
            for k in range(nz+1):
                x = origin[0] + i * dx
                y = origin[1] + j * dy
                z = origin[2] + k * dz

                phi = kk * x + ll * y - sigma * t

                sinphi = np.sin(phi)
                cosphi = np.cos(phi)
                cosmz = np.cos(mm * z)

                xi[i, j, k] = s2f2i * what * cosmz * (f * kk * N2s2si * cosphi - ll * N2f2 * sinphi)

                eta[i, j, k] = s2f2i * what * cosmz * (f * ll * N2s2si * cosphi + kk * N2f2 * sinphi)

                zeta[i, j, k] = f * mm * what * sigi * np.sin(mm * z) * sinphi

    return xi, eta, zeta


mpl.rcParams['font.size'] = 12

fig = plt.figure(figsize=(8, 6.5), dpi=200)
grid = ImageGrid(fig, 111,
                 nrows_ncols=(2, 3),
                 aspect=True,
                 axes_pad=(0.3, 0.33),
                 direction='row',
                 share_all=True,
                 cbar_location="right",
                 cbar_mode='single',
                 cbar_size="4%",
                 cbar_pad=0.05)

fnames = ['epic_iw_128x128x32_fields.nc', 'ps3d_iw_128x128x32_fields.nc']
labels = [r'EPIC', r'PS3D']
steps = [9, 18, 36]
times = [r'\pi\sqrt{2}', r'2\pi\sqrt{2}', r'4\pi/\sqrt{2}']
loc = 64

yticks   = np.pi * np.array([-0.5, 0.0, 0.5])                                            
yticklab = [r'$-\pi/2$', r'$0$', r'$\pi/2$']
    
xticks   = np.pi * np.array([-2, -1, 0.0, 1, 2])
xticklab = [r'$-2\pi$', r'$-\pi$', r'$0$', r'$\pi$', r'$2\pi$']

vmin = 10000
vmax = -10000

def get_rms(data):
    nx, ny, nz = data.shape

    # exclude periodic layers when computing the rms
    nx = nx - 1
    ny = ny - 1

    rms = (data[0:nx, 0:ny, 1:nz-1] ** 2).sum()
    rms = rms + 0.5 * (data[0:nx, 0:ny, 0] ** 2).sum()
    rms = rms + 0.5 * (data[0:nx, 0:ny, nz-1] ** 2).sum()
    return np.sqrt(rms / (nx * ny * (nz-1)))

for j, fname in enumerate(fnames):
    ncr.open(fname)
    nsteps = ncr.get_num_steps()
    t = ncr.get_all('t')
    for step in range(nsteps):
        xi, eta, zeta = get_exact(t[step])
        vmag = np.sqrt(xi ** 2 + eta ** 2 + zeta ** 2)
        rms = get_rms(vmag)
        
        xi = ncr.get_dataset(step=step, name='x_vorticity') - xi
        eta = ncr.get_dataset(step=step, name='y_vorticity') - eta
        zeta = ncr.get_dataset(step=step, name='z_vorticity') - zeta

        vmag = np.sqrt(xi ** 2 + eta ** 2 + zeta ** 2) / rms
        vmag = get_plane('xz', loc, vmag)
        vmin = min(vmin, vmag.min())
        vmax = max(vmax, vmag.max())

        
for j, fname in enumerate(fnames):

    ncr.open(fname)

    t = ncr.get_all('t')

    for i, step in enumerate(steps):
        ax = grid[i+j*3]

        colorbar = (i == 0)
    
        xi, eta, zeta = get_exact(t[step])

        vmag = np.sqrt(xi ** 2 + eta ** 2 + zeta ** 2)
        rms = get_rms(vmag)

        xi = ncr.get_dataset(step=step, name='x_vorticity') - xi
        eta = ncr.get_dataset(step=step, name='y_vorticity') - eta
        zeta = ncr.get_dataset(step=step, name='z_vorticity') - zeta

        vmag = np.sqrt(xi ** 2 + eta ** 2 + zeta ** 2) / rms

        im, cbar = make_imshow(
            ax=ax,
            plane='xz',
            loc=loc,
            fdata=vmag,
            vmin=vmin,
            vmax=vmax,
            ncr=ncr,
            cmap='Reds',
            cmap_norm=None,
            clabel=r'$|\bm{\omega} - \bm{\omega}^{\star}|/\omega^{\star}_{\mathrm{rms}}$',
            colorbar=colorbar,
            xticks=xticks,
            xticklab=xticklab,
            yticks=yticks,
            yticklab=yticklab)

        xg, yg, zg = ncr.get_meshgrid()
        vmag = get_plane('xz', loc, vmag)
        ax.contour(xg[:, loc, :], zg[:, loc, :], vmag, 10, colors='black', linewidths=0.2)

        if j == 0:
            remove_xticks(ax)
    
        if i == 0 or i == 3:
            pass
        else:
            remove_yticks(ax)
            
        add_annotation(ax, r'$t=' + str(round(t[step], 2)) + r'\approx' + times[i] + r'$',
                       xy=(0.03, 1.25)) #, fmt="%.2f")

        add_annotation(grid[j*3], labels[j], xy=(-0.52, 1.25))

    add_annotation(grid[2], r'$y=0$', xy=(0.9, 1.7))

    ncr.close()


ax = fig.add_subplot(3, 1, 3)

ax.grid(which='both', linestyle='dashed')

for i, fname in enumerate(fnames):
    ncr.open(fname)
    nsteps = ncr.get_num_steps()
    t = ncr.get_all('t')
    rms = np.zeros(nsteps)

    for step in range(nsteps):
        xi, eta, zeta = get_exact(t[step])

        vmag = np.sqrt(xi ** 2 + eta ** 2 + zeta ** 2)
        ex_rms = get_rms(vmag)
        
        xi = ncr.get_dataset(step=step, name='x_vorticity') - xi
        eta = ncr.get_dataset(step=step, name='y_vorticity') - eta
        zeta = ncr.get_dataset(step=step, name='z_vorticity') - zeta

        rms[step] = get_rms(np.sqrt(xi ** 2 + eta ** 2 + zeta ** 2)) / ex_rms

    ncr.close()

    print("Max. normalised r.m.s. error for", labels[i], rms.max())

    ax.plot(t, rms, label=labels[i])

ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
pisqrt2 = np.pi/np.sqrt(2.0)
ax.set_xticks([0.0, pisqrt2, 2.0*pisqrt2, 3.0*pisqrt2, 4.0*pisqrt2],
              [r'$0$', r'$\pi/\sqrt{2}$', r'$2\pi/\sqrt{2}$', r'$3\pi/\sqrt{2}$', r'$4\pi/\sqrt{2}$'])
ax.set_xlabel(r'time, $t$')
ax.set_ylabel(r'$|\bm{\omega}-\bm{\omega}^{\star}|_{\mathrm{rms}} / \omega^{\star}_{\mathrm{rms}}$')
ax.legend(loc='best')

plt.subplots_adjust(hspace=0.4)
 
plt.savefig('iw_cross_sections.pdf', bbox_inches='tight')
plt.close()
