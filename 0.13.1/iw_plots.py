import numpy as np
from tools.nc_reader import nc_reader
import matplotlib.pyplot as plt
from tools.mpl_style import *
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib as mpl
from matplotlib.legend_handler import HandlerBase
from tools.utils import *
from tools.mpl_beautify import *
from linear_interpl import get_dataset, time_interpl, find_bounds
from iw_exact import get_exact

data_dir = '../data/iw'


ncr = nc_reader() 

ke_ref  = 5.0e-7
ape_ref = 2.5e-7
en_ref  = 7.5e-7

#
# IW -- timings
#

# EPIC ran with 8 MPI ranks and PS3D ran with 8 OpenMP threads
#==> ../data/iw/epic_iw_48x48x12.csv <==
#function name,#calls,percentage,min time in seconds,mean time in seconds,max time in seconds
#epic,1,100.00,36.824,36.824,36.824
#
#==> ../data/iw/epic_iw_64x64x16.csv <==
#function name,#calls,percentage,min time in seconds,mean time in seconds,max time in seconds
#epic,1,100.00,82.360,82.360,82.361
#
#==> ../data/iw/epic_iw_96x96x24.csv <==
#function name,#calls,percentage,min time in seconds,mean time in seconds,max time in seconds
#epic,1,100.00,272.144,272.144,272.147
#
#==> ../data/iw/epic_iw_128x128x32.csv <==
#function name,#calls,percentage,min time in seconds,mean time in seconds,max time in seconds
#epic,1,100.00,687.792,687.793,687.798
#
#==> ../data/iw/epic_iw_256x256x64.csv <==
#function name,#calls,percentage,min time in seconds,mean time in seconds,max time in seconds
#epic,1,100.00,10499.892,10499.897,10499.934


#==> ../data/iw/ps3d_iw_48x48x12.csv <==
#function name,#calls,total time,percentage
#ps,1,41.767,100.00

#==> ../data/iw/ps3d_iw_64x64x16.csv <==
#function name,#calls,total time,percentage
#ps,1,79.272,100.00

#==> ../data/iw/ps3d_iw_96x96x24.csv <==
#function name,#calls,total time,percentage
#ps,1,199.616,100.00

#==> ../data/iw/ps3d_iw_128x128x32.csv <==
#function name,#calls,total time,percentage
#ps,1,486.119,100.00

#==> ../data/iw/ps3d_iw_256x256x64.csv <==
#function name,#calls,total time,percentage
#ps,1,4539.040,100.00

grid = np.array([12, 16, 24, 32, 64])
xticklabs = [r'$48^2\times 12$', r'$64^2\times 16$',
             r'$96^2\times 24$', r'$128^2\times 32$', r'$256^2\times 64$']
epic_timings = [36.824, 82.361, 272.147, 687.798, 10499.934]
ps3d_timings = [41.767, 79.272, 199.616, 486.119, 4539.040]

plt.figure(figsize=(8, 3), dpi=200)

mpl.rcParams['font.size'] = 12

plt.plot(grid, epic_timings, label=r'EPIC', marker='o')
plt.plot(grid, ps3d_timings, label=r'PS3D', marker='o')
plt.plot(grid, 0.025*grid**3, label=r'$\propto n_z^3$', linestyle='dashed', color='black')
plt.yscale('log', base=10)
plt.xscale('log', base=2)
plt.grid(which='both', linestyle='dashed', zorder=-2)
plt.xticks(grid, xticklabs)
plt.ylabel(r'time (s)')
plt.xlabel(r'grid resolution, $n_x\times n_y\times n_z$')
plt.legend(loc='upper center', ncols=3, bbox_to_anchor=(0.5, 1.2))
plt.savefig('iw_timings.pdf', bbox_inches='tight')
plt.close()

exit()

#
# IW - energy loss
#

te_ref = ape_ref + ke_ref

t_ref = 4.0*np.pi/np.sqrt(2.0)

grids = ['48x48x12', '64x64x16', '96x96x24', '128x128x32', '256x256x64']

for prog in ['epic_iw', 'ps3d_iw']:
    print(prog.upper())
    for i, grid in enumerate(grids):
        if 'epic' in prog:
            fname = os.path.join(data_dir,
                                 prog + '_' + grid + '_parcel_stats.nc')
        else:
            fname = os.path.join(data_dir, prog + '_' + grid + '_field_stats.nc')

        if not os.path.exists(fname):
            print("File does not exist", grid)
            continue

        ncr.open(fname)
        t = ncr.get_all('t')
        lo, hi = find_bounds(t_ref, t)
        if 'ps3d' in prog:
            te = ncr.get_all('ke') + ncr.get_all('ape')
        else:
            te = ncr.get_all('te')

        te_init = te[0]
        te = time_interpl(t_ref, te[lo], t[lo], te[hi], t[hi])
        ncr.close()
        loss = abs(te - te_ref) / te_ref * 1000.0
        print("grid", grid, "theoretical energy loss in TE (permille)", loss, "rounded:",round(loss, 2))

        loss = abs(te - te_init) / te_init * 1000.0
        print("grid", grid, "energy loss in TE (permille)", loss, "rounded:", round(loss, 2))

#
# IW - evolution of raio of energy and enstrophy to their initial values
#

mpl.rcParams['font.size'] = 13
mpl.rcParams['lines.linewidth'] = 1.5

units['time'] = None


##grids = ['48x48x12', '64x64x16', '96x96x24', '128x128x32', '256x256x64']
#grids = ['64x64x16', '128x128x32', '256x256x64']
#annotate = [
##    r'$48^2\times 12$',
#    r'$64^2\times 16$',
##    r'$96^2\times 24$',
#    r'$128^2\times 32$',
#    r'$256^2\times 64$'
#]
#

#n = len(grids)
#
#fig, axs = plt.subplots(nrows=n, ncols=3, figsize=(9, n*2.6), dpi=200, sharex=True)#

#colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#linestyles = ['solid', 'dashed']

#for j in range(3):
#    axs[n-1, j].set_xlabel(r'time, $t$')
#    for i in range(0, n-1):
#        remove_xticks(axs[i, j])
#
#    for i in range(n):
#        axs[i, j].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
#        axs[i, j].grid(linestyle='dashed', zorder=-2)
#        axs[i, j].set_xticks([0, 2.0*np.pi/np.sqrt(2.0), 4.0*np.pi/np.sqrt(2.0)],
#                             [r'$0$', r'$2\pi/\sqrt{2}$', r'$4\pi/\sqrt{2}$'])
#
#labels = [r'EPIC', r'EPIC fields', r'PS3D']
#        
#for i, r in enumerate(grids):
#    
#    ncr.open(os.path.join(data_dir, 'epic_iw_' + r + '_parcel_stats.nc'))
#    t = ncr.get_all('t')
#    ke = ncr.get_all('ke')
#    ape = ncr.get_all('ape')
#    en = ncr.get_all('en')
#    ncr.close()
#
#    label = None
#    if i == 0:
#        label = labels[0]
#
#    axs[i, 0].plot(t, ke - ke_ref, color=colors[0], label=label)
#    axs[i, 1].plot(t, ape - ape_ref, color=colors[0])
#    axs[i, 2].plot(t, en - en_ref, color=colors[0])
#
#    label = None
#    if i == 0:
#        label = labels[1]
#
#    ncr.open(os.path.join(data_dir, 'epic_iw_' + r + '_field_stats.nc'))
#    t = ncr.get_all('t')
#    ke = ncr.get_all('ke')
#    #ape = ncr.get_all('ape')
#    en = ncr.get_all('en')
#    ncr.close()
#
#    axs[i, 0].plot(t, ke - ke_ref, color=colors[1], label=label)
#    #axs[i, 1].plot(t, ape - ape_ref, color=colors[1], label=label)
#    axs[i, 2].plot(t, en - en_ref, color=colors[1])
#    
#    
#    ncr.open(os.path.join(data_dir, 'ps3d_iw_' + r + '_field_stats.nc'))
#    t = ncr.get_all('t')
#    ke = ncr.get_all('ke')
#    ape = ncr.get_all('ape')
#    en = ncr.get_all('en')
#    ncr.close()
#
#    label = None
#    if i == 0:
#        label = labels[2]
#
#    axs[i, 0].plot(t, ke - ke_ref, color=colors[2], linestyle='dashed', label=label)
#    axs[i, 1].plot(t, ape - ape_ref, color=colors[2], linestyle='dashed')
#    axs[i, 2].plot(t, en - en_ref, color=colors[2], linestyle='dashed')
#
#    add_annotation(axs[i, 0], annotate[i], xy=(-0.52, 1.0)) 
#
#for i in range(n):
#    #axs[i, 0].set_ylabel(r'kinetic energy error, $\mathcal{K}(t) - \mathcal{K}^{*}$')
#    axs[i, 0].set_ylabel(r'kin. energy error, $\mathcal{K} - \mathcal{K}^{*}$')
#    #axs[i, 1].set_ylabel(r'avail. pot. energy error, $\mathcal{A}(t) - \mathcal{A}^{*}$')
#    axs[i, 1].set_ylabel(r'APE error, $\mathcal{A} - \mathcal{A}^{*}$')
#    axs[i, 2].set_ylabel(r'enstrophy error, $\Upsilon - \Upsilon^{*}$')
#
##axs[0, 1]
#fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.02), ncols=3)

#fig.tight_layout()

#plt.savefig('iw_ke_ape_en.pdf', bbox_inches='tight')
#plt.close()


#
# IW - time-averaged values of kinetic energy, ape and enstrophy
#


grids = ['48x48x12', '64x64x16', '96x96x24', '128x128x32', '256x256x64']

labels = [r'EPIC', r'PS3D']

kes = np.zeros((len(labels), len(grids)))
apes = np.zeros((len(labels), len(grids)))
ens = np.zeros((len(labels), len(grids)))

for i, r in enumerate(grids):
    ncr.open(os.path.join(data_dir, 'epic_iw_' + r + '_parcel_stats.nc'))
    ke = ncr.get_all('ke')
    ape = ncr.get_all('ape')
    en = ncr.get_all('en')
    ncr.close()

    kes[0, i] = (abs(ke - ke_ref) / ke_ref).mean()
    apes[0, i] = (abs(ape - ape_ref) / ape_ref).mean()
    ens[0, i] = (abs(en - en_ref) / en_ref).mean()
    
    ncr.open(os.path.join(data_dir, 'ps3d_iw_' + r + '_field_stats.nc'))
    ke = ncr.get_all('ke')
    ape = ncr.get_all('ape')
    en = ncr.get_all('en')
    ncr.close()

    kes[1, i] = (abs(ke - ke_ref) / ke_ref).mean()
    apes[1, i] = (abs(ape - ape_ref) / ape_ref).mean()
    ens[1, i] = (abs(en - en_ref) / en_ref).mean()


mpl.rcParams['font.size'] = 13
    
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(9, 3), dpi=200, sharex=True, sharey=False)

z_grids = np.array([12, 16, 24, 32, 64])

for i in range(len(labels)):
    axs[0].plot(z_grids, kes[i, :], color=colors[i], label=labels[i], marker='o')
    axs[1].plot(z_grids, apes[i, :], color=colors[i], marker='o')
    axs[2].plot(z_grids, ens[i, :], color=colors[i], marker='o')

shift = [5.0, 1.5, 1.5]
    
for i in range(3):
    label = None
    if i == 0:
        label = r'$\propto n_z^{-2}$'
    axs[i].plot(z_grids, shift[i]/z_grids**2, label=label, color='black', linestyle='dashed')
    
    axs[i].set_xscale('log', base=2)
    axs[i].set_yscale('log', base=10)
    axs[i].set_xlabel('vertical grid resolution, $n_z$')
    axs[i].grid(which='both', linestyle='dashed', zorder=-2)

axs[0].set_ylabel(
    r'$\langle|\mathcal{K} - \mathcal{K}^{*}|/\mathcal{K}^{*}\rangle_{t}$')
axs[1].set_ylabel(
    r'$\langle|\mathcal{A} - \mathcal{A}^{*}|/\mathcal{A}^{*}\rangle_{t}$')
axs[2].set_ylabel(r'$\langle|\Upsilon - \Upsilon^{*}|/\Upsilon^{*}\rangle_{t}$')

fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncols=4)
fig.tight_layout()
plt.savefig('iw_ke_ape_en_error_conv.pdf', bbox_inches='tight')
plt.close() 

exit()

#
# IW - difference cross section plot of enstrophy at 3 times
#

mpl.rcParams['font.size'] = 12

fig = plt.figure(figsize=(8, 6.5), dpi=200)
grid = ImageGrid(fig, 111,
                 nrows_ncols=(2, 3),
                 aspect=True,
                 axes_pad=(0.3, 0.2),
                 direction='row',
                 share_all=True,
                 cbar_location="right",
                 cbar_mode='single',
                 cbar_size="4%",
                 cbar_pad=0.05)

fnames = ['epic_iw_256x256x64_fields.nc', 'ps3d_iw_256x256x64_fields.nc']
labels = [r'EPIC', r'PS3D']
times = [np.pi/np.sqrt(2.0), 2.0*np.pi/np.sqrt(2.0), 4.0*np.pi/np.sqrt(2.0)]
times_lab = [r'\pi/\sqrt{2}', r'2\pi/\sqrt{2}', r'4\pi/\sqrt{2}']
loc = 128

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

# find global vmin and vmax
for j, fname in enumerate(fnames):
    ncr.open(os.path.join(data_dir, fname))
    nsteps = ncr.get_num_steps()
    origin = ncr.get_box_origin()
    extent = ncr.get_box_extent()
    (nx, ny, nz) = ncr.get_box_ncells()
    t = ncr.get_all('t')
    for step in range(nsteps):
        xi, eta, zeta = get_exact(t[step], origin, extent, nx, ny, nz,
                                  kk=0.5, ll=0.5, mm=1.0, N=2.0, f=1.0,
                                  what=0.001, copy_periodic=True)
        vmag = np.sqrt(xi ** 2 + eta ** 2 + zeta ** 2)
        rms = get_rms(vmag)

        xi = ncr.get_dataset(step=step, name='x_vorticity', copy_periodic=True) - xi
        eta = ncr.get_dataset(step=step, name='y_vorticity', copy_periodic=True) - eta
        zeta = ncr.get_dataset(step=step, name='z_vorticity', copy_periodic=True) - zeta
        
        vmag = np.sqrt(xi ** 2 + eta ** 2 + zeta ** 2) / rms
        vmag = get_plane('xz', loc, vmag)
        vmin = min(vmin, vmag.min())
        vmax = max(vmax, vmag.max())

for j, fname in enumerate(fnames):

    ncr.open(os.path.join(data_dir, fname))
    nsteps = ncr.get_num_steps()
    origin = ncr.get_box_origin()
    extent = ncr.get_box_extent()
    (nx, ny, nz) = ncr.get_box_ncells()
    
    for i, t in enumerate(times):
        ax = grid[i+j*3]

        colorbar = (i == 0)
    
        xi, eta, zeta = get_exact(t, origin, extent, nx, ny, nz,
                                  kk=0.5, ll=0.5, mm=1.0, N=2.0, f=1.0,
                                  what=0.001, copy_periodic=True)

        vmag = np.sqrt(xi ** 2 + eta ** 2 + zeta ** 2)
        rms = get_rms(vmag)

        xi = get_dataset(ncr, 'x_vorticity', t, True) - xi
        eta = get_dataset(ncr, 'y_vorticity', t, True) - eta
        zeta = get_dataset(ncr, 'z_vorticity', t, True) - zeta

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
            add_annotation(ax, r'$t=' + times_lab[i] + r'$', xy=(0.03, 1.25))
    
        if i == 0 or i == 3:
            pass
        else:
            remove_yticks(ax)
            
        #add_annotation(ax, r'$t=' + times_lab[i] + r'$',
        #               xy=(0.03, 1.25)) #, fmt="%.2f")

        add_annotation(grid[j*3], labels[j], xy=(-0.52, 1.0))

    add_annotation(grid[2], r'$y=0$', xy=(0.9, 1.7))

    ncr.close()


ax = fig.add_subplot(3, 1, 3)

ax.grid(which='both', linestyle='dashed')

for i, fname in enumerate(fnames):
    ncr.open(os.path.join(data_dir, fname))
    nsteps = ncr.get_num_steps()
    t = ncr.get_all('t')
    rms = np.zeros(nsteps)
    nsteps = ncr.get_num_steps()
    origin = ncr.get_box_origin()
    extent = ncr.get_box_extent()
    (nx, ny, nz) = ncr.get_box_ncells()

    for step in range(nsteps):
        xi, eta, zeta = get_exact(t[step], origin, extent, nx, ny, nz,
                                  kk=0.5, ll=0.5, mm=1.0, N=2.0, f=1.0,
                                  what=0.001, copy_periodic=True)

        vmag = np.sqrt(xi ** 2 + eta ** 2 + zeta ** 2)
        ex_rms = get_rms(vmag)
        
        xi = ncr.get_dataset(step=step, name='x_vorticity', copy_periodic=True) - xi
        eta = ncr.get_dataset(step=step, name='y_vorticity', copy_periodic=True) - eta
        zeta = ncr.get_dataset(step=step, name='z_vorticity', copy_periodic=True) - zeta

        vmag = np.sqrt(xi ** 2 + eta ** 2 + zeta ** 2)
        rms[step] = get_rms(vmag) / ex_rms

    ncr.close()

    print("Max. normalised r.m.s. error for", labels[i], rms.max())
    print("Mean of normalised r.m.s. error for", labels[i], rms.mean())

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
