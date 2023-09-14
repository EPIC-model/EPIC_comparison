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

plot_iw_zeta_cross = True


ncr = nc_reader() 

ke_ref  = 5.0e-7
ape_ref = 2.5e-7
en_ref  = 7.5e-7

#
# IW -- timings
#

# EPIC ran with 8 MPI ranks and PS3D ran with 8 OpenMP threads
#head -n 2 ../data/iw/epic_iw_*.csv
#==> ../data/iw/epic_iw_48x48x12.csv <==
#function name,#calls,percentage,min time in seconds,mean time in seconds,max time in seconds
#epic,1,100.00,37.221,37.221,37.221
#
#==> ../data/iw/epic_iw_64x64x16.csv <==
#function name,#calls,percentage,min time in seconds,mean time in seconds,max time in seconds
#epic,1,100.00,84.474,84.474,84.475
#
#==> ../data/iw/epic_iw_96x96x24.csv <==
#function name,#calls,percentage,min time in seconds,mean time in seconds,max time in seconds
#epic,1,100.00,272.543,272.543,272.545
#
#==> ../data/iw/epic_iw_128x128x32.csv <==
#function name,#calls,percentage,min time in seconds,mean time in seconds,max time in seconds
#epic,1,100.00,689.483,689.483,689.489
#
#==> ../data/iw/epic_iw_256x256x64.csv <==
#function name,#calls,percentage,min time in seconds,mean time in seconds,max time in seconds
#epic,1,100.00,10673.100,10673.104,10673.129

#head -n 2 ../data/iw/ps3d_iw_*.csv
#==> ../data/iw/ps3d_iw_48x48x12.csv <==
#function name,#calls,total time,percentage
#ps,1,404.282,100.00
#
#==> ../data/iw/ps3d_iw_64x64x16.csv <==
#function name,#calls,total time,percentage
#ps,1,779.599,100.00
#
#==> ../data/iw/ps3d_iw_96x96x24.csv <==
#function name,#calls,total time,percentage
#ps,1,1949.638,100.00
#
#==> ../data/iw/ps3d_iw_128x128x32.csv <==
#function name,#calls,total time,percentage
#ps,1,4744.732,100.00
#
#==> ../data/iw/ps3d_iw_256x256x64.csv <==
#function name,#calls,total time,percentage
#ps,1,42869.971,100.00


grid = np.array([12, 16, 24, 32, 64])
xticklabs = [r'$48^2\times 12$', r'$64^2\times 16$',
             r'$96^2\times 24$', r'$128^2\times 32$', r'$256^2\times 64$']
epic_timings = [37.221, 84.475, 272.545, 689.489, 10673.129]
ps3d_timings = [404.282, 779.599, 1949.638, 4744.732, 42869.971]

plt.figure(figsize=(8, 2), dpi=200)

mpl.rcParams['font.size'] = 11

plt.plot(grid, epic_timings, label=r'EPIC', marker='o')
plt.plot(grid, ps3d_timings, label=r'PS3D', marker='o')
plt.plot(grid, 0.05*grid**3, label=r'$\propto n_z^3$', linestyle='dashed', color='black')
plt.yscale('log', base=10)
plt.xscale('log', base=2)
plt.grid(which='both', linestyle='dashed', zorder=-2)
plt.xticks(grid, xticklabs)
plt.ylabel(r'time (s)')
plt.xlabel(r'grid resolution, $n_x\times n_y\times n_z$')
plt.legend(loc='upper center', ncols=3, bbox_to_anchor=(0.5, 1.3))
plt.savefig('iw_timings.pdf', bbox_inches='tight')
plt.close()

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
        print("grid", grid, "theoretical energy loss in TE (permille)", loss, "rounded:",round(loss, 3))

        loss = abs(te - te_init) / te_init * 1000.0
        print("grid", grid, "energy loss in TE (permille)", loss, "rounded:", round(loss, 3))

#
# IW -- energy, enstrophy convergence plot
#

mpl.rcParams['font.size'] = 14
mpl.rcParams['lines.linewidth'] = 1.5

units['time'] = None

grids = ['48x48x12', '64x64x16', '96x96x24', '128x128x32', '256x256x64']
glabs = [r'$16$', r'$32$', r'$64$']
xticks = [16, 32, 64]
z_grids = np.array([12, 16, 24, 32, 64])

labels = [r'EPIC', r'PS3D']

kes = np.zeros((len(labels), len(grids)))
apes = np.zeros((len(labels), len(grids)))
tes = np.zeros((len(labels), len(grids)))
ens = np.zeros((len(labels), len(grids)))

t_ref = 4.0*np.pi/np.sqrt(2.0)

for i, r in enumerate(grids):
    for j, fname in enumerate(['epic_iw_' + r + '_parcel_stats.nc',
                               'ps3d_iw_' + r + '_field_stats.nc']):
        
        ncr.open(os.path.join(data_dir, fname))
        ke = ncr.get_all('ke')
        ape = ncr.get_all('ape')
        en = ncr.get_all('en')

        if 'ps3d' in fname:
            te = ke + ape
        else:
            te = ncr.get_all('te')
        
        t = ncr.get_all('t')
        ncr.close()

        kes[j, i]  = abs(1.0 - ke  / ke.mean()).mean()
        apes[j, i] = abs(1.0 - ape / ape.mean()).mean()
        tes[j, i]  = abs(1.0 - te  / te.mean()).mean()
        ens[j, i]  = abs(1.0 - en  / en.mean()).mean()
    
fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(9.5, 2.25), dpi=200, sharex=True, sharey=False)

axs = axs.flatten()

for i in range(len(labels)):
    axs[0].plot(z_grids, kes[i, :], color=colors[i], label=labels[i], marker='o')
    axs[1].plot(z_grids, apes[i, :], color=colors[i], marker='o')
    axs[2].plot(z_grids, tes[i, :], color=colors[i], marker='o')
    axs[3].plot(z_grids, ens[i, :], color=colors[i], marker='o')

shift = [0.15, 0.4, 0.2, 0.15]
power = [2, 2, 3, 2]
labels = [r'$\propto n_z^{-2}$', None, r'$\propto n_z^{-3}$', None]
linestyles = {
    2: 'dashed',
    3: 'dashdot'
}

for i in range(4):
    axs[i].plot(z_grids, shift[i]/z_grids**power[i], label=labels[i],
                color='black', linestyle=linestyles[power[i]])
    
    axs[i].set_xscale('log', base=2)
    axs[i].set_yscale('log', base=10)
    axs[i].set_xticks(xticks, glabs)
    axs[i].set_xlabel(r'grid resolution, $n_z$')
    axs[i].grid(which='both', linestyle='dashed', zorder=-2)

axs[0].set_ylabel(r'$\langle\lvert 1 - \mathcal{K}/\langle\mathcal{K}\rangle\rvert\rangle$')
axs[1].set_ylabel(r'$\langle\lvert 1 -  \mathcal{A}/\langle\mathcal{A}\rangle\rvert\rangle$')
axs[2].set_ylabel(r'$\langle\lvert 1 - \mathcal{T}/\langle\mathcal{T}\rangle\rvert\rangle$')
axs[3].set_ylabel(r'$\langle\lvert 1 - \Upsilon/\langle\Upsilon\rangle\rvert\rangle$')

fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.14), ncols=4)
plt.tight_layout()
plt.savefig('iw_ke_ape_te_en_error_conv.pdf', bbox_inches='tight')
plt.close()


#
# IW - Relative error to theoretical value
#
mpl.rcParams['font.size'] = 13
mpl.rcParams['lines.linewidth'] = 1.5

grids = ['48x48x12', '64x64x16', '96x96x24', '128x128x32', '256x256x64']
annotate = [
    r'$48^2\times 12$',
    r'$64^2\times 16$',
    r'$96^2\times 24$',
    r'$128^2\times 32$',
    r'$256^2\times 64$'
]


n = len(grids)

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(9, 5), dpi=200, sharex=True)#

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
linestyles = ['solid', 'dashed']

for i in range(3):
    axs[1, i].set_xlabel(r'time, $t$')
    for j in range(2):
        axs[j, i].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        axs[j, i].grid(linestyle='dashed', zorder=-2)
    remove_xticks(axs[0, i])
    axs[1, i].set_xticks([0, 2.0*np.pi/np.sqrt(2.0), 4.0*np.pi/np.sqrt(2.0)],
                         [r'$0$', r'$2\pi/\sqrt{2}$', r'$4\pi/\sqrt{2}$'])

labels = [r'EPIC', r'PS3D']

for i, r in enumerate(grids):
    for j, fname in enumerate(['epic_iw_' + r + '_parcel_stats.nc', 'ps3d_iw_' + r + '_field_stats.nc']):
        ncr.open(os.path.join(data_dir, fname))
        t = ncr.get_all('t')
        ke = ncr.get_all('ke')
        ape = ncr.get_all('ape')
        en = ncr.get_all('en')
        ncr.close()

        lab = None
        if j == 0:
            lab = annotate[i]
        
        axs[j, 0].plot(t, 1.0 - ke / ke_ref, color=colors[i], label=lab)
        axs[j, 1].plot(t, 1.0 - ape / ape_ref, color=colors[i])
        axs[j, 2].plot(t, 1.0 - en / en_ref, color=colors[i])

        #print(fname, r, "max KE", (1.0 - ke / ke_ref).max())
        #print(fname, r, "max AP", (1.0 - ape / ape_ref).max())
        #print(fname, r, "max EN", (1.0 - en / en_ref).max())

add_annotation(axs[0, 0], labels[0], xy=(-0.52, 1.0))
add_annotation(axs[1, 0], labels[1], xy=(-0.52, 1.0))

for i in range(len(labels)):
    axs[i, 0].set_ylabel(r'$1 - \mathcal{K} / \mathcal{K}^{\star}$')
    axs[i, 1].set_ylabel(r'$1 - \mathcal{A} / \mathcal{A}^{\star}$')
    axs[i, 2].set_ylabel(r'$1 - \Upsilon /  \Upsilon^{\star}$')

fig.legend(loc='upper center', bbox_to_anchor=(0.55, 1.05), ncols=5)
fig.tight_layout()
plt.savefig('iw_ke_ape_en_diff_theo_val.pdf', bbox_inches='tight')
plt.close()


#
# IW - difference cross section plot of zeta at 3 times
#

if plot_iw_zeta_cross:
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


    ax = fig.add_subplot(3, 1, 3)
    ax.grid(which='both', linestyle='dashed')

def get_rms(data):
    nx, ny, nz = data.shape

    # exclude periodic layers when computing the rms
    nx = nx - 1
    ny = ny - 1

    arr = data[0:nx, 0:ny, :]

    #rms = (data[0:nx, 0:ny, 1:nz-1] ** 2).sum()
    #rms = rms + 0.5 * (data[0:nx, 0:ny, 0] ** 2).sum()
    #rms = rms + 0.5 * (data[0:nx, 0:ny, nz-1] ** 2).sum()
    #return np.sqrt(rms / (nx * ny * (nz-1)))

    return np.sqrt((arr ** 2).mean())


if plot_iw_zeta_cross:
    # find global vmin and vmax, and fill rms arrays
    for j, fname in enumerate(fnames):
        ncr.open(os.path.join(data_dir, fname))
        nsteps = ncr.get_num_steps()
        origin = ncr.get_box_origin()
        extent = ncr.get_box_extent()
        (nx, ny, nz) = ncr.get_box_ncells()
        t = ncr.get_all('t')

        rmse = np.zeros(nsteps)

        for tt in times:
            _, _, zeta, _, _, _, _ = get_exact(tt, origin, extent, nx, ny, nz,
                                               kk=0.5, ll=0.5, mm=1.0, N=2.0, f=1.0,
                                               what=0.001, copy_periodic=True)
            rms = get_rms(zeta)

            zeta = get_dataset(ncr, 'z_vorticity', tt, True) - zeta

            zeta = zeta / rms
            zeta = get_plane('xz', loc, zeta)
            vmin = min(vmin, zeta.min())
            vmax = max(vmax, zeta.max())


        for step in range(nsteps):
            _, _, zeta, _, _, _, _ = get_exact(t[step], origin, extent, nx, ny, nz,
                                               kk=0.5, ll=0.5, mm=1.0, N=2.0, f=1.0,
                                               what=0.001, copy_periodic=True)

            rms = get_rms(zeta)

            zeta = ncr.get_dataset(step=step, name='z_vorticity', copy_periodic=True) - zeta
        
            rmse[step] = get_rms(zeta) / rms

        print("Max. normalised r.m.s. error for", labels[j], rmse.max())
        print("Mean of normalised r.m.s. error for", labels[j], rmse.mean())

        ax.plot(t, rmse, label=labels[j])

    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    pisqrt2 = np.pi/np.sqrt(2.0)
    ax.set_xticks([0.0, pisqrt2, 2.0*pisqrt2, 3.0*pisqrt2, 4.0*pisqrt2],
                  [r'$0$', r'$\pi/\sqrt{2}$', r'$2\pi/\sqrt{2}$',
                   r'$3\pi/\sqrt{2}$', r'$4\pi/\sqrt{2}$'])
    ax.set_xlabel(r'time, $t$')
    ax.set_ylabel(r'$(\zeta-\zeta^{\star})_{\mathrm{rms}} / \zeta^{\star}_{\mathrm{rms}}$')
    ax.legend(loc='best')

    for j, fname in enumerate(fnames):

        ncr.open(os.path.join(data_dir, fname))
        nsteps = ncr.get_num_steps()
        origin = ncr.get_box_origin()
        extent = ncr.get_box_extent()
        (nx, ny, nz) = ncr.get_box_ncells()
    
        for i, t in enumerate(times):
            ax = grid[i+j*3]

            colorbar = (i == 0)
    
            _, _, zeta, _, _, _, _ = get_exact(t, origin, extent, nx, ny, nz,
                                               kk=0.5, ll=0.5, mm=1.0, N=2.0, f=1.0,
                                               what=0.001, copy_periodic=True)

            rms = get_rms(zeta)

            zeta = get_dataset(ncr, 'z_vorticity', t, True) - zeta

            zeta = zeta / rms

            im, cbar = make_imshow(
                ax=ax,
                plane='xz',
                loc=loc,
                fdata=zeta,
                vmin=vmin,
                vmax=vmax,
                ncr=ncr,
                cmap='bwr',
                cmap_norm=None,
                clabel=r'$(\zeta - \zeta^{\star})/\zeta^{\star}_{\mathrm{rms}}$',
                colorbar=colorbar,
                xticks=xticks,
                xticklab=xticklab,
                yticks=yticks,
                yticklab=yticklab)

            xg, yg, zg = ncr.get_meshgrid()
            zeta = get_plane('xz', loc, zeta)
            ax.contour(xg[:, loc, :], zg[:, loc, :], zeta, 10, colors='black', linewidths=0.2)

            if j == 0:
                remove_xticks(ax)
                add_annotation(ax, r'$t=' + times_lab[i] + r'$', xy=(0.03, 1.25))
    
            if i == 0 or i == 3:
                pass
            else:
                remove_yticks(ax)
            
                add_annotation(grid[j*3], labels[j], xy=(-0.52, 1.0))

        ncr.close()

    add_annotation(grid[2], r'$y=0$', xy=(0.9, 1.7))

    plt.subplots_adjust(hspace=0.4)
 
    plt.savefig('iw_cross_sections_zeta.pdf', bbox_inches='tight')
    plt.close()
