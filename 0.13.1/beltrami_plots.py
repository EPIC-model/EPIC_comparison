import numpy as np
from tools.nc_reader import nc_reader
import matplotlib.pyplot as plt
from tools.mpl_style import *
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib as mpl
from tools.utils import *
import colorcet as cc
from tools.mpl_beautify import *
import os
from numpy.polynomial import polynomial as pp

mpl.rcParams['font.size'] = 13
mpl.rcParams['lines.linewidth'] = 1.5

ncr = nc_reader()

#
# Beltrami -- convergence of initial field error
#
data_dir = '../data/beltrami/'
file_bases = ['epic_beltrami_', 'epic_beltrami_', 'ps3d_beltrami_']
labels = [r'$\xi$', r'$\eta$', r'$\zeta$']
linestyles = ['solid', 'dashed', 'solid']
annotations = [r'EPIC (all)', r'EPIC (w/o surfaces)', r'PS3D']
xy = [(0.05, 1.07), (0.05, 1.07), (0.05, 1.07)]

grids = np.array([32, 48, 64, 96, 128])
xticks = [r'$32^3$', r'$48^3$', r'$64^3$', r'$96^3$', r'$128^3$']

def get_rms(data):
    return np.sqrt((data ** 2).mean())

n = len(grids)

ncr_ref = nc_reader()

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(9, 3), dpi=200)

for k, fb in enumerate(file_bases):
    rms = np.zeros((3, n))
    for j, grid in enumerate(grids):
        ncr.open(os.path.join(data_dir, fb + str(grid) + '_fields.nc'))
        ncr_ref.open(os.path.join(data_dir,
                                  'beltrami_' + str(grid) + 'x' + str(grid) + 'x' + str(grid) + '.nc'))
        for i, name in enumerate(['x_vorticity', 'y_vorticity', 'z_vorticity']):
            dset_ref = ncr_ref.get_dataset(step=0, name=name, copy_periodic=False)
            dset = ncr.get_dataset(step=0, name=name, copy_periodic=False)

            if 'w/o surfaces' in annotations[k]:
                rms[i, j] = get_rms(dset[:, :, 1:-1] - dset_ref[:, :, 1:-1])
            else:
                rms[i, j] = get_rms(dset - dset_ref)

            if 'epic' in fb:
                rms_ex = get_rms(dset[:, :, 1:-1] - dset_ref[:, :, 1:-1])
                print("EPIC", grid, labels[i], "rms error incl. bndry", round(rms[i, j], 4),
                      ", excl. bndry", round(rms_ex, 4),
                      "error reduction [%]", round((rms[i, j] - rms_ex) / rms[i, j] * 100.0, 2))
            
        ncr.close()
        ncr_ref.close()


    for i, label in enumerate(labels):
        if not k == 1:
            label = None
        if k == 2 and i == 2:
           ax = axs[k].twinx()
           ax.plot(grids, rms[i, :], label=label, color=colors[i], linestyle=linestyles[i], marker='o')
           ax.set_xscale('log', base=2)
           ax.set_yscale('log', base=10)
           ax.set_ylabel(r'r.m.s. error in $\zeta$')
           continue
        axs[k].plot(grids, rms[i, :], label=label, color=colors[i], linestyle=linestyles[i], marker='o')
        axs[k].set_xscale('log', base=2)
        axs[k].set_yscale('log', base=10)

    add_annotation(axs[k], annotations[k], xy=xy[k])
    axs[k].set_xlabel(r'grid resolution, $n_x\times n_y\times n_z$')
    axs[k].set_ylabel(r'r.m.s. error')
    axs[k].grid(which='both', zorder=-2, linestyle='dashed')

axs[0].plot(grids, 20.0 / grids ** 2, label=r'$\propto n_x^{-2}$', color='black', linestyle='dashed')
axs[1].plot(grids, 15.0 / grids ** 2, color='black', linestyle='dashed') 
axs[2].plot(grids, 0.5 / grids ** 2, color='black', linestyle='dashed')

for i in range(3):
    axs[i].set_xticks(grids, xticks)

fig.legend(loc='upper center', ncols=4, bbox_to_anchor=(0.5, 1.1))
fig.tight_layout()
plt.savefig('beltrami_initial_field_convergence.pdf', bbox_inches='tight')
plt.close()

ncr_ref = None


#
# Beltrami -- energy and enstrophy
#

mpl.rcParams['font.size'] = 13
mpl.rcParams['lines.linewidth'] = 1.5

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(8, 5), dpi=200, sharex=True, sharey='row')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

annotations = [r'EPIC (parcels)', r'EPIC (fields)', r'PS3D']


def get_slope(ts, ene):
    
    emid = 0.5 * (ene.max() + ene.min())
    elo = 0.75 * emid
    ehi = 1.25 * emid

    i = np.argmin(abs(ene - ehi))
    j = np.argmin(abs(ene - elo))

    x = ts[i:j]
    y = ene[i:j]
    
    # y = a * x + b
    [b, a] = pp.polyfit(x, y, deg=1)

    return a

for i, grid in enumerate([32, 48, 64, 96, 128]):
    ncr.open(os.path.join(data_dir, 'epic_beltrami_' + str(grid) + '_parcel_stats.nc'))
    t = ncr.get_all('t')
    ke = ncr.get_all('ke')
    ncr.close()
    get_slope(t, ke)
    print("Slope for", grid, "EPIC (parcels)", round(get_slope(t, ke), 4), end=' ')

    ncr.open(os.path.join(data_dir, 'epic_beltrami_' + str(grid) + '_field_stats.nc'))
    t = ncr.get_all('t')
    ke = ncr.get_all('ke')
    en = ncr.get_all('en')
    ncr.close()
    print("EPIC (fields)", round(get_slope(t, ke), 4), end=' ')

    ncr.open(os.path.join(data_dir, 'ps3d_beltrami_' + str(grid) + '_field_stats.nc'))
    t = ncr.get_all('t')
    ke = ncr.get_all('ke')
    en = ncr.get_all('en')
    ncr.close()
    print("PS3D", round(get_slope(t, ke), 4))


for i, grid in enumerate([32, 48, 64, 96, 128]):

    ncr.open(os.path.join(data_dir, 'epic_beltrami_' + str(grid) + '_parcel_stats.nc'))
    t = ncr.get_all('t')
    ke = ncr.get_all('ke')
    en = ncr.get_all('en')
    ncr.close()

    print("Enstrophy peak for", grid, "EPIC (parcels)", round(en.max(), 2), end=' ')
    
    axs[0, 0].plot(t, ke, color=colors[i])
    axs[1, 0].plot(t, en, color=colors[i])

    ncr.open(os.path.join(data_dir, 'epic_beltrami_' + str(grid) + '_field_stats.nc'))
    t = ncr.get_all('t')
    ke = ncr.get_all('ke')
    en = ncr.get_all('en')
    ncr.close()

    print("EPIC (fields)", round(en.max(), 2), end=' ')

    axs[0, 1].plot(t, ke, color=colors[i], label=r'$' + str(grid) + '^3$')
    axs[1, 1].plot(t, en, color=colors[i])
    
    ## PS3D version 0.0.6:
    #t, ke, en = np.loadtxt('../../ps3d_runs/paper_runs/beltrami_' + str(grid) + '_ecomp.asc',
    #                       unpack=True, skiprows=1)
    #ncelli = 1.0 / (grid ** 3)
    #voli = 1.0 / np.pi ** 3 
    #ke = ke * ncelli * voli
    #en = en * ncelli * voli

    # PS3D version 0.0.10:
    ncr.open(os.path.join(data_dir, 'ps3d_beltrami_' + str(grid) + '_field_stats.nc'))
    t = ncr.get_all('t')
    ke = ncr.get_all('ke')
    en = ncr.get_all('en')
    ncr.close()

    print("PS3D", round(en.max(), 2))

    axs[0, 2].plot(t, ke, color=colors[i])
    axs[1, 2].plot(t, en, color=colors[i])

for i in range(2):
    for j in range(3):
        axs[i, j].grid(linestyle='dashed', zorder=-1)

for j in range(3):
    remove_xticks(axs[0, j])
    axs[1, j].set_xlabel(r'time, $t$')
    add_annotation(axs[0, j], annotations[j], xy=(0.05, 1.07))

for j in range(1, 3):
    remove_yticks(axs[0, j])
    remove_yticks(axs[1, j])

axs[0, 0].set_ylabel(r'kinetic energy, $\mathcal{K}$')
axs[1, 0].set_ylabel(r'enstrophy, $\Upsilon$')

fig.legend(loc='upper center', ncols=5, bbox_to_anchor=(0.5, 1.07))
fig.tight_layout()
plt.savefig('beltrami_ke_en.pdf', bbox_inches='tight')
plt.close()


#
# Beltrami -- cross sections
#
mpl.rcParams['font.size'] = 12

fnames = ['epic_beltrami_128_fields.nc', 'ps3d_beltrami_128_fields.nc']
labels = [r'EPIC', r'PS3D']
times = [9.0, 10.0, 11.0]
planes = ['xy', 'xz', 'xy', 'xy']
locations = ['centre', 'centre', 'bottom', 'top']
vmin_limit = [None, None, 1.0e-2, 1.0e-2]
locs = [64, 64, 0, 128]
loc_labels = [r'$z = 0$', r'$y = 0$', r'$z=-\pi/2$', r'$z = \pi/2$']
cmap = cc.cm['rainbow4']

ticks = [-np.pi/2.0, -np.pi/4.0, 0.0, np.pi/4.0, np.pi/2.0]
ticklabs = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$'] 

# ts has times of "*_fields.nc" files (has lower write frequency than stats files)
def get_closest_steps(fn, ts):
    ncr2 = nc_reader()
    ncr2.open(fn)
    ke = ncr2.get_all('ke')
    tstats = ncr2.get_all('t')
    ncr2.close()
    steps = [0, 0, 0]
    for i, tt in enumerate(times):
        step = find_nearest(ts, tt)
        idx = find_nearest(tstats, tt)
        print("Field step", step, "is closest time", round(ts[step], 4), "to", tt,
              "Energy decay:", round((ke[0] - ke[idx]) / ke[0] * 100, 2), "(stat file step:", idx, ")")
        steps[i] = step
    return steps

for l, plane in enumerate(planes):

    fig = plt.figure(figsize=(8, 5), dpi=200)
    grid = ImageGrid(fig, 111,
                     nrows_ncols=(2, 3),
                     aspect=True,
                     axes_pad=(0.4, 0.3),
                     direction='row',
                     share_all=True,
                     cbar_location="right",
                     cbar_mode='single',
                     cbar_size="4%",
                     cbar_pad=0.05)

    vmin = 1000
    vmax = -1000
    for j, fname in enumerate(fnames):
        ncr.open(os.path.join(data_dir, fname))
        t = ncr.get_all('t')
        if 'epic' in fname:
            steps = get_closest_steps(os.path.join(data_dir, 'epic_beltrami_128_parcel_stats.nc'), t)
        else:
            steps = get_closest_steps(os.path.join(data_dir, 'ps3d_beltrami_128_field_stats.nc'), t)
        for i, step in enumerate(steps):
            vmag = ncr.get_dataset(step, name='vorticity_magnitude')
            vmag = get_plane(plane, locs[l], vmag)
            vmin = min(vmin, vmag.min())
            vmax = max(vmax, vmag.max())
        ncr.close()

    print("Global vmax:", vmax, "Global vmin:", vmin)

    if not vmin_limit[l] is None:
        vmin = vmin_limit[l]

    for j, fname in enumerate(fnames):
        ncr.open(os.path.join(data_dir, fname))

        t = ncr.get_all('t')

        if 'epic' in fname:
            steps = get_closest_steps(os.path.join(data_dir, 'epic_beltrami_128_parcel_stats.nc'), t)
        else:
            steps = get_closest_steps(os.path.join(data_dir, 'ps3d_beltrami_128_field_stats.nc'), t)

        for i, step in enumerate(steps):
            ax = grid[i+j*3]

            vmag = ncr.get_dataset(step, name='vorticity_magnitude')

            im, cbar = make_imshow(ax=ax,
                                   plane=plane,
                                   loc=locs[l],
                                   fdata=vmag,
                                   step=step,
                                   ncr=ncr,
                                   interpolation='bilinear',
                                   vmin=vmin,
                                   vmax=vmax,
                                   cmap_norm='log',
                                   cmap=cmap,
                                   colorbar=True)

            if i+j*3 < 3:
                remove_xticks(ax)
                #add_timestamp(ax, t[step], xy=(0.03, 1.06), fmt="%.2f")
            else:
                ax.set_xticks(ticks, ticklabs)

            if i == 0: # or i == 3:
                ax.set_yticks(ticks, ticklabs)
            else:
                remove_yticks(ax)

            add_timestamp(ax, t[step], xy=(0.03, 1.06), fmt="%.2f")

        ncr.close()


    cbar.set_label(r'vorticity magnitude, $|\bm{\omega}|$')
    add_annotation(grid[2], loc_labels[l], xy=(0.8, 1.2))

    add_annotation(grid[0], labels[0], xy=(-0.6, 1.06))
    add_annotation(grid[3], labels[1], xy=(-0.6, 1.06))

    plt.savefig('beltrami_vmag_cs_plane_' + plane + '_loc_' + locations[l] + '.pdf',
                bbox_inches='tight', dpi=200)
    plt.close()
