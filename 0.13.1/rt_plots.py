import numpy as np
import matplotlib.pyplot as plt
from tools.nc_reader import nc_reader
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
from tools.utils import *
from tools.mpl_style import *
import argparse
import os
from tools.mpl_beautify import *
import matplotlib as mpl
import colorcet as cc
from tools.units import *
from linear_interpl import time_interpl, find_bounds, get_dataset

units['time'] = None

ncr = nc_reader()


data_dir = '/local_raid2/mf248/data/rt'

#
# RT -- buoyancy and zeta r.m.s. plots compared to coarsend 384^3 simulation
#    
grids = np.array([32, 48, 64, 96, 128])
gridlabs = [r'$32^3$', r'$48^3$', r'$64^3$', r'$96^3$', r'$128^3$']
vmins = [10, 20, 30, 40]

mpl.rcParams['font.size'] = 16
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

ncr_crse = nc_reader()

# assumes that periodic layers are not copied;
# this rms counts the vertical boundaries together as 1.
def get_rt_rms(data):
    # nz is the number of grid points and not cells!
    nx, ny, nz = data.shape

    rms = (data[:, :, 1:nz-1] ** 2).sum()
    rms = rms + 0.5 * (data[:, :, 0] ** 2).sum()
    rms = rms + 0.5 * (data[:, :, nz-1] ** 2).sum()
    return np.sqrt(rms / (nx * ny * (nz-1)))

t_refs = [4.0, 6.0, 10.0]

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(11.5, 6), dpi=200, sharex=True)

for k, name in enumerate(['buoyancy', 'z_vorticity']):

    for l, t_ref in enumerate(t_refs):

        for i, vmin in enumerate(vmins):

            rms = np.zeros(len(grids))
        
            for j, grid in enumerate(grids):
                ncr_crse.open(os.path.join(data_dir,
                                           'coarsened',
                                           'crse_' + str(grid) + '_epic_rt_384_fields.nc'))

                dset_crse = get_dataset(ncr_crse, name, t_ref, copy_periodic=False)
            
                ncr_crse.close()
            
                ncr.open(os.path.join(data_dir,
                                      'epic_rt_' + str(grid) + '_vmin_' + str(vmin) + '_fields.nc'))
            
                dset = get_dataset(ncr, name, t_ref, copy_periodic=False)
            
                ncr.close()
            
                rms[j] = get_rt_rms(dset - dset_crse)

            lab = None
            if k == 0 and l == 0:
                lab = r'$\tilde{V}_{\min} = 1/' + str(vmin) + '$'

            axs[k, l].plot(grids, rms, label=lab, marker='o')

for k in range(2):
    for l in range(len(t_refs)):
        axs[k, l].set_xscale('log', base=2)
        axs[k, l].set_yscale('log', base=10)
        axs[k, l].grid(which='both', linestyle='dashed', zorder=-10)
        
        if k < 1:
            remove_xticks(axs[k, l])
        else:
            axs[k, l].set_xticks(grids, gridlabs)
            axs[k, l].set_xlabel(r'grid resolution, $n_x\times n_y\times n_z$')

for l in range(len(t_refs)):
    add_timestamp(axs[0, l], t_refs[l], xy=(0.03, 1.06), fmt="%.1f") 


l1, = axs[0, 0].plot(grids, 8.0/grids**1.2, linestyle='dashed', color='black')
l2, = axs[1, 0].plot(grids, 17.0/grids**1.2, linestyle='dashed', color='black')
axs[0, 0].legend([l1], [r'$\propto n_z^{-1.2}$'], loc='upper right')
axs[1, 0].legend([l2], [r'$\propto n_z^{-1.2}$'], loc='upper right')

axs[0, 0].set_ylabel(r'r.m.s. buoyancy error') #, $b_{\mathrm{rms}}$')
axs[1, 0].set_ylabel(r'r.m.s. z-vorticity error') #, $\zeta_{\mathrm{rms}}$')

fig.legend(loc='upper center', ncols=len(vmins), bbox_to_anchor=(0.55, 1.08))

fig.tight_layout()
plt.savefig('rt_buoyancy_zeta_rms.pdf', bbox_inches='tight')
plt.close()

ncr_crse = None


#
# RT -- energy convergence:
#
t_ref = 4.0

ape_ref = 4.0 / np.pi
ke_ref = 0.0

te_ref = ape_ref + ke_ref

grids = ['32', '48', '64', '96', '128', '384']

for prog in ['epic', 'ps3d']:
    print(prog.upper())
    for i, grid in enumerate(grids):
        if prog == 'epic':
            if grid == '384':
                fname = os.path.join(data_dir,
                                     prog + '_rt_' + grid + '_parcel_stats.nc')
            else:
                fname = os.path.join(data_dir,
                                     prog + '_rt_' + grid + '_vmin_20_parcel_stats.nc')
        else:
            fname = os.path.join(data_dir, prog + '_rt_' + grid + '_field_stats.nc')

        if not os.path.exists(fname):
            print("File does not exist", grid)
            continue
            
        ncr.open(fname)
        t = ncr.get_all('t')
        lo, hi = find_bounds(t_ref, t)
        if prog == 'ps3d':
            te = ncr.get_all('ke') + ncr.get_all('ape')
        else:
            te = ncr.get_all('te')
        te_init = te[0]
        te = time_interpl(t_ref, te[lo], t[lo], te[hi], t[hi])
        ncr.close()
        loss = abs(te - te_init) / te_init * 100.0
        print("grid", grid, "energy loss in TE (percent)", loss, "rounded:", round(loss, 2))

#
# RT -- energy and enstrophy plot
#
mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.5

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 4), dpi=200, sharex=True)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
linestyles = ['solid', 'dashed']

labels = [r'EPIC', r'PS3D']
fnames = [
    os.path.join(data_dir, 'epic_rt_128_vmin_20_parcel_stats.nc'),
    os.path.join(data_dir, 'ps3d_rt_128_field_stats.nc')
]

linestyles = ['solid', 'dashed']

for i, fname in enumerate(fnames):

    ncr.open(fname)
    t = ncr.get_all('t')
    ke = ncr.get_all('ke')
    ape = ncr.get_all('ape')
    if 'ps3d' in fname:
        te = ke + ape
    else:
        te = ncr.get_all('te')
    en = ncr.get_all('en')
    ncr.close()

    axs[0, 0].plot(t, ke, color=colors[i], linestyle=linestyles[i], label=labels[i])
    axs[0, 1].plot(t, ape, color=colors[i], linestyle=linestyles[i])
    axs[1, 0].plot(t, te, color=colors[i], linestyle=linestyles[i])
    axs[1, 1].plot(t, en, color=colors[i], linestyle=linestyles[i])

remove_xticks(axs[0, 0])
remove_xticks(axs[0, 1])
axs[0, 0].grid(linestyle='dashed', zorder=-1)
axs[0, 1].grid(linestyle='dashed', zorder=-1)
axs[1, 0].grid(linestyle='dashed', zorder=-1)
axs[1, 1].grid(linestyle='dashed', zorder=-1)

axs[1, 0].set_xlabel(r'time, $t$')
axs[1, 1].set_xlabel(r'time, $t$')

axs[0, 0].set_ylabel(r'kinetic energy, $\mathcal{K}$')
axs[0, 1].set_ylabel(r'avail. pot. energy, $\mathcal{A}$')
axs[1, 0].set_ylabel(r'total energy, $\mathcal{T}$')
axs[1, 1].set_ylabel(r'enstrophy, $\Upsilon$')

fig.legend(loc='upper center', ncols=3, bbox_to_anchor=(0.5, 1.07))

fig.tight_layout()

plt.savefig('rt_ke_ape_en.pdf', bbox_inches='tight')
plt.close()

#
# RT -- buoyancy extrema plot
#
fnames = [
    os.path.join(data_dir, 'epic_rt_128_vmin_20_parcel_stats.nc'),
    os.path.join(data_dir, 'epic_rt_128_vmin_20_field_stats.nc'),
    os.path.join(data_dir, 'ps3d_rt_128_field_stats.nc')
]


labels = [r'EPIC', r'PS3D']

fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(9, 4), dpi=200, sharex=True)

mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['font.size'] = 14

for i, fname in enumerate(fnames):
    ncr.open(fname)

    t = ncr.get_all('t')

    bmin = ncr.get_all('min_buoyancy')
    bmax = ncr.get_all('max_buoyancy')

    j = 1
    k = 0
    llabel = [None] * 2
    if i == 0:
        j = 0
        k = 2
        llabel = [r'$b_{\min}$', r'$b_{\max}$']
    if i == 1:
        j = 0
        llabel = [r'$\bar{b}_{\min}$', r'$\bar{b}_{\max}$']

    ncr.close()

    axs[j].plot(t, bmin, color=colors[0+k], label=llabel[0])
    axs[j].plot(t, bmax, color=colors[1+k], label=llabel[1])
    if i > 0:
        axs[j].axhline(-1, linestyle='dashed', color='black', zorder=-1)
        axs[j].axhline(1, linestyle='dashed', color='black', zorder=-1)
        axs[j].grid(linestyle='dashed', zorder=-2)
        axs[j].set_ylabel(r'buoyancy')
        axs[j].set_xlim([0, 10])

        bbox = dict(boxstyle="round", facecolor="wheat", linewidth=0.5)
        axs[j].annotate(
            labels[j], xy=(0.02, 0.78), xycoords="axes fraction", bbox=bbox
        )

    if i == 2:
        axs[1].set_xlabel(r'time, $t$')

remove_xticks(axs[0])
axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.35), ncols=4)
plt.subplots_adjust(hspace=0.08)

plt.savefig('rt_buoyancy_extrema.pdf', bbox_inches='tight')
plt.close()


#
# RT -- buoyancy mean profile
#
steps = np.zeros((2, 5), dtype='int')
steps[0, :] = np.array([0, 62, 270, 600, 1000])
steps[1, :] = np.array([0, 4, 5, 6, 10])

fnames = [
    os.path.join(data_dir, 'epic_rt_128_vmin_20_fields.nc'),
    os.path.join(data_dir, 'ps3d_rt_128_fields.nc')
]

labels = [r'EPIC', r'PS3D']

fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(12, 3), dpi=200, sharey=True)

mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['font.size'] = 14

for i, fname in enumerate(fnames):
    ncr.open(fname)

    t = ncr.get_all('t')
    z = ncr.get_all('z')

    for j, step in enumerate(steps[i, :]):

        buoyg = ncr.get_dataset(step=step, name='buoyancy', copy_periodic=False)

        prof = buoyg.mean(axis=(0, 1))

        if j == 2:
            label = labels[i]
        else:
            label = None

        axs[j].plot(prof, z, color=colors[i], label=label)

#        print(t[step])

        if i == 0:
            axs[j].grid(linestyle='dashed', zorder=-2)
            axs[j].set_xlabel(r'$\langle\bar{b}\rangle$')
            if j == 0:
                axs[j].set_ylabel(r'$z$')
            add_timestamp(axs[j], t[step], xy=(0.03, 1.06), fmt="%.1f")

        if j > 0:
            remove_yticks(axs[j])

    ncr.close()

yticks = [-np.pi/2.0, -np.pi/4.0, 0.0, np.pi/4.0, np.pi/2.0]
yticklabs = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']
axs[0].set_yticks(yticks, yticklabs)

axs[2].legend(loc='upper center', bbox_to_anchor=(0.5, 1.34), ncols=2)

plt.savefig('rt_buoyancy_mean_profile.pdf', bbox_inches='tight')
plt.close()


#
# RT -- cross sections
#

mpl.rcParams['font.size'] = 12
fnames = [
    os.path.join(data_dir, 'epic_rt_128_vmin_20_fields.nc'),
    os.path.join(data_dir, 'ps3d_rt_128_fields.nc')
]

labels = [r'EPIC', r'PS3D']
fignames = ['rt_epic_buoyancy_cs.pdf', 'rt_ps3d_buoyancy_cs.pdf']
steps = np.zeros((2, 6), dtype='int')
steps[0, :] = np.array([0, 20, 62, 270, 600, 1000])
steps[1, :] = np.array([0, 3, 4, 5, 6, 10])
plane = 'xz'
loc = 64
loc_label = r'$y = 0$'
cmap = cc.cm['coolwarm']
cmap_norm = None
vmin = -1.0
vmax = 1.0

ticks = [-np.pi/2.0, -np.pi/4.0, 0.0, np.pi/4.0, np.pi/2.0]
ticklabs = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']

for j, fname in enumerate(fnames):
    ncr.open(fname)

    t = ncr.get_all('t')

    cbar_ext = (labels[j] == r'PS3D')

    fig = plt.figure(figsize=(8, 5), dpi=200)
    grid = ImageGrid(fig, 111,
                     nrows_ncols=(2, 3),
                     aspect=True,
                     axes_pad=(0.38, 0.3),
                     direction='row',
                     share_all=True,
                     cbar_location="right",
                     cbar_mode='single',
                     cbar_size="4%",
                     cbar_pad=0.05)

    for i, step in enumerate(steps[j, :]):

        ax = grid[i]

        buoyg = ncr.get_dataset(step, name='buoyancy')

        im, cbar = make_imshow(ax=ax,
                               plane=plane,
                               loc=loc,
                               fdata=buoyg,
                               step=step,
                               ncr=ncr,
                               cmap=cmap,
                               cmap_norm=cmap_norm,
                               vmin=vmin,
                               vmax=vmax,
                               colorbar=(i == 0), # this is fine because we explicity set vmin and vmax
                               cbar_ext=cbar_ext)

        if i < 3:
            remove_xticks(ax)
        else:
            ax.set_xticks(ticks, ticklabs)

        if i == 0 or i == 3:
            ax.set_yticks(ticks, ticklabs)
        else:
            remove_yticks(ax)

        add_timestamp(ax, t[step], xy=(0.03, 1.06), fmt="%.2f")

        if not cbar is None:
            cbar.set_label(r'buoyancy')

    add_annotation(grid[2], loc_label, xy=(0.8, 1.2))

    add_annotation(grid[0], labels[j], xy=(-0.3, 1.2))


    plt.savefig(fignames[j], bbox_inches='tight')
    plt.close()

    ncr.close()

#
# RT -- rms volume error for different lambda_max
#
fnames = ['epic_rt_64_vmin_20_lmax_3_field_stats.nc',
          'epic_rt_64_vmin_20_field_stats.nc',
          'epic_rt_64_vmin_20_lmax_5_field_stats.nc']

labels = [
    r'$\lambda_{\max} = 3$',
    r'$\lambda_{\max} = 4$',
    r'$\lambda_{\max} = 5$'
]

plt.figure(figsize=(8, 3), dpi=200)
mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.5
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

for i, fname in enumerate(fnames):
    ncr.open(os.path.join(data_dir, fname))

    t = ncr.get_all('t')
    rms_v = ncr.get_all('rms_v')

    ncr.close()

    plt.plot(t, rms_v, label=labels[i])

plt.grid(which='both', linestyle='dashed')
plt.xlabel('time, $t$')
plt.ylabel('r.m.s.\ volume error')
plt.legend(loc='upper center', ncol=len(fnames), bbox_to_anchor=(0.5, 1.27), columnspacing=0.95)
plt.tight_layout()
plt.savefig('rt_64_rms_v_error_lambda.pdf', bbox_inches='tight')
plt.close()


#
# RT -- rms volume error for different vmin
#
grids = ['32', '48', '64', '96', '128']

for grid in grids:
    fnames = ['epic_rt_' + grid + '_vmin_10_field_stats.nc',
              'epic_rt_' + grid + '_vmin_20_field_stats.nc',
              'epic_rt_' + grid + '_vmin_30_field_stats.nc',
              'epic_rt_' + grid + '_vmin_40_field_stats.nc']

    labels = [
        r'$\tilde{V}_{\min} = 1/10$',
        r'$\tilde{V}_{\min} = 1/20$',
        r'$\tilde{V}_{\min} = 1/30$',
        r'$\tilde{V}_{\min} = 1/40$'
    ]

    plt.figure(figsize=(8, 3), dpi=200)
    mpl.rcParams['font.size'] = 12
    mpl.rcParams['lines.linewidth'] = 1.5
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    for i, fname in enumerate(fnames):
        ncr.open(os.path.join(data_dir, fname))

        t = ncr.get_all('t')
        rms_v = ncr.get_all('rms_v')

        ncr.close()

        plt.plot(t, rms_v, label=labels[i])

    plt.grid(which='both', linestyle='dashed')
    plt.xlabel('time, $t$')
    plt.ylabel('r.m.s.\ volume error')
    plt.legend(loc='upper center', ncol=len(fnames), bbox_to_anchor=(0.5, 1.27), columnspacing=0.95)
    plt.tight_layout()
    plt.savefig('rt_' + grid + '_rms_v_error_vmin.pdf', bbox_inches='tight')
    plt.close()


#
# RT -- rms volume error across grid resolutions
#
fnames = ['epic_rt_32_vmin_20_field_stats.nc',
          'epic_rt_48_vmin_20_field_stats.nc',
          'epic_rt_64_vmin_20_field_stats.nc',
          'epic_rt_96_vmin_20_field_stats.nc',
          'epic_rt_128_vmin_20_field_stats.nc',
          'epic_rt_384_field_stats.nc']

labels = [r'$32^3$',
          r'$48^3$',
          r'$64^3$',
          r'$96^3$',
          r'$128^3$',
          r'$384^3$']

plt.figure(figsize=(8, 3), dpi=200)
mpl.rcParams['font.size'] = 13
mpl.rcParams['lines.linewidth'] = 1.5
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

for i, fname in enumerate(fnames):
    ncr.open(os.path.join(data_dir, fname))

    t = ncr.get_all('t')
    rms_v = ncr.get_all('rms_v')

    ncr.close()

    plt.plot(t, rms_v, label=labels[i])

plt.grid(which='both', linestyle='dashed')
plt.xlabel('time, $t$')
plt.ylabel('r.m.s.\ volume error')
plt.legend(loc='upper center', ncol=len(fnames), bbox_to_anchor=(0.5, 1.27), columnspacing=0.95)
plt.tight_layout()
plt.savefig('rt_rms_v_error_grid_resolution.pdf', bbox_inches='tight')
plt.close()


exit()

#
# RT -- histogram plots
#
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.size'] = 13


nbins = 100
dmin = -1.0
dmax =  1.0
as_density = False
name = 'buoyancy'

# returns hist and bin_centres
def get_histogram(x, bins, lo, hi, density=False):
    hist, bin_edges = np.histogram(x.flatten(), bins=bins, range=(lo, hi), density=density)
    # 7 August 2023
    # https://stackoverflow.com/a/72689634
    bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
    return hist, bin_centres

def do_plot(files, labels, times, fname, **kwargs):
    nc_ref = nc_reader()
    nc_run = nc_reader()

    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(9, 4), dpi=200, sharex=True, sharey=False)
    fig.tight_layout() # 10 August 2023, https://stackoverflow.com/a/9827848
    grid = axs.flatten()

    try:

        if not len(times) == 6:
            raise ValueError("6 times must be provided!")

        for k, f in enumerate(files):

            print("Processing", f)

            nc_run.open(os.path.join(data_dir, f))



            for i, t in enumerate(times):

                print("Time", t)

                lab = None
                if i == 1:
                    lab = labels[k]

                ax = grid[i]

                dset = get_dataset(nc_run, name, t)

                hist, bin_centres = get_histogram(dset, nbins, dmin, dmax, as_density)
                ax.step(x=bin_centres, y=hist / hist.sum(), where='mid', label=lab)
                ax.set_yscale('log')

            nc_run.close()

        for i, t in enumerate(times):

            ax = grid[i]

            ax.grid(which='both', linestyle='dashed', zorder=-10)

            add_timestamp(ax, t, xy=(0.03, 1.08), fmt="%.2f")

            if i < 3:
                remove_xticks(ax)
            else:
                ax.set_xlabel(r'buoyancy, $\bar{b}$')

            if i == 0 or i == 3:
                ax.set_ylabel(r'normalised bin count')

        grid[1].legend(loc='upper center', ncol=len(files), bbox_to_anchor=(0.5, 1.52))
        plt.subplots_adjust(wspace=0.34)
        plt.savefig(fname + '.pdf', bbox_inches='tight')

    except Exception as err:
        print(err)





time = [5, 5.5, 6, 7.0, 8.0, 10.0]

files = [
    'epic_rt_64_vmin_10_fields.nc',
    'epic_rt_64_vmin_20_fields.nc',
    'epic_rt_64_vmin_30_fields.nc',
    'epic_rt_64_vmin_40_fields.nc',
    'extracted_steps_from_epic_rt_384_fields.nc'
    ]

labels = [
    r'$\tilde{V}_{\min} = 1/10$',
    r'$\tilde{V}_{\min} = 1/20$',
    r'$\tilde{V}_{\min} = 1/30$',
    r'$\tilde{V}_{\min} = 1/40$',
    r'$384^3$'
]

do_plot(files, labels, time, 'rt_buoyancy_histogram_64_vmin')


files = [
    'epic_rt_32_vmin_20_fields.nc',
    'epic_rt_48_vmin_20_fields.nc',
    'epic_rt_64_vmin_20_fields.nc',
    'epic_rt_96_vmin_20_fields.nc',
    'epic_rt_128_vmin_20_fields.nc',
    'extracted_steps_from_epic_rt_384_fields.nc'
    ]

labels = [
    r'$32^3$',
    r'$48^3$',
    r'$64^3$',
    r'$96^3$',
    r'$128^3$',
    r'$384^3$'
]

do_plot(files, labels, time, 'rt_buoyancy_histogram_grid_resolution')


#
# RT -- number of parcels per cell
#
files = [
    'epic_rt_64_vmin_10_field_stats.nc',
    'epic_rt_64_vmin_15_field_stats.nc',
    'epic_rt_64_vmin_20_field_stats.nc',
    'epic_rt_64_vmin_30_field_stats.nc',
    'epic_rt_64_vmin_40_field_stats.nc',
    'epic_rt_64_vmin_40_lmax_3_field_stats.nc',
    'epic_rt_64_vmin_40_lmax_5_field_stats.nc',
    'epic_rt_64_vmin_20_lmax_3_field_stats.nc',
    'epic_rt_64_vmin_20_lmax_5_field_stats.nc'
]


print("File, number of samples, average number of parcels per cell:")
for f in files:
    ncr.open(os.path.join(data_dir, f))
    avg_npar = ncr.get_all('avg_npar')
#    min_npar = ncr.get_all('min_npar')
#    max_npar = ncr.get_all('max_npar')

    print(f, len(avg_npar), avg_npar.mean()) #, min_npar.min(), max_npar.max())

    ncr.close()

