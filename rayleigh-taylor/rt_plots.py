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

units['time'] = None

ncr = nc_reader()

#
# RT -- rms volume error
#
ncr.open('epic_rt_field_stats.nc')

t = ncr.get_all('t')
rms_v = ncr.get_all('rms_v')

ncr.close()

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.5 
plt.figure(figsize=(8, 3), dpi=200)
plt.grid(which='both', linestyle='dashed')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.plot(t, rms_v)
plt.xlabel('time, $t$')
plt.ylabel('r.m.s.\ volume error')
plt.tight_layout()
plt.savefig('rt_rms_v_error.pdf', bbox_inches='tight')
plt.close()

#
# RT -- energy and enstrophy plot
#
mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.5

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 4), dpi=200, sharex=True)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
linestyles = ['solid', 'dashed']

labels = [r'EPIC', r'PS3D']
fnames = ['epic_rt_parcel_stats.nc', 'ps3d_rt_field_stats.nc']

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

fig.legend(loc='upper center', ncols=2, bbox_to_anchor=(0.5, 1.07))

plt.tight_layout()

plt.savefig('rt_ke_ape_en.pdf', bbox_inches='tight')
plt.close()


#
# RT -- buoyancy extrema plot
#

fnames = ['epic_rt_parcel_stats.nc', 'epic_rt_field_stats.nc', 'ps3d_rt_field_stats.nc']
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
steps = [0, 4, 5, 6, 10]

fnames = ['epic_rt_fields.nc', 'ps3d_rt_fields.nc']
labels = [r'EPIC', r'PS3D']

fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(10, 3), dpi=200, sharey=True)

mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['font.size'] = 12

for i, fname in enumerate(fnames):
    ncr.open(fname)

    t = ncr.get_all('t')
    z = ncr.get_all('z')

    for j, step in enumerate(steps):
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

axs[2].legend(loc='upper center', bbox_to_anchor=(0.5, 1.32), ncols=2)

plt.savefig('rt_buoyancy_mean_profile.pdf', bbox_inches='tight')
plt.close()

#
# RT -- cross sections
#

mpl.rcParams['font.size'] = 12
fmames = ['epic_rt_fields.nc', 'ps3d_rt_fields.nc']
labels = [r'EPIC', r'PS3D']
fignames = ['rt_epic_buoyancy_cs.pdf', 'rt_ps3d_buoyancy_cs.pdf']
steps = [0, 2, 4, 5, 6, 10]
plane = 'xz'
loc = 32
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

    for i, step in enumerate(steps):
        ax = grid[i]

        colorbar = (i == 0)

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
                               colorbar=colorbar,
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
