import numpy as np
from tools.nc_reader import nc_reader
import matplotlib.pyplot as plt
from tools.mpl_style import *
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib as mpl
from tools.utils import *
import colorcet as cc
from tools.mpl_beautify import *

#
# Beltrami -- energy and enstrophy
#

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.5

ncr = nc_reader()

# 14 Nov 2022
# https://stackoverflow.com/a/41765095
from matplotlib.legend_handler import HandlerBase
class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle,
                       x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0,y0+width], [0.7*height,0.7*height], color=orig_handle[0])
        l2 = plt.Line2D([x0,y0+width], [0.3*height,0.3*height], color=orig_handle[0],
                        linestyle=orig_handle[1])
        return [l1, l2]

class MyHandler(HandlerBase):
    def create_artists(self, legend, orig_handle,
                       x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0,y0+width], [0.5*height,0.5*height], color=orig_handle[0],
                        linestyle=orig_handle[1])
        l2 = plt.Line2D([x0,y0+width], [0.1*height,0.1*height], color=orig_handle[2],
                        linestyle=orig_handle[3])
        return [l1, l2]

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(9, 4), dpi=200)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
linestyles = ['solid', 'dashed', 'dashdot']

for i, grid in enumerate([32, 64]):

    ncr.open('epic_beltrami_' + str(grid) + '_parcel_stats.nc')
    t = ncr.get_all('t')
    ke = ncr.get_all('ke')
    en = ncr.get_all('en')
    ncr.close()
    
    axs[0].plot(t, ke, color=colors[0], linestyle=linestyles[i])
    axs[1].plot(t, en, color=colors[0], linestyle=linestyles[i])

    ## PS3D version 0.0.6:
    #t, ke, en = np.loadtxt('../../ps3d_runs/paper_runs/beltrami_' + str(grid) + '_ecomp.asc',
    #                       unpack=True, skiprows=1)
    #ncelli = 1.0 / (grid ** 3)
    #voli = 1.0 / np.pi ** 3 
    #ke = ke * ncelli * voli
    #en = en * ncelli * voli

    # PS3D version 0.0.10:
    ncr.open('ps3d_beltrami_' + str(grid) + '_field_stats.nc')
    t = ncr.get_all('t')
    ke = ncr.get_all('ke')
    en = ncr.get_all('en')
    ncr.close()

    axs[0].plot(t, ke, color=colors[1], linestyle=linestyles[i])
    axs[1].plot(t, en, color=colors[1], linestyle=linestyles[i])

axs[0].grid(linestyle='dashed', zorder=-1)
axs[1].grid(linestyle='dashed', zorder=-1)

axs[0].set_xlabel(r'$t$')
axs[1].set_xlabel(r'$t$')

axs[0].legend([(colors[0], 'dashed'), (colors[1], 'dashed'),
               [colors[0], 'solid', colors[1], 'solid'],
               [colors[0], 'dashed', colors[1], 'dashed']],
              [r'EPIC', r'PS3D', r'$32^3$' , r'$64^3$'],
              handler_map={tuple: AnyObjectHandler(), list: MyHandler()},
              loc='upper center', bbox_to_anchor=(1.1, 1.15), ncols=4)

axs[0].set_ylabel(r'kinetic energy, $\mathcal{K}$')
axs[1].set_ylabel(r'enstrophy, $\Upsilon$')

plt.savefig('beltrami_ke_en.pdf', bbox_inches='tight')
plt.close()


#
# Beltrami -- cross sections
#

fnames = ['epic_beltrami_64_fields.nc', 'ps3d_beltrami_64_fields.nc']
labels = [r'EPIC', r'PS3D']
fignames = ['beltrami_epic_vmag_cs.pdf', 'beltrami_ps3d_vmag_cs.pdf']
ke_decay = [0.95, 0.85, 0.75]
plane = 'xy'
loc = 32
loc_label = r'$z = 0$'
cmap = cc.cm['rainbow4_r']

ticks = [-np.pi/2.0, -np.pi/4.0, 0.0, np.pi/4.0, np.pi/2.0]
ticklabs = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$'] 

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

def get_closest_steps(fn, times):
    ncr2 = nc_reader()
    ncr2.open(fn)
    ts = ncr2.get_all('t')
    ke = ncr2.get_all('ke')
    ncr2.close()
    steps = [0, 0, 0]
    for i, d in enumerate(ke_decay):
        idx = find_nearest(ke, d * ke[0])
        step = find_nearest(times, ts[idx])
        print(times[step], ts[idx])
        steps[i] = step
    return steps

vmin = 1000
vmax = -1000
for j, fname in enumerate(fnames):
    ncr.open(fname)
    t = ncr.get_all('t')
    if 'epic' in fname:
        steps = get_closest_steps('epic_beltrami_64_parcel_stats.nc', t)
    else:
        steps = get_closest_steps('ps3d_beltrami_64_field_stats.nc', t)
    for i, step in enumerate(steps):
        vmag = ncr.get_dataset(step, name='vorticity_magnitude')
        vmag = get_plane(plane, loc, vmag)
        vmin = min(vmin, vmag.min())
        vmax = max(vmax, vmag.max())
    ncr.close()

for j, fname in enumerate(fnames):
    ncr.open(fname)

    t = ncr.get_all('t')

    if 'epic' in fname:
        steps = get_closest_steps('epic_beltrami_64_parcel_stats.nc', t)
    else:
        steps = get_closest_steps('ps3d_beltrami_64_field_stats.nc', t)

    for i, step in enumerate(steps):
        ax = grid[i+j*3]

        vmag = ncr.get_dataset(step, name='vorticity_magnitude')

        im, cbar = make_imshow(ax=ax,
                               plane=plane,
                               loc=loc,
                               fdata=vmag,
                               step=step,
                               ncr=ncr,
                               vmin=vmin,
                               vmax=vmax,
                               cmap_norm='log',
                               cmap=cmap,
                               colorbar=True)

        if i+j*3 < 3:
            remove_xticks(ax)
        else:
            ax.set_xticks(ticks, ticklabs)

        if i == 0: # or i == 3:
            ax.set_yticks(ticks, ticklabs)
        else:
            remove_yticks(ax)

        add_timestamp(ax, t[step], xy=(0.03, 1.06), fmt="%.2f")

    ncr.close()


cbar.set_label(r'vorticity magnitude, $|\bm{\omega}|$')
add_annotation(grid[2], loc_label, xy=(0.8, 1.2))

add_annotation(grid[0], labels[0], xy=(-0.6, 1.06))
add_annotation(grid[3], labels[1], xy=(-0.6, 1.06))

plt.savefig('beltrami_vmag_cs.pdf', bbox_inches='tight', dpi=200)
plt.close()


#
# Restart cross sections
#
fnames = ['restart/epic_restart_64_fields.nc', 'ps3d_beltrami_64_fields.nc']
labels = [r'EPIC', r'PS3D']
plane = 'xy'
loc = 32
loc_label = r'$z = 0$'
cmap = cc.cm['rainbow4_r']

ticks = [-np.pi/2.0, -np.pi/4.0, 0.0, np.pi/4.0, np.pi/2.0]
ticklabs = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']

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

ncrestart = nc_reader()
ncrestart.open("restart/restart_step_60.nc")
restart_time = ncrestart.get_all('t')[0]
ncrestart.close()

for j, fname in enumerate(fnames):
    ncr.open(fname)
    t = ncr.get_all('t')
    if 'epic' in fname:
        t = t + restart_time
        steps = [3, 4, 5]
    else:
        steps = [63, 64, 65]
    for i, step in enumerate(steps):
        vmag = ncr.get_dataset(step, name='vorticity_magnitude')
        vmag = get_plane(plane, loc, vmag)
        vmin = min(vmin, vmag.min())
        vmax = max(vmax, vmag.max())
    ncr.close()


for j, fname in enumerate(fnames):
    ncr.open(fname)

    t = ncr.get_all('t')

    if 'epic' in fname:
        t = t + restart_time
        steps = [3, 4, 5]
    else:
        steps = [63, 64, 65]

    for i, step in enumerate(steps):
        ax = grid[i+j*3]

        vmag = ncr.get_dataset(step, name='vorticity_magnitude')

        im, cbar = make_imshow(ax=ax,
                               plane=plane,
                               loc=loc,
                               fdata=vmag,
                               step=step,
                               ncr=ncr,
                               vmin=vmin,
                               vmax=vmax,
                               cmap_norm='log',
                               cmap=cmap,
                               colorbar=True)

        if i+j*3 < 3:
            remove_xticks(ax)
        else:
            ax.set_xticks(ticks, ticklabs)

        if i == 0:
            ax.set_yticks(ticks, ticklabs)
        else:
            remove_yticks(ax)

        add_timestamp(ax, t[step], xy=(0.03, 1.06), fmt="%.2f")

    ncr.close()

cbar.set_label(r'vorticity magnitude, $|\bm{\omega}|$')
add_annotation(grid[2], loc_label, xy=(0.8, 1.2))

add_annotation(grid[0], labels[0], xy=(-0.6, 1.06))
add_annotation(grid[3], labels[1], xy=(-0.6, 1.06))

plt.savefig('beltrami_vmag_cs_restart.pdf', bbox_inches='tight')
plt.close()
