import argparse
import os
from tools.nc_reader import nc_reader
from mpl_toolkits.axes_grid1 import ImageGrid
from tools.mpl_style import *
from tools.units import units
from tools.utils import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import colorcet as cc
import matplotlib.colors as cls

import netCDF4 as nc

try:

    def do_zoom(ncr, step, dset, loc, plane, ax, bounds, xlo, xhi, ylo, yhi, fname, **kwargs):
        
        # inset axes....
        axins = ax.inset_axes(bounds, alpha=1.0)

        show_axes = kwargs.pop('show_axes', False)
        scale = kwargs.pop('scale', 1.0)

        ncfile = nc.Dataset(filename, "r", format="NETCDF4")
        data = np.array(ncfile.variables['humidity']).squeeze()
        data = data * scale
        axins.pcolormesh(np.array(ncfile.variables['X']).squeeze(),
                         np.array(ncfile.variables['Z']).squeeze(),
                         data,
                         vmin=vmin,
                         vmax=vmax,
                         cmap=my_cmap,
                         shading='gouraud',
                         rasterized=True)
        axins.set_aspect(1)
        
        encr.open('intersected_ellipses_step_12_from_moist_0000000012_parcels.nc')
        dx = encr.get_box_extent() / encr.get_box_ncells()
        x_pos = encr.get_dataset(11, 'x_position')
        z_pos = encr.get_dataset(11, 'z_position')
        ind = np.argwhere((x_pos >= xlo - dx[0]) & (x_pos <= xhi + dx[0]) &
                          (z_pos >= ylo - dx[1]) & (z_pos <= yhi + dx[1]))
        ind = ind.squeeze()
        hum = encr.get_dataset(11, 'humidity', indices=ind)
        hum = hum * scale
        isort = np.argsort(hum)
        ell = encr.get_ellipses(step=11, indices=ind[isort])
        axins.add_collection(ell)
        ell.set_offset_transform(axins.transData)
        ell.set_clip_box(axins.bbox)
        ell.set_alpha(1.0)
        ell.set_rasterized(True)
        norm = cls.Normalize(vmin=vmin, vmax=vmax)
        ell.set_facecolor(my_cmap(norm(hum[isort])))
        encr.close()  
        
        axins.set_xlabel('')
        axins.set_ylabel('')
        axins.set_rasterized(True)

        axins.set_xlim(xlo, xhi)
        axins.set_ylim(ylo, yhi)

        if not show_axes:
            axins.set_xticks([])
            axins.set_xticklabels([])
            axins.set_yticks([])
            axins.set_yticklabels([])
        else:
            axins.set_xticks(kwargs.pop('xticks', []))
            axins.set_yticks(kwargs.pop('yticks', []))
        ax.indicate_inset_zoom(axins, edgecolor="black", alpha=0.5)
        return axins



    filename = 'epic_unsc_256_lim3r_t6.nc'
    loc = 128
    plane = 'xz'
    step = 0

    # 13 Dec 2022
    # https://stackoverflow.com/questions/67605719/displaying-lowest-values-as-white
    my_cmap = mpl.colors.LinearSegmentedColormap.from_list('', ['white',
                                                                *plt.cm.Blues(np.arange(255))])


    mpl.rcParams['figure.figsize'] = 9, 4
    mpl.rcParams['font.size'] = 12
    legend_dict['ncol'] = 3


    encr = nc_reader()

    dpi = 960

    left = 1250
    right = 5250
    top = 5000
    bottom = 2000

    # 7 Feb 2022
    # https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/figure_size_units.html
    #cm = 1/2.54  # centimeters in inches
    fig = plt.figure(figsize=(9, 4), dpi=dpi) #13*cm, 5*cm), dpi=dpi)

    ax = plt.gca()

    encr.open('intersected_ellipses_step_12_from_moist_0000000012_parcels.nc')
    dx = encr.get_box_extent() / encr.get_box_ncells()
    x_pos = encr.get_dataset(11, 'x_position')
    z_pos = encr.get_dataset(11, 'z_position')
    ind = np.argwhere((x_pos >= left - dx[0]) & (x_pos <= right + dx[0]) &
                      (z_pos >= bottom - dx[1]) & (z_pos <= top + dx[1]))
    ind = ind.squeeze()
    hum = encr.get_dataset(11, 'humidity', indices=ind)

    dmin = hum.min()
    dmax = hum.max()
    
    encr.close()

    ncfile = nc.Dataset(filename, "r", format="NETCDF4")
    data = np.array(ncfile.variables['humidity']).squeeze()
    origin = [0, 0, 0]
    extent = [6280, 6280, 6280]

    dmin = min(dmin, data.min())
    dmax = max(dmax, data.max())

    # scale data
    # note: vmin < dmin, vmax > dmax
    vmin = 0.0
    vmax = 0.09
    scale = vmax / dmax

    data = data * scale
    hum = hum * scale
    
    ax.pcolormesh(np.array(ncfile.variables['X']).squeeze(),                                      
                  np.array(ncfile.variables['Z']).squeeze(),
                  data,
                  vmin=vmin,
                  vmax=vmax,
                  cmap=my_cmap,
                  shading='gouraud',
                  rasterized=True)
    ax.set_aspect(1)

    ax.set_xlim([left, right])
    ax.set_ylim([bottom, top])

    #ax.set_yticks([])
    #ax.set_xticks([])

    ax.set_xlabel(r'$x$ (m)')
    ax.set_ylabel(r'$z$ (m)')

    encr.close()

    encr.open('intersected_ellipses_step_12_from_moist_0000000012_parcels.nc')
    isort = np.argsort(hum)   
    ell = encr.get_ellipses(step=11, indices=ind[isort])
    ax.add_collection(ell)
    ell.set_offset_transform(ax.transData)
    ell.set_clip_box(ax.bbox)
    ell.set_alpha(1.0)
    ell.set_rasterized(True)
    norm = cls.Normalize(vmin=vmin, vmax=vmax)
    ell.set_facecolor(my_cmap(norm(hum[isort])))
    encr.close()

    # corresponds to 32 cells
    axins = do_zoom(
        ncr = encr,
        step = step,
        dset = 'humidity',
        loc=loc,
        plane=plane,
        cmap=my_cmap,
        vmin=vmin,
        vmax=vmax,
        scale=scale,
        ax = ax,
        bounds = [1., 0.47, 0.6, 0.6],
        xlo = 3700,
        xhi = 4485,
        ylo = 3400,
        yhi = 4185,
        fname = filename)

    # sub region of the sub-region image, corresponds to 16 cells
    axins2 = do_zoom(
        ncr = encr,
        step = step,
        dset = 'humidity',
        loc=loc,
        plane=plane,
        cmap=my_cmap,
        vmin=vmin,
        vmax=vmax,
        scale=scale,
        ax = axins,
        bounds = [1.2, 0.06, 0.9, 0.9],
        xlo = 3900,
        xhi = 4292.5,
        ylo = 3600,
        yhi = 3992.5,
        fname = filename)

    # zoom of sub-sub-region, corresponds to 8 cells
    axins3 = do_zoom(
        ncr = encr,
        step = step,
        dset = 'humidity',
        loc=loc,
        plane=plane,
        cmap=my_cmap,
        vmin=vmin,
        vmax=vmax,
        scale=scale,
        ax = axins2,
        #bounds = [1.2, 0.05, 0.9, 0.9],
        bounds = [0.1, -1.04, 0.9, 0.9], 
        xlo = 4000,
        xhi = 4196.25,
        ylo = 3700,
        yhi = 3896.25,
        fname = filename)

    # zoom of sub-sub-sub-region, corresponds to 4 cels
    do_zoom(
        ncr = encr,
        step = step,
        dset = 'humidity',
        loc=loc,
        plane=plane,
        cmap=my_cmap,
        vmin=vmin,
        vmax=vmax,
        scale=scale,
        ax = axins3,
        bounds = [-1.3, 0.06, 0.9, 0.9],
        #bounds = [-1.4, 0.06, 0.9, 0.9],
        xlo = 4050,
        xhi = 4148.125,
        ylo = 3750,
        yhi = 3848.125,
        fname = filename,
        show_axes=True,
        xticks=[4070, 4100, 4130],
        yticks=[3770, 3800, 3830])

    ax.set_rasterized(True)


    plt.savefig('moist_bubble_zooms.pdf',
                bbox_inches='tight', dpi=450) #dpi)

#    os.system('convert moist_bubble_zooms.png eps3:moist_bubble_zooms.eps')

    #os.chdir(cwd)

except Exception as ex:
    print(ex)
