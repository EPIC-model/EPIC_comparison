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

    field_file = 'extracted_steps_from_epic_rt_384_fields.nc'
    field_step = 5
    parcel_file = 'intersected_ellipses_step_2_from_epic_rt_384_0000000001_parcels.nc'
    parcel_step = 0

    def do_zoom(ncr, step, dset, loc, plane, ax, bounds, xlo, xhi, ylo, yhi, fname, **kwargs):
        
        # inset axes....
        axins = ax.inset_axes(bounds, alpha=1.0)

        show_axes = kwargs.pop('show_axes', False)
        scale = kwargs.pop('scale', 1.0)

        data = ncr.get_dataset(step=step, name='buoyancy')
        _, _ = make_imshow(ax=axins,
                           plane=plane,
                           loc=loc,
                           fdata=data,
                           step=step,
                           ncr=ncr,
                           cmap=my_cmap,
                           vmin=vmin,
                           vmax=vmax,
                           interpolation='bilinear',
                           interpolation_stage='rgba',
                           colorbar=False)

        encr.open(parcel_file)
        dx = encr.get_box_extent() / encr.get_box_ncells()
        x_pos = encr.get_dataset(parcel_step, 'x_position')
        z_pos = encr.get_dataset(parcel_step, 'z_position')
        ind = np.argwhere((x_pos >= xlo - dx[0]) & (x_pos <= xhi + dx[0]) &
                          (z_pos >= ylo - dx[1]) & (z_pos <= yhi + dx[1]))
        ind = ind.squeeze()
        buo = encr.get_dataset(0, 'buoyancy', indices=ind)
        #buo = buo * scale
        isort = np.argsort(buo)
        ell = encr.get_ellipses(step=0, indices=ind[isort])
        axins.add_collection(ell)
        ell.set_offset_transform(axins.transData)
        ell.set_clip_box(axins.bbox)
        ell.set_alpha(1.0)
        ell.set_rasterized(True)
        norm = cls.Normalize(vmin=vmin, vmax=vmax)
        ell.set_facecolor(my_cmap(norm(buo[isort])))
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
            axins.set_xticks(kwargs.pop('xticks', []), kwargs.pop('xticklabs', []))
            axins.set_yticks(kwargs.pop('yticks', []), kwargs.pop('yticklabs', []))
            if kwargs.pop('ytick_right', False):
               axins.yaxis.tick_right()
            if kwargs.pop('xtick_top', False):
                axins.xaxis.tick_top()
        ax.indicate_inset_zoom(axins, edgecolor="black", alpha=0.5)
        return axins



    filename = field_file
    loc = 192
    plane = 'xz'
    step = field_step

    my_cmap = cc.cm['coolwarm']
    

    mpl.rcParams['figure.figsize'] = 9, 4
    mpl.rcParams['font.size'] = 12
    legend_dict['ncol'] = 3


    encr = nc_reader()

    dpi = 960

    left = -np.pi/2.0
    right = np.pi/2.0
    top = np.pi/2.0
    bottom = -np.pi/2.0

    # 7 Feb 2022
    # https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/figure_size_units.html
    #cm = 1/2.54  # centimeters in inches
    fig = plt.figure(figsize=(9, 4), dpi=dpi) #13*cm, 5*cm), dpi=dpi)

    ax = plt.gca()

    encr.open(parcel_file)

    tp = encr.get_dataset(step=parcel_step, name='t')

    print("Parcel time:", tp)
    
    dx = encr.get_box_extent() / encr.get_box_ncells()
    x_pos = encr.get_dataset(parcel_step, 'x_position')
    z_pos = encr.get_dataset(parcel_step, 'z_position')
    ind = np.argwhere((x_pos >= left - dx[0]) & (x_pos <= right + dx[0]) &
                      (z_pos >= bottom - dx[1]) & (z_pos <= top + dx[1]))
    ind = ind.squeeze()
    buo = encr.get_dataset(parcel_step, 'buoyancy', indices=ind)

    dmin = buo.min()
    dmax = buo.max()
    
    encr.close()

    ncr = nc_reader()
    ncr.open(filename)

    tf = ncr.get_all('t')

    print("Field time:", tf[field_step])
    
    data = ncr.get_dataset(step=field_step, name='buoyancy')
    origin = ncr.get_box_origin()
    extent = ncr.get_box_extent()

    dmin = min(dmin, data.min())
    dmax = max(dmax, data.max())

    # scale data
    # note: vmin < dmin, vmax > dmax
    vmin = -1.0
    vmax =  1.0
    scale = 1.0
    #scale = vmax / dmax

    #data = data * scale
    #buo = buo * scale

    im, _ = make_imshow(ax=ax,
                        plane=plane,
                        loc=loc,
                        fdata=data,
                        vmin=vmin,
                        vmax=vmax,
                        step=step,
                        ncr=encr,
                        cmap=my_cmap,
                        interpolation='bilinear',
                        interpolation_stage='rgba',
                        colorbar=False)
    
    ax.set_xlim([left, right])
    ax.set_ylim([bottom, top])

    xticks = [-np.pi/2.0, -np.pi/4.0, 0.0, np.pi/4.0, np.pi/2.0]
    xticklabs = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']

    yticks = [-np.pi/2.0, -np.pi/4.0, 0.0, np.pi/4.0, np.pi/2.0]
    yticklabs = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']
    

    ax.set_xticks(xticks, xticklabs)
    ax.set_yticks(yticks, yticklabs)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$z$')

   # encr.close()

    encr.open(parcel_file)
    isort = np.argsort(buo)   
    ell = encr.get_ellipses(step=parcel_step, indices=ind[isort])
    ax.add_collection(ell)
    ell.set_offset_transform(ax.transData)
    ell.set_clip_box(ax.bbox)
    ell.set_alpha(1.0)
    ell.set_rasterized(True)
    norm = cls.Normalize(vmin=vmin, vmax=vmax)
    ell.set_facecolor(my_cmap(norm(buo[isort])))
    encr.close()

    # corresponds to 192 cells
    axins = do_zoom(
        ncr = ncr,
        step = step,
        dset = 'buoyancy',
        loc=loc,
        plane=plane,
        cmap=my_cmap,
        vmin=vmin,
        vmax=vmax,
        scale=scale,
        ax = ax,
        bounds = [1.05, 0.47, 0.6, 0.6],
        xlo = -np.pi/4.0,
        xhi =  np.pi/4.0,
        ylo = -np.pi/4.0,
        yhi =  np.pi/4.0,
        #show_axes=True,
        #xticks=[-np.pi/4.0, np.pi/4.0],
        #yticks=[-np.pi/4.0, np.pi/4.0],
        #xticklabs = [r'$-\pi/4$', r'$\pi/4$'],
        #yticklabs = [r'$-\pi/4$', r'$\pi/4$'],
        #ytick_right=True,
        #xtick_top=True,
        fname = filename)

    # sub region of the sub-region image, corresponds to 96 cells
    axins2 = do_zoom(
        ncr = ncr,
        step = step,
        dset = 'buoyancy',
        loc=loc,
        plane=plane,
        cmap=my_cmap,
        vmin=vmin,
        vmax=vmax,
        scale=scale,
        ax = axins,
        bounds = [1.1, 0.06, 0.9, 0.9],
        xlo = -np.pi/8.0,
        xhi =  np.pi/8.0,
        ylo = -np.pi/8.0,
        yhi =  np.pi/8.0,
        #show_axes=True,
        #xticks=[-np.pi/8.0, np.pi/8.0],
        #yticks=[-np.pi/8.0, np.pi/8.0],
        #xticklabs = [r'$-\pi/8$', r'$\pi/8$'],
        #yticklabs = [r'$-\pi/8$', r'$\pi/8$'],
        #ytick_right=True,
        #xtick_top=True,
        fname = filename)

    # zoom of sub-sub-region, corresponds to 48 cells
    axins3 = do_zoom(
        ncr = ncr,
        step = step,
        dset = 'buoyancy',
        loc=loc,
        plane=plane,
        cmap=my_cmap,
        vmin=vmin,
        vmax=vmax,
        scale=scale,
        ax = axins2,
        bounds = [0.05, -1.0, 0.9, 0.9], 
        xlo = -np.pi/16.0,
        xhi =  np.pi/16.0,
        ylo = -np.pi/16.0,
        yhi =  np.pi/16.0,
        #show_axes=True,
        #xticks=[-np.pi/16.0, np.pi/16.0],
        #yticks=[-np.pi/16.0, np.pi/16.0],
        #xticklabs = [r'$-\pi/16$', r'$\pi/16$'],
        #yticklabs = [r'$-\pi/16$', r'$\pi/16$'],
        #ytick_right=True,
        fname = filename)

    # zoom of sub-sub-sub-region, corresponds to 24 cels
    do_zoom(
        ncr = ncr,
        step = step,
        dset = 'buoyancy',
        loc=loc,
        plane=plane,
        cmap=my_cmap,
        vmin=vmin,
        vmax=vmax,
        scale=scale,
        ax = axins3,
        bounds = [-1.0, 0.0, 0.9, 0.9],
        xlo = -np.pi/32.0,
        xhi =  np.pi/32.0,
        ylo = -np.pi/32.0,
        yhi =  np.pi/32.0,
        fname = filename,
        show_axes=True,
        xticks=[-np.pi/32, 0.0, np.pi/32],
        yticks=[-np.pi/32, 0.0, np.pi/32],
        xticklabs = [r'$-\pi/32$', r'$0$', r'$\pi/32$'],
        yticklabs = [r'$-\pi/32$', r'$0$', r'$\pi/32$']
    )

    ax.set_rasterized(True)

    ncr.close()

    plt.savefig('rt_zooms.pdf',
                bbox_inches='tight', dpi=450) #dpi)

#    os.system('convert moist_bubble_zooms.png eps3:moist_bubble_zooms.eps')

    #os.chdir(cwd)

except Exception as ex:
    print(ex)
