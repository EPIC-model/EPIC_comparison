#python graphical_abstract.py --save-dir figures/ --data-dir epic_vs_ps/Straka
import argparse
import os
from tools.nc_reader import nc_reader
#from tools.pmpic_reader import pmpic_reader
from mpl_toolkits.axes_grid1 import ImageGrid
from tools.mpl_style import *
from tools.units import units
from tools.utils import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import colorcet as cc
import matplotlib.colors as cls

try:
    data_dir = '../data/moist/'
    parcel_intersections = 'intersected_ellipses_step_4_from_moist_0000000004_parcels.nc'
    field_file = 'extracted_steps_from_moist_fields.nc'

    def do_zoom(ncr, step, dset, loc, plane, ax, bounds, xlo, xhi, ylo, yhi, fname, **kwargs):
        ncr.open(fname)

        data = ncr.get_dataset(step, name=dset)

        cmap = kwargs.pop('cmap', 'Blues')
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)
        show_axes = kwargs.pop('show_axes', False)
        scale = kwargs.pop('scale', 1.0)

        # 4 Feb 2022
        # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/zoom_inset_axes.html
        # inset axes....
        axins = ax.inset_axes(bounds, alpha=1.0)

        data = data * scale

        _, _ = make_imshow(ax=axins,
                           plane=plane,
                           loc=loc,
                           fdata=data,
                           step=step,
                           ncr=ncr,
                           cmap=cmap,
                           vmin=vmin,
                           vmax=vmax,
                           interpolation='bilinear',
                           interpolation_stage='rgba',
                           colorbar=False)

        ncr.close()

        encr.open(os.path.join(data_dir, parcel_intersections))
        dx = encr.get_box_extent() / encr.get_box_ncells()
        x_pos = encr.get_dataset(3, 'x_position')
        z_pos = encr.get_dataset(3, 'z_position')
        ind = np.argwhere((x_pos >= xlo - dx[0]) & (x_pos <= xhi + dx[0]) &
                          (z_pos >= ylo - dx[1]) & (z_pos <= yhi + dx[1]))
        ind = ind.squeeze()
        hum = encr.get_dataset(3, 'humidity', indices=ind)
        hum = hum * scale
        isort = np.argsort(hum)
        ell = encr.get_ellipses(step=3, indices=ind[isort])
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


    filename = os.path.join(data_dir, field_file)
    loc = 256
    plane = 'xz'
    step = 0

    # 13 Dec 2022
    # https://stackoverflow.com/a/46778420
    my_cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["white", "darkblue"], 255)

    # 13 Dec 2022
    # https://stackoverflow.com/questions/67605719/displaying-lowest-values-as-white
    my_cmap = mpl.colors.LinearSegmentedColormap.from_list('', ['white',
                                                                *plt.cm.Blues(np.arange(255))])


    mpl.rcParams['figure.figsize'] = 9, 4
    mpl.rcParams['font.size'] = 12
    legend_dict['ncol'] = 3


    encr = nc_reader()

    dpi = 960.0

    left = 1000
    right = 5000
    top = 5000
    bottom = 2000

    # 7 Feb 2022
    # https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/figure_size_units.html
    cm = 1/2.54  # centimeters in inches
    fig = plt.figure(figsize=(9, 4), dpi=dpi)

    ax = plt.gca()

    encr.open(os.path.join(data_dir, parcel_intersections))
    dx = encr.get_box_extent() / encr.get_box_ncells()
    x_pos = encr.get_dataset(3, 'x_position')
    z_pos = encr.get_dataset(3, 'z_position')
    ind = np.argwhere((x_pos >= left - dx[0]) & (x_pos <= right + dx[0]) &
                      (z_pos >= bottom - dx[1]) & (z_pos <= top + dx[1]))
    ind = ind.squeeze()
    hum = encr.get_dataset(3, 'humidity', indices=ind)

    dmin = hum.min()
    dmax = hum.max()

    encr.close()

    encr.open(filename)
    data = encr.get_dataset(step, name='humidity')

    cs = get_plane(plane, loc, data)
    dmin = min(dmin, cs.min())
    dmax = max(dmax, cs.max())

    # vmin < dmin, vmax > dmax
    vmin = 0.0
    vmax = 0.0012
    # scale data such that they are in vmin, vmax
    scale = 1.0 #vmax / dmax
    data = data * scale
    hum = hum * scale

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

    ax.set_xlabel(r'$x$ (m)')
    ax.set_ylabel(r'$z$ (m)')

    encr.close()

    encr.open(os.path.join(data_dir, parcel_intersections))
    isort = np.argsort(hum)
    ell = encr.get_ellipses(step=3, indices=ind[isort])
    ax.add_collection(ell)
    ell.set_offset_transform(ax.transData)
    ell.set_clip_box(ax.bbox)
    ell.set_alpha(1.0)
    ell.set_rasterized(True)
    norm = cls.Normalize(vmin=vmin, vmax=vmax)
    ell.set_facecolor(my_cmap(norm(hum[isort])))
    encr.close()

    # corresponds to 64 cells
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
        xlo = 3200,
        xhi = 3985,
        ylo = 3200,
        yhi = 3985,
        fname = filename)

    # sub region of the sub-region image, corresponds to 32 cells
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
        xlo = 3400,
        xhi = 3792.5,
        ylo = 3400,
        yhi = 3792.5,
        fname = filename)

    # zoom of sub-sub-region, corresponds to 16 cells
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
        xlo = 3500,
        xhi = 3696.25,
        ylo = 3500,
        yhi = 3696.25,
        fname = filename)

    # zoom of sub-sub-sub-region, corresponds to 8 cels
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
        bounds = [-1.4, 0.06, 0.9, 0.9],
        xlo = 3550,
        xhi = 3648.125,
        ylo = 3550,
        yhi = 3648.125,
        show_axes=True,
        xticks=[3570, 3600, 3630],
        yticks=[3570, 3600, 3630],
        fname = filename)

    ax.set_rasterized(True)

    plt.savefig('moist_bubble_zooms.pdf',
                bbox_inches='tight', dpi=450)

    #os.system('convert graphical_abstract.png eps3:graphical_abstract.eps')

    #os.chdir(cwd)

except Exception as ex:
    print(ex)
