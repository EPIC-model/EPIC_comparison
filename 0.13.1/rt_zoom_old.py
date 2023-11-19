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

    def do_zoom(ncr, step, dset, loc, plane, ax, bounds, xlo, xhi, ylo, yhi, fname, **kwargs):
        ncr.open(fname)

        data = ncr.get_dataset(step, name=dset)

        cmap = kwargs.pop('cmap', 'Blues')
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)
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

        encr.open('../mnt/archer2/epic3d-paper-runs/rayleigh-taylor/intersected_ellipses_step_2_from_epic_rt_384_0000000001_parcels.nc')
        dx = encr.get_box_extent() / encr.get_box_ncells()
        x_pos = encr.get_dataset(0, 'x_position')
        z_pos = encr.get_dataset(0, 'z_position')
        ind = np.argwhere((x_pos >= xlo - dx[0]) & (x_pos <= xhi + dx[0]) &
                          (z_pos >= ylo - dx[1]) & (z_pos <= yhi + dx[1]))
        ind = ind.squeeze()
        hum = encr.get_dataset(0, 'buoyancy', indices=ind)
        hum = hum * scale
        isort = np.argsort(hum)
        ell = encr.get_ellipses(step=0, indices=ind[isort])
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
        axins.set_xticks([])
        axins.set_xticklabels([])
        axins.set_yticks([])
        axins.set_yticklabels([])
        ax.indicate_inset_zoom(axins, edgecolor="black", alpha=0.5)
        return axins


    filename = '../mnt/archer2/epic3d-paper-runs/rayleigh-taylor/extracted_steps_from_epic_rt_384_fields.nc'
    loc = 192
    plane = 'xz'
    step = 5

    # 13 Dec 2022
    # https://stackoverflow.com/a/46778420
    my_cmap = cc.cm['coolwarm'] #mpl.colors.LinearSegmentedColormap.from_list("", ["white", "darkblue"], 255)

    # 13 Dec 2022
    # https://stackoverflow.com/questions/67605719/displaying-lowest-values-as-white
#    my_cmap = mpl.colors.LinearSegmentedColormap.from_list('', ['white',
 #                                                               *plt.cm.Blues(np.arange(255))])


    mpl.rcParams['figure.figsize'] = 9, 4
    mpl.rcParams['font.size'] = 12
    legend_dict['ncol'] = 3


    encr = nc_reader()

    dpi = 960.0

    left = -np.pi/2
    right = np.pi/2
    top = np.pi/2
    bottom = -np.pi/2

    # 7 Feb 2022
    # https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/figure_size_units.html
    cm = 1/2.54  # centimeters in inches
    fig = plt.figure(figsize=(13*cm, 5*cm), dpi=dpi)

    ax = plt.gca()

    encr.open('../mnt/archer2/epic3d-paper-runs/rayleigh-taylor/intersected_ellipses_step_2_from_epic_rt_384_0000000001_parcels.nc')
    dx = encr.get_box_extent() / encr.get_box_ncells()
    x_pos = encr.get_dataset(0, 'x_position')
    z_pos = encr.get_dataset(0, 'z_position')
    ind = np.argwhere((x_pos >= left - dx[0]) & (x_pos <= right + dx[0]) &
                      (z_pos >= bottom - dx[1]) & (z_pos <= top + dx[1]))
    ind = ind.squeeze()
    hum = encr.get_dataset(0, 'buoyancy', indices=ind)

    dmin = hum.min()
    dmax = hum.max()

    encr.close()

    encr.open(filename)
    data = encr.get_dataset(step, name='buoyancy')

    cs = get_plane(plane, loc, data)
    dmin = min(dmin, cs.min())
    dmax = max(dmax, cs.max())

    # vmin < dmin, vmax > dmax
    vmin = -1.0
    vmax = 1.0
    # scale data such that they are in vmin, vmax
    scale = vmax / dmax
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

    ax.set_yticks([])
    ax.set_xticks([])

    ax.set_xlabel('')
    ax.set_ylabel('')

    encr.close()

    encr.open('../mnt/archer2/epic3d-paper-runs/rayleigh-taylor/intersected_ellipses_step_2_from_epic_rt_384_0000000001_parcels.nc')
    isort = np.argsort(hum)
    ell = encr.get_ellipses(step=0, indices=ind[isort])
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
        dset = 'buoyancy',
        loc=loc,
        plane=plane,
        cmap=my_cmap,
        vmin=vmin,
        vmax=vmax,
        scale=scale,
        ax = ax,
        bounds = [1., 0.47, 0.6, 0.6],
        xlo = -np.pi/4.0,
        xhi =  np.pi/4.0,
        ylo = -np.pi/4.0,
        yhi =  np.pi/4.0,
        fname = filename)

    # sub region of the sub-region image, corresponds to 16 cells
    axins2 = do_zoom(
        ncr = encr,
        step = step,
        dset = 'buoyancy',
        loc=loc,
        plane=plane,
        cmap=my_cmap,
        vmin=vmin,
        vmax=vmax,
        scale=scale,
        ax = axins,
        bounds = [1.2, 0.06, 0.9, 0.9],
        xlo = -np.pi/8.0,
        xhi =  np.pi/8.0,
        ylo = -np.pi/8.0,
        yhi =  np.pi/8.0,
        fname = filename)

    # zoom of sub-sub-region, corresponds to 8 cells
    axins3 = do_zoom(
        ncr = encr,
        step = step,
        dset = 'buoyancy',
        loc=loc,
        plane=plane,
        cmap=my_cmap,
        vmin=vmin,
        vmax=vmax,
        scale=scale,
        ax = axins2,
        #bounds = [1.2, 0.05, 0.9, 0.9],
        bounds = [0.1, -1.04, 0.9, 0.9],
        xlo = -np.pi/16.0,
        xhi =  np.pi/16.0,
        ylo = -np.pi/16.0,
        yhi =  np.pi/16.0,
        fname = filename)

    # zoom of sub-sub-sub-region, corresponds to 4 cels
    do_zoom(
        ncr = encr,
        step = step,
        dset = 'buoyancy',
        loc=loc,
        plane=plane,
        cmap=my_cmap,
        vmin=vmin,
        vmax=vmax,
        scale=scale,
        ax = axins3,
        bounds = [-1.4, 0.06, 0.9, 0.9],
        xlo = -np.pi/32.0,
        xhi =  np.pi/32.0,
        ylo = -np.pi/32.0,
        yhi =  np.pi/32.0,
        fname = filename)

    ax.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)


    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_rasterized(True)

    plt.savefig('rt_zoom.png',
                bbox_inches='tight', dpi=dpi)

    os.system('convert rt_zoom.png eps3:rt_zoom.eps')

    #os.chdir(cwd)

except Exception as ex:
    print(ex)
