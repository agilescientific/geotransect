#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for plotting.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import os

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import hsv_to_rgb
import matplotlib.transforms as transforms
from matplotlib import gridspec
from shapely.ops import transform

from autowrap import on_draw
import utils
from feature_plot import plot_feature_well


def plot_line(m, line, colour='b', lw=1, alpha=1):
    """
    Plots a line given a line with lon,lat coordinates.

    Note:
        This means you probably have to call shapely `transform` on your
        line before passing it to this function.

        There is a helper partial function in utils called `utm2lola` which
        makes this easy.

    Args:
        m (Basemap): A matplotlib Basemap.
        line (shape): A shapely geometry.
        colour (str): A colour from the matplotlib dictionary.

    Returns:
        List: A list of matplotlib lines.
    """
    lo, la = line.xy
    x, y = m(lo, la)
    return m.plot(x, y,
                  color=colour,
                  linewidth=lw,
                  alpha=alpha,
                  solid_capstyle='round')


def plot_point(m, point, colour='b', shape='o', alpha=1, zorder=None):
    """
    Plots a point given a point with lon,lat coordinates.

    Note:
        This means you probably have to call shapely `transform` on your
        point before passing it to this function.

        There is a helper partial function in utils called `utm2lola` which
        makes this easy.

    Args:
        m (Basemap): A matplotlib Basemap.
        point (shape): A shapely geometry.
        colour (str): A colour from the matplotlib dictionary.

    Returns:
        List: A list of matplotlib points.
    """
    lo, la = point.xy
    x, y = m(lo, la)
    return m.scatter(x, y,
                     s=20,
                     color=colour,
                     alpha=alpha,
                     zorder=zorder)


def draw_basemap(m):
    """
    Puts some standard bits of decoration on a matplotlib Basemap.

    Args:
        m (Basemap): A matplotlib Basemap.

    Returns:
        m (Basemap): The newly decorated Basemap.
    """
    m.drawcoastlines(color='#9caf9c')
    m.drawcountries()
    m.fillcontinents(color='#d8e3d8')
    m.drawmapboundary(color='gray')
    m.drawmeridians(np.arange(0, 360, 1), color='gray')
    m.drawparallels(np.arange(-90, 90, 1), color='gray')
    return m


def add_subplot_axes(ax, rect, axisbg='w'):
    """
    Facilitates the addition of a small subplot within another plot.

    From: http://stackoverflow.com/users/2309442/pablo
    License: CC-BY-SA

    Args:
        ax (axis): A matplotlib axis.
        rect (list): A rect specifying [left pos, bottom pos, with, height]

    Returns:
        axis: The sub-axis in the specified position.
    """
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x, y, width, height], axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax


def plot(tc):
    """
    Constructs a multi-subplot matplotlib figure.

    Args:
        transect (TransectContainer): A transect container.
    """
    h = 15
    mw = 16  # width of main section (inches) must be div by 4
    fw = 5   # width of the feature plot (inches) must be div by 5
    n_grids = len(tc.potfield.data)

    fig = plt.figure(figsize=(mw + fw + 1, 15),
                     facecolor='w',
                     edgecolor='k',
                     dpi=None,
                     frameon=True)

    gs = gridspec.GridSpec(h, mw + fw + 1)

    # Left-hand column.
    header = fig.add_subplot(gs[0:1, 0:mw/2])
    description = fig.add_subplot(gs[1:3, 0:mw/2])
    locmap = fig.add_subplot(gs[0:3, mw/2:mw])  # Aspect = 8:3
    elevation = fig.add_subplot(gs[3, :mw])
    xsection = fig.add_subplot(gs[4:h-n_grids, :mw])
    xsec_logs = fig.add_subplot(gs[4:h-n_grids, :mw])
    potfield = fig.add_subplot(gs[h-n_grids:, :mw])

    # Right-hand column.
    log_header = fig.add_subplot(gs[0:1, -1*fw:])
    # log_plot is dealt with by passing gs to feature_plot.plot_feature_well()
    #logo = fig.add_subplot(gs[-1:, -1*fw:])

    # Adjust white space between subplots (maybe put in function?)
    # ------------------------------------------------------------ #
    left = 0.05     # left side of the subplots of the figure
    right = 0.95    # right side of the subplots of the figure
    bottom = 0.05   # bottom of the subplots of the figure
    top = 0.95      # top of the subplots of the figure
    wspace = 0.05   # blank w space between subplots
    hspace = 0.1   # blank h space between subplots

    fig.subplots_adjust(left, bottom, right, top, wspace, hspace)

    bbox = {'fc': 'w', 'pad': 0, 'ec': 'none', 'alpha': 0.5}
    props = {'ha': 'left', 'va': 'center', 'bbox': bbox}

    # Header
    # ---------------------------------------------------------#
    header.axis("off")
    dy = 0.2
    header.text(0.0, 0.5 + dy, tc.title,
                props,
                fontsize=30,
                horizontalalignment='left',
                verticalalignment='bottom'
                )

    # horizontal line
    header.axhline(y=0.5,
                   xmin=0,
                   xmax=1.25,
                   linewidth=2,
                   color='k')

    # Meta description
    header.text(1.0, 0.5 + dy,
                (tc.subtitle),
                fontsize=14,
                horizontalalignment='right',
                verticalalignment='bottom')

    # Wrap text
    fig.canvas.mpl_connect('draw_event', lambda event: on_draw(event))

    # Paragraph description
    # -----------------------------------------------------#
    description.text(0, 1.0,
                     tc.description,
                     horizontalalignment='left',
                     verticalalignment='top',
                     fontsize=12
                     )

    description.axis('off')

    # Map
    # ---------------------------------------------------------#
    tx, ty = tc.data.coords.xy

    res = 'h'   # c, l, i, h, f

    # Generate Basemap object.
    m = Basemap(projection='tmerc',
                lon_0=tc.locmap.mid.x, lat_0=tc.locmap.mid.y,
                resolution=res,
                llcrnrlon=tc.locmap.ll.x, llcrnrlat=tc.locmap.ll.y,
                urcrnrlon=tc.locmap.ur.x, urcrnrlat=tc.locmap.ur.y,
                ax=locmap)

    # Finish drawing the basemap.
    draw_basemap(m)

    for layer, details in tc.locmap.layers.items():
        data = getattr(tc.locmap, layer)
        for l in data:
            line_t = utils.utm2lola(l)
            params = details.get('params', None)
            if params:
                plot_line(m, line_t, **params)
            else:
                plot_line(m, line_t, colour='k', alpha=0.5)
        if layer == "wells":
            for p in data:
                point_t = utils.utm2lola(p)
                params = details.get('params', None)
                if params:
                    plot_point(m,
                               point_t,
                               zorder=100,
                               **params)
                else:
                    plot_point(m,
                               point_t,
                               zorder=100,
                               colour='k',
                               alpha=0.5)
                lo, la = point_t.xy
                x, y = m(lo, la)
                # Plot the names of wells in the xsection
                # When we have names we can also plot special symbol
                # for the feature well.

    # Plot this transect line
    line_t = utils.utm2lola(tc.data)
    plot_line(m, line_t, colour='r', lw=3)

    # Elevation and bedrock plot
    # -----------------------------------------------------------#
    for i, geo in enumerate(tc.bedrock.data[:-1]):
        lim1 = tc.bedrock.coords[i]
        lim2 = tc.bedrock.coords[i + 1]
        idx = np.where(np.logical_and(tc.elevation.coords >= lim1,
                                      tc.elevation.coords <= lim2))[0]
        if len(idx) > 1:
            if idx[-1] < tc.elevation.coords.size - 1:
                idx = np.append(idx, (idx[-1] + 1))
        hsv = np.array([[geo["AV_HUE"], geo["AV_SAT"],
                         geo["AV_VAL"]]]).reshape(1, 1, 3)
        color = hsv_to_rgb(hsv / 255.)

        elevation.bar(tc.elevation.coords[idx],
                      tc.elevation.data[idx],
                      width=1.0,
                      linewidth=0,
                      color=color.flatten())

    max_height = np.amax(tc.elevation.all_data)

    elevation.set_ylim((0, max_height))
    elevation.set_xlim(tc.extents[:2])
    elevation.set_yticks([0, int(max_height),
                          int(np.amax(tc.elevation.data))])
    elevation.set_xticks([])
    elevation.tick_params(axis='y', which='major', labelsize=8)
    elevation.patch.set_alpha(0.1)
    elevation.set_ylabel("Elevation [m]", fontsize=8)
    elevation.grid(True)
    elevation.tick_params(axis='x', which='major', labelsize=0)
    elevation.xaxis.grid(True, which='major')
    elevation.text(0.0, .5*max_height,
                   "Surface geology",
                   props,
                   fontsize=10)
    elevation.set_frame_on(False)

    # Seismic cross section
    # ------------------------------------------------------------#
    for coords, data in zip(tc.seismic.coords, tc.seismic.data):
        xsection.imshow(data,
                        extent=[np.amin(coords) / 1000.0,
                                np.amax(coords) / 1000.0,
                                tc.range[-1], 0],
                        aspect="auto", cmap=tc.seismic_cmap)

    for horizon, data in tc.horizons.data.items():
        print horizon, data.shape, np.nanmin(data), np.nanmax(data)
        coords = tc.horizons.coords[horizon]
        xsection.plot(coords, data, 'o', ms=50)

    plot_axis = [tc.extents[0] / 1000., tc.extents[1] / 1000.,
                 tc.extents[2], tc.extents[3]]
    xsection.axis(plot_axis)
    xsection.set_xticks([])
    xsection.set_ylabel("Depth [m]", fontsize=8)
    xsection.tick_params(axis='y', which='major', labelsize=8)
    xsection.tick_params(axis='x', which='major', labelsize=0)
    xsection.xaxis.grid(True, which='major')
    xsection.grid(True)
    xsection.set_frame_on(False)

    # Potential field data
    # -----------------------------------------------------------#
    for i, (field, payload) in enumerate(tc.potfield.data.items()):
        bot = 1 - (i+1.)/n_grids
        height = (1./n_grids) - 0.05
        rect = [0, bot, 1, height]
        this_ax = add_subplot_axes(potfield, rect)

        this_ax.scatter(payload['coords'],
                        payload['data'],
                        c=payload['colour'],
                        cmap=payload['cmap'],
                        s=1,
                        edgecolor='')

        this_ax.set_xlim(tc.extents[:2])
        this_ax.set_frame_on(False)
        this_ax.set_xticks([])
        this_ax.tick_params(axis='y', which='major', labelsize=8)
        this_ax.grid(True)
        scale = payload['scale']
        if scale:
            this_ax.set_ylim(float(scale[1]), float(scale[0]))
        ymax = this_ax.get_ylim()[1]
        this_ax.text(1000, ymax, field,
                     va='top',
                     fontsize=10,
                     color='k')

    potfield.axis(plot_axis)
    potfield.set_xlim(tc.extents[:2])
    potfield.set_frame_on(False)
    potfield.set_yticks([])
    potfield.xaxis.grid(True, which='major')
    potfield.tick_params(axis='x', which='major', labelsize=10)
    potfield.set_xlabel("Transect range [km]", fontsize=10)

    # Log overlays
    # --------------------------------------------------------#
    if tc.locmap.layers.get('wells'):
        if tc.locmap.layers['wells'].get('colour'):
            well_colour = tc.locmap.layers['wells']['colour']
        else:
            well_colour = tc.settings.get('default_colour')
            if not well_colour:
                well_colour = 'k'

    for name, las, pos in zip(tc.log.names, tc.log.data, tc.log.coords):
        c = tc.seismic_log_colour

        if name == tc.feature_well:
            alpha, lw = 0.5, 1.5
            weight = 'bold'
        else:
            alpha, lw = 0.25, 1.0
            weight = 'normal'

        if las:
            data = np.nan_to_num(las.data[tc.seismic_log])
            data /= np.amax(data)
            lgsc = 0.015  # hack to compress the log width
            data *= lgsc * (tc.extents[1] - tc.extents[0])
            data += pos - 0.5 * np.amax(data)
            z = las.data['DEPT']
            xsec_logs.axvline(x=pos,
                              color=well_colour,
                              alpha=alpha,
                              lw=lw,
                              ymin=0,
                              ymax=z[-1])
            xsec_logs.plot(data, z, c, lw=0.5, alpha=0.75)
            xsec_logs.set_xlim((tc.extents[0], tc.extents[1]))
            xsec_logs.set_ylim((tc.extents[2], tc.extents[3]))
            xsec_logs.axis("off")
        else:
            # Need to get TD from SHP or well header sheet.
            z = [tc.extents[2]-40]
            xsec_logs.axvline(x=pos, color=well_colour, alpha=0.25)

        elevation.axvline(x=pos, color=well_colour, alpha=0.25)
        potfield.axvline(x=pos, color=well_colour, alpha=0.25)

        # Well name annotation
        xsec_logs.text(pos, 20,
                       name,
                       color=well_colour,
                       va='top',
                       ha='center',
                       fontsize=10,
                       weight=weight)

    # Log type annotation, top left
    xsec_logs.text(500, 20,
                   tc.seismic_log+' log',
                   color=c,
                   va='top',
                   fontsize=12)

    #  Feature plot
    # -----------------------------------------------------#
    if tc.feature_well:
        log_header.text(0.0, 1.0,
                        ("Well " + tc.feature_well),
                        verticalalignment='top',
                        horizontalalignment='left',
                        fontsize=14,
                        fontweight='bold'
                        )
        log_header.axis("off")

        plot_feature_well(tc, gs)

    # Logo
    # --------------------------------------------------------------#
    # path = os.path.join(tc.data_dir, tc.settings['images_dir'], 'logo.png')
    # im = Image.open(path)
    # im = np.array(im).astype(np.float) / 255

    # logo.axhline(y=0.7,
    #              xmin=0.1,
    #              xmax=0.9,
    #              linewidth=1,
    #              color='k')
    # logo.imshow(im)
    # logo.axis("off")

    path = os.path.join(tc.data_dir, tc.settings['images_dir'], 'logo.png')
    im = Image.open(path)
    #im = im.thumbnail((100,100), Image.ANTIALIAS)
    w, h = im.size

    # We need a float array between 0-1, rather than
    # a uint8 array between 0-255
    im = np.array(im).astype(np.float) / 255

    # With newer (1.0) versions of matplotlib, you can
    # use the "zorder" kwarg to make the image overlay
    # the plot, rather than hide behind it... (e.g. zorder=10)
    fig.figimage(im, fig.bbox.xmax - w, fig.bbox.ymax - h)

    # Finish
    # --------------------------------------------------------------#
    plt.show()
