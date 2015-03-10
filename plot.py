#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for plotting.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import csv
import os
from functools import partial

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import hsv_to_rgb
import matplotlib.transforms as transforms
from matplotlib import gridspec
from shapely.ops import transform
import pyproj as pp

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


def plot_point(m, point, colour='b', shape='o', alpha=1):
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
                     s=10,
                     color=colour,
                     alpha=alpha)


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
    grids = 2

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
    xsection = fig.add_subplot(gs[4:h-grids, :mw])
    potfield = fig.add_subplot(gs[h-grids:, :mw])

    # Right-hand column.
    log_header = fig.add_subplot(gs[0:1, -1*fw:])
    log = fig.add_subplot(gs[6:h-1, :mw])
    logo = fig.add_subplot(gs[-1:, -1*fw:])

    # ------------------------------------------------------------ #
    # Adjust white space between subplots (maybe put in function?)
    # ------------------------------------------------------------ #
    left = 0.05     # left side of the subplots of the figure
    right = 0.95    # right side of the subplots of the figure
    bottom = 0.1   # bottom of the subplots of the figure
    top = 0.9      # top of the subplots of the figure
    wspace = 0.1   # width reserved for blank space between subplots
    hspace = 0.05   # height reserved for white space between subplots

    fig.subplots_adjust(left, bottom, right, top, wspace, hspace)

    bbox = {'fc': 'w', 'pad': 0, 'ec': 'none', 'alpha': 0.5}
    props = {'ha': 'left', 'va': 'center', 'bbox': bbox}

    # ---------------------------------------------------------#
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
                (tc.basin),
                fontsize=14,
                horizontalalignment='right',
                verticalalignment='bottom')

    # -----------------------------------------------------#
    # Paragraph description
    # -----------------------------------------------------#
    description.text(0, 1.0,
                     tc.description,
                     horizontalalignment='left',
                     verticalalignment='top',
                     fontsize=12
                     )

    description.axis('off')

    # --------------------------------------------------------#
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

    # Plot the seismic lines
    for l in tc.locmap.seismic:
        line_t = utils.utm2lola(l)
        plot_line(m, line_t, colour='black', alpha=0.3)

    # Plot the wells
    for p in tc.locmap.wells:
        point_t = utils.utm2lola(p)
        plot_point(m, point_t, colour='black', alpha=0.6)

    # Plot this transect line
    line_t = utils.utm2lola(tc.data)
    plot_line(m, line_t, colour='r', lw=3)

    # Finish drawing the basemap.
    draw_basemap(m)

    # --------------------------------------------------------#
    # Small scale map
    # ---------------------------------------------------------#
    # small_scale.plot(tx, ty)
    # small_scale.patch.set_facecolor("0.95")
    # small_scale.patch.set_edgecolor("0.90")
    # small_scale.set_xticks([])
    # small_scale.set_yticks([])

    # ----------------------------------------------------------#
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
        elevation.plot(tc.elevation.coords[idx], tc.elevation.data[idx],
                       linewidth=3, color=color.flatten())

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
    elevation.text(0.0, .5 * (max_height),
                   "Surface geology",
                   props,
                   fontsize=8,
                   rotation=0)
    elevation.set_frame_on(False)

    # ------------------------------------------------------------#
    # Seismic cross section
    # ------------------------------------------------------------#
    for coords, data in zip(tc.seismic.coords, tc.seismic.data):

        z0 = 0
        depth = 2500
        xsection.imshow(data,
                        extent=[np.amin(coords) / 1000.0,
                                np.amax(coords) / 1000.0,
                                depth, z0],
                        aspect="auto", cmap="Greys")

        plot_axis = [tc.extents[0] / 1000., tc.extents[1] / 1000.,
                     tc.extents[2], tc.extents[3]]
        xsection.axis(plot_axis)

    xsection.set_ylabel("Depth [m]", fontsize=8)
    xsection.set_xlabel("Transect range [km]", fontsize=8)

    xsection.tick_params(axis='y', which='major', labelsize=8)
    xsection.tick_params(axis='y', which='minor', labelsize=8)

    xsection.tick_params(axis='x', which='major', labelsize=8)
    xsection.tick_params(axis='x', which='minor', labelsize=8)

    xsection.grid(True)
    xsection.set_frame_on(False)

    # --------------------------------------------------------#
    # Log overlays
    # --------------------------------------------------------#
    for las, pos in zip(tc.log.data, tc.log.coords):

        data = np.nan_to_num(las.data["GR"])

        # normalize
        data /= np.amax(data)

        # scale position
        lgsc = 0.015  # hack to compress the log width
        data *= lgsc * (tc.extents[1] - tc.extents[0])
        data += pos

        log.plot(data, las.data['DEPT'],
                 'g', lw=0.5, alpha=0.5)

        log.set_xlim((tc.extents[0], tc.extents[1]))
        log.set_ylim((tc.extents[2], tc.extents[3]))

        log.axis("off")

    log.axis('off')

    # -----------------------------------------------------------#
    # Potential field data
    # -----------------------------------------------------------#
    for field in tc.potfield.data:
        potfield.plot(tc.potfield.coords[field], tc.potfield.data[field])
        potfield.set_xlim(tc.extents[:2])
        potfield.set_frame_on(False)
        potfield.set_xticks([])
        potfield.tick_params(axis='y', which='major', labelsize=8)


    # def plot(self, ax):

    #     ax.plot(self.coords, self.data)
    #     ax.set_yticks([-4, 0, 4])
    #     ax.set_xticks([])

    #     ax.ylabel("Anomaly [mGal]", fontsize=8)
    #     ax.tick_params(axis='y', which="major", labelsize=8)

    #     ax.grid(True)
    #     ax.xlim((self.extents[0], self.extents[1]))
    #     ax.ylim((-4, 4))

    #     return ax


    # -----------------------------------------------------#
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

        plot_feature_well(log)

    # --------------------------------------------------------------#
    # Logo
    # --------------------------------------------------------------#
    path = os.path.join(tc.data_dir, tc.settings['images_dir'], 'logo.png')
    im = Image.open(path)
    im = np.array(im).astype(np.float) / 255
    # width, height = im.size

    # Place in axis:
    logo.imshow(im)
    logo.axis('off')

    # Place on figure:
    # fig.figimage(im, 0, fig.bbox.ymax - height)

    # logo.text(0.1, 0.7,
    #           ("Department of Energy \n" +
    #            "Nova Scotia, Canada"),
    #           verticalalignment='top',
    #           horizontalalignment='left')

    # Horizontal line.
    logo.axhline(y=0.7,
                 xmin=0.1,
                 xmax=0.9,
                 linewidth=1,
                 color='k')
    logo.axis("off")

    # --------------------------------------------------------------#
    # Finish
    # --------------------------------------------------------------#

    # Wrap text
    fig.canvas.mpl_connect('draw_event', on_draw)

    plt.show()
