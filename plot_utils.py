#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for plotting.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
import matplotlib.pyplot as plt


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
        list: A list of matplotlib lines.
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
        list: A list of matplotlib points.
    """
    lo, la = point.xy
    x, y = m(lo, la)
    return m.scatter(x, y,
                     s=20,
                     color=colour,
                     alpha=alpha,
                     zorder=zorder)


def draw_basemap(m, tc):
    """
    Puts some standard bits of decoration on a matplotlib Basemap.

    Args:
        m (Basemap): A matplotlib Basemap.
        tc (TransectContainer): We need the map edges from the
            active container. Ideally we'd get them out of the 
            Basemap object, but I can't see how to do this.

    Returns:
        Basemap: The newly decorated Basemap.
    """
    m.drawmapboundary(fill_color='#ccddee')
    m.drawcoastlines(color='#9caf9c')
    m.drawcountries()
    m.fillcontinents(color='#d8e3d8', lake_color='#ccddee')
    m.drawmeridians(np.arange(0, 360, 0.5),
                    color='gray',
                    labels=[False, False, True, False])
    m.drawparallels(np.arange(-90, 90, 0.5),
                    color='gray',
                    labels=[True, False, False, False])
    m.drawmapscale(tc.locmap.ll.x + 0.2, tc.locmap.ll.y + 0.075,
                   tc.locmap.ll.x, tc.locmap.ll.y,
                   20.,
                   barstyle='fancy', labelstyle='simple',
                   fillcolor1='w', fillcolor2='#555555',
                   fontcolor='#555555',
                   zorder=5)

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
