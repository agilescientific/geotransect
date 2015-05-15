#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for plotting.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import os
import datetime

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import hsv_to_rgb
import matplotlib.transforms as transforms
from matplotlib import gridspec

from autowrap import on_draw
import utils
from feature_plot import plot_feature_well
from plot_utils import *
from notice import Notice


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

    # We will save the same figure we make, to ensure the saved figure
    # has everything in the right places. For example, see this discussion:
    # http://stackoverflow.com/questions/7906365/
    save_dpi = getattr(tc, 'save_dpi', tc.settings.get('default_dpi'))
    dpi = save_dpi or 80
    fig = plt.figure(figsize=(mw + fw + 1, 15),
                     facecolor='w',
                     edgecolor='k',
                     dpi=dpi,
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

    # Adjust white space between subplots
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
    print "Header"
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
                   linewidth=1.5,
                   color='k')

    # Subtitle
    header.text(1.0, 0.5 + dy,
                (tc.subtitle),
                fontsize=15,
                horizontalalignment='right',
                verticalalignment='bottom', weight='bold')

    descr = tc.description
    if tc.meta:
        description.text(0, 1.0,
                         tc.domain.upper(),
                         horizontalalignment='left',
                         verticalalignment='bottom',
                         fontsize=14
                         )

        description.text(1.0, 1.0,
                         tc.velocity,
                         horizontalalignment='right',
                         verticalalignment='bottom',
                         fontsize=12
                         )
        descr_pos = 0.8  # Where to position the rest.
    else:
        descr_pos = 1.0

    # Paragraph description
    # -----------------------------------------------------#
    description.text(0, descr_pos,
                     descr,
                     horizontalalignment='left',
                     verticalalignment='top',
                     fontsize=12, family='serif'
                     )

    description.axis('off')

    # Wrap text
    fig.canvas.mpl_connect('draw_event', lambda event: on_draw(event, description))

    # Map
    # ---------------------------------------------------------#
    print "Locmap"
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
    draw_basemap(m, tc)

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

    # Adjust border thickness
    [i.set_linewidth(8) for i in locmap.spines.itervalues()]
    [i.set_color("white") for i in locmap.spines.itervalues()]

    # Elevation and bedrock plot
    # -----------------------------------------------------------#
    print "Elevation"
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

        elevation.plot(tc.elevation.coords[idx],
                       tc.elevation.data[idx],
                       lw=0.5, color='k')

    max_height = np.amax(tc.elevation.all_data)

    elevation.set_ylim((0, max_height))
    elevation.set_xlim(tc.extents[:2])
    print "Elevation X-limits, ", tc.extents[:2]
    elevation.set_yticks([0, int(max_height),
                          int(np.amax(tc.elevation.data))])
    # elevation.set_xticklabels([])
    elevation.tick_params(axis='y', which='major', labelsize=8)
    elevation.patch.set_alpha(0.1)
    elevation.set_ylabel("Elevation [m]", fontsize=8)
    elevation.grid(True)
    elevation.tick_params(axis='x', which='major', labelsize=0)
    elevation.xaxis.grid(True, which='major')
    elevation.text(0.01, 0.8,
                   "Elevation",
                   props, va='bottom',
                   fontsize=11, weight='bold',
                   transform=elevation.transAxes)
    elevation.text(0.01, 0.75,
                   "Surface geology",
                   props, va='top',
                   fontsize=10,
                   transform=elevation.transAxes)
    elevation.set_frame_on(False)
    for tick in elevation.get_xaxis().get_major_ticks():
        tick.set_pad(-8.)
        tick.label1 = tick._get_text1()

    # Seismic cross section
    # ------------------------------------------------------------#
    print "Seismic"
    print "seismic trace coords", tc.seismic.coords
    for coords, data in zip(tc.seismic.coords, tc.seismic.data):

        max_z = data["traces"].shape[0] * data["sample_interval"]
        im = xsection.imshow(data["traces"],
                             extent=[np.amin(coords),  # / 1000.0,
                                     np.amax(coords),  # / 1000.0,
                                     max_z, 0],
                             aspect="auto", cmap=tc.seismic_cmap)

    # Horizons
    colours = ['b', 'g', 'orange', 'c', 'magenta', 'pink']
    for i, (horizon, data) in enumerate(tc.horizons.data.items()):

        coords = tc.horizons.coords[horizon]
        print "Horizon COORDS", coords
        xsection.scatter(coords,  # /1000.0,
                         data,
                         color=colours[i],
                         marker='.')

        # labels
        xsection.text(0.01, 0.025*(i+1),
                      horizon,
                      transform=xsection.transAxes,
                      ha='left', color=colours[i],
                      va='center', fontsize=12)

    # Axes etc.
    plot_axis = [tc.extents[0], tc.extents[1],
                 tc.extents[2], tc.extents[3]]
    xsection.axis(plot_axis)
    xsection.set_xticklabels([])
    if tc.domain.lower() in ['depth', 'd', 'z']:
        xsection.set_ylabel("Depth [m]", fontsize=8)
    else:
        xsection.set_ylabel("TWTT [ms]", fontsize=8)
    xsection.tick_params(axis='y', which='major', labelsize=8)
    xsection.tick_params(axis='x', which='major', labelsize=0)
    xsection.yaxis.grid(True, which='major')
    xsection.xaxis.grid(True, which='major')
    xsection.set_frame_on(False)

    # Seismic colorbar
    colorbar_ax = add_subplot_axes(xsection, [0.975, 0.025, 0.01, 0.1])
    fig.colorbar(im, cax=colorbar_ax)
    colorbar_ax.text(0.5, 0.9, "+",
                     transform=colorbar_ax.transAxes,
                     ha='center', color='white',
                     va='center', fontsize=12)
    colorbar_ax.text(0.5, 0.15, "-",
                     transform=colorbar_ax.transAxes, color='k',
                     ha='center', va='center', fontsize=16)
    colorbar_ax.set_axis_off()

    # Title
    xsection.text(0.01, 0.99,
                  "Seismic",
                  color='k',
                  ha='left', va='top',
                  fontsize=12, weight='bold',
                  transform=xsec_logs.transAxes)

    # Potential field data
    # -----------------------------------------------------------#
    print "Potfields"
    for i, (field, payload) in enumerate(tc.potfield.data.items()):
        bot = 1 - (i + 1.) / n_grids
        height = (1. / n_grids) - 0.05
        rect = [0, bot, 1, height]
        this_ax = add_subplot_axes(potfield, rect)
        sc = this_ax.scatter(payload['coords'],
                             payload['data'],
                             c=payload['colour'],
                             cmap=payload['cmap'],
                             s=1,
                             edgecolor='',
                             vmin=-50, vmax=150)

        this_ax.set_xlim(tc.extents[:2])
        this_ax.set_frame_on(False)
        this_ax.set_xticks([])
        this_ax.tick_params(axis='y', which='major', labelsize=8)
        this_ax.grid(True)
        scale = payload['scale']
        if scale:
            this_ax.set_ylim(float(scale[1]), float(scale[0]))

        if payload['colour_is_file']:
            tcol = '#555555'
        else:
            tcol = payload['colour']

        this_ax.text(0.01, 0.01, field,
                     ha='left', va='bottom',
                     fontsize=10, color=tcol,
                     transform=this_ax.transAxes)

        # potfield colorbar
        # TODO: This doesn't work.
        if payload.get('cmap'):
            # pf_cbar_ax = add_subplot_axes(this_ax, [0.975, 0.1, 0.01, 0.8])
            # fig.colorbar(sc, cax=pf_cbar_ax)
            # pf_cbar_ax.set_axis_off()
            pass

    potfield.axis(plot_axis)
    potfield.set_frame_on(False)
    potfield.set_yticks([])
    potfield.xaxis.grid(True, which='major')
    potfield.tick_params(axis='x', which='major', labelsize=10)

    # Display x-ticklabels in km instead of metres
    labels = potfield.get_xticks().tolist()
    xlabels = [int(label) / 1000 for label in labels]
    potfield.set_xticklabels(xlabels)

    potfield.set_xlabel("Transect range [km]", fontsize=10, ha='center')

    potfield.text(0.01, 1.0, "Potential fields",
                  ha='left', va='top',
                  fontsize=11, weight='bold', color='k',
                  transform=potfield.transAxes)

    # Log overlays
    # --------------------------------------------------------#
    print "Logs"
    if tc.locmap.layers.get('wells'):
        if tc.locmap.layers['wells'].get('colour'):
            well_colour = tc.locmap.layers['wells']['colour']
        else:
            well_colour = tc.settings.get('default_colour')
            if not well_colour:
                well_colour = 'k'

    c = tc.seismic_log_colour
    for name, las, pos in zip(tc.log.names, tc.log.data, tc.log.coords):

        if name == tc.feature_well:
            alpha, lw = 0.5, 1.5
            weight = 'bold'
        else:
            alpha, lw = 0.25, 1.0
            weight = 'normal'

        if name == "P-129":  # This is a hack because we only have one file with tops
            tops = utils.get_tops(tc.tops_file)
            print tops
            # plot tops on the well track:
            for mkr, md in tops.iteritems():
                # Get domain conversion working
                # if tc.domain.lower() in ['time', 'twt', 'twtt', 't']:
                #    md = tc.seismic.velocity.depth2time(depth, pos, dz=depth, dt=dt)
                #    print 'after time conv:', mkr + ': ' + md
                xsec_logs.scatter(x=pos,
                                  y=md,
                                  marker='_',
                                  s=50,
                                  color='b',
                                  alpha=1.0,
                                  zorder=2000)

                # draw text box at the right edge of the last track
                pad = 300.0  # pushes the labels slightly off to the right (map units)
                xsec_logs.text(x=pos + pad,
                               y=md,
                               s=mkr,  # add precluding space for asthetics
                               alpha=0.75,
                               color='k',
                               fontsize='8',
                               horizontalalignment='left',
                               verticalalignment='center',
                               zorder=1000,
                               bbox=dict(facecolor='white',
                                         edgecolor='k',
                                         alpha=0.5,
                                         lw=0.25),
                               weight='light')

        if las:
            data = np.nan_to_num(las.data[tc.seismic_log])
            data /= np.amax(data)
            z = las.data['DEPT']

            if tc.domain.lower() in ['time', 'twt', 'twtt', 't']:
                dt = 0.001
                data = tc.seismic.velocity.depth2time(data, pos, dz=z, dt=dt)
                start = tc.seismic.velocity.depth2timept(las.start, pos)
                z = np.arange(0, len(data), 1) + 1000.0 * start  # ms

            # Some post-processing for display
            lgsc = 0.015  # hack to compress the log width
            data *= lgsc * (tc.extents[1] - tc.extents[0])
            if data.size > 0:
                data += pos - 0.5 * np.amax(data)

            xsec_logs.axvline(x=pos,
                              color=well_colour,
                              alpha=alpha,
                              lw=lw,
                              ymin=0,
                              ymax=z[-1])
            xsec_logs.plot(data, z, c, lw=0.5, alpha=0.75)
        else:
            # Need to get TD from SHP or well header sheet.
            z = [tc.extents[2] - 40]
            xsec_logs.axvline(x=pos, color=well_colour, alpha=0.25)

        elevation.axvline(x=pos, color=well_colour, alpha=alpha, lw=lw)
        potfield.axvline(x=pos,  #/1000.,
                         color=well_colour, alpha=alpha, lw=lw)

        # Well name annotation
        elevation.text(pos, max_height-10,
                       name,
                       color=well_colour,
                       va='top',
                       ha='center',
                       fontsize=10,
                       weight=weight)

    xsec_logs.set_xticks([])

    # Log type annotation, top left
    xsec_logs.text(0.01, 0.965,
                   tc.seismic_log + ' log',
                   color=c,
                   ha='left', va='top',
                   fontsize=12,
                   transform=xsec_logs.transAxes)

    #  Feature plot
    # -----------------------------------------------------#
    if tc.feature_well:
        print "Feature well:", tc.feature_well
        log_header.text(0.0, 1.0,
                        ("Well " + tc.feature_well),
                        verticalalignment='top',
                        horizontalalignment='left',
                        fontsize=14,
                        fontweight='bold'
                        )

        # horizontal line
        log_header.axhline(y=0.5,
                           xmin=0,
                           xmax=1.25,
                           linewidth=1,
                           color='k')

        plot_feature_well(tc, gs)
    else:
        Notice.warning("No feature well")

    log_header.axis("off")

    # Logo etc.
    # --------------------------------------------------------------#
    print "Logo"

    try:
        path = os.path.join(tc.data_dir, tc.settings['logo_file'])
        im = Image.open(path)
        im.thumbnail((fig.dpi, fig.dpi), Image.ANTIALIAS)
    except IOError:
        print "Image is missing", path
        im = np.zeros((1, 1, 3))

    w, h = im.size

    # We need a float array between 0-1, rather than
    # a uint8 array between 0-255
    im = np.array(im).astype(np.float) / 255

    # With newer (1.0) versions of matplotlib, you can
    # use the "zorder" kwarg to make the image overlay
    # the plot, rather than hide behind it... (e.g. zorder=10)
    fig.figimage(im, (0.95 * fig.bbox.xmax) - w, (0.96 * fig.bbox.ymax) - h)

    # Annotate config file and creation date
    plt.figtext(0.950, 0.030, tc.time,
                ha="right", va="bottom", color="gray", size=8)
    plt.figtext(0.950, 0.04, tc.config_file,
                ha="right", va="bottom", color="gray", size=8)

    year = datetime.datetime.now().year
    text = "$\copyright$ {} Department of Energy".format(year)
    plt.figtext(0.750, 0.03, text,
                ha="left", va="bottom", color="gray", size=8)

    # Finish
    # --------------------------------------------------------------#
    save_file = getattr(tc, 'save_file', None)
    if save_file:
        if type(save_file) != 'str':
            # Then it's probably a bool from 'yes' or 'true' in the config
            save_file = tc.config_file.split('.')[0] + '.png'
            if save_file == 'config.png':
                save_file = 'temp.png'
        Notice.ok("Saving file " + save_file + "...", hold=True)
        plt.savefig(save_file, dpi=fig.dpi)
        Notice.ok("Done")
        Notice.warning("The displayed image is not identical to the saved one")
        plt.show()
    else:
        plt.show()
