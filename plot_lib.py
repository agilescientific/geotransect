import numpy as np

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

from matplotlib import pyplot as plt

from matplotlib.colors import hsv_to_rgb, hex2color

import matplotlib.transforms as transforms
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle

from matplotlib import pyplot as plt, gridspec

from lithology.legends import legend

import csv
import StringIO
import re

from PIL import Image

def float_formatter(x):
    """
    two decimal places of precision for a floating pt number
    """
    x = float(x)
    return lambda x: "%.2f" % x

def get_tops(fname):
    """
    takes a tops_dictionary for plotting
    the all petrophysical tracks and wells
    """
    tops = {}
    with open(fname) as f:
        for line in f.readlines():
            if not line.startswith('#'):
                temp = line.strip().split(',')
                tops[temp[0]] = float(temp[1])
    return tops


def plot_striplog(ax, striplog, width=1,
                  ladder=False, summaries=False, minthick=1,
                  alpha=0.75):

    for t, b, l in zip(striplog['tops'], striplog['bases'],
                       striplog['liths']):

        origin = (0, t)
        colour = '#' + l['map colour'].strip()
        thick = b - t

        if ladder:
            w = legend[colour[1:]]['width']
        else:
            w = width

        rect = Rectangle(origin, w, thick, color=colour,
                         edgecolor='k',
                         linewidth=1.0, alpha=alpha)
        ax.add_patch(rect)

        if summaries:
            if thick >= minthick:
                ax.text(6, t + thick / 2, summarize(l),
                        verticalalignment='center')

        ax.set_ylim([striplog['bases'][-1], 0])


def get_curve_params(abbrev, fname):
    """
    returns a dictionary of petrophysical parameters for
    plotting purposes
    """
    params = {}
    with open(fname, 'rU') as csvfile:
        reader = csv.DictReader(csvfile) 
        for row in reader:
            if row['acronymn'] == abbrev:
                params['acronymn'] = row['acronymn']
                params['track'] = int(row['track'])
                params['units'] = row['units']
                params['xleft'] = float(row['xleft'])
                params['xright'] = float(row['xright'])
                params['logarithmic'] = row['logarithmic']
                params['hexcolor'] = (row['hexcolor'])
                params['fill_left_cond'] = bool(row['fill_left_cond'])
                params['fill_left'] = (row['fill_left'])
                params['fill_right_cond'] = bool(row['fill_right_cond'])
                params['fill_right'] = (row['fill_right'])
                params['xticks'] = row['xticks'].split(',')
    return params


def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    rolled = np.lib.stride_tricks.as_strided(a,
                                             shape=shape,
                                             strides=strides)
    rolled = np.median(rolled, -1)
    rolled = np.pad(rolled, window / 2, mode='edge')
    return rolled

def despike(curve, curve_sm, max_clip):
    spikes = np.where(curve - curve_sm > max_clip)[0]
    spukes = np.where(curve_sm - curve > max_clip)[0]
    out = np.copy(curve)
    out[spikes] = curve_sm[spikes] + max_clip  # Clip at the max allowed diff
    out[spukes] = curve_sm[spukes] - max_clip  # Clip at the min allowed diff
    return out




def uberPlot(transect, seismic_container, elevation_container,
             log_container, bedrock_container, striplog_container,
             em_container, extents):

    # -------------------------------------------------#
    # get the data from the containers
    # --------------------------------------------------#
    transectx, transecty = transect.coords.xy

    seismic_data = seismic_container.data  # list of 2D arrays
    seismic_x = seismic_container.coords  # list of 1D arrays Range along transect

    elevation_data = elevation_container.data  # 1D array
    elevation_x = elevation_container.coords  # Range along transect
    max_height = np.amax(elevation_container.elevation_profile)

    log_data = log_container.data  # List of LASReader objs
    log_x = log_container.coords  # List of locations for each log

    bedrock_data = bedrock_container.data  # List of dictionaries
    bedrock_x = bedrock_container.coords  # Range along transect

    striplog_data = striplog_container.data  # List of dicts
    striplog_x = striplog_container.coords  # Range along transect
    striplog = striplog_container.get("P-129")

    em_data = em_container.data
    em_x = em_container.coords

    # ------------------------------------------------------------#
    # Figure Layout
    # ------------------------------------------------------------#

    h = 15
    mw = 16  # width of main section (inches) must be div by 4
    fw = 5   # width of the feature plot (inches) must be div by 5

    gs = gridspec.GridSpec(h, mw + fw + 1)
    fig = plt.figure(facecolor='w', figsize=(mw + fw + 1, 15),
                     edgecolor='k',
                     dpi=None,
                     frameon=True)

    header = fig.add_subplot(gs[0:1, 0: mw])
    description = fig.add_subplot(gs[1:4, 0:2 * mw / 4])
    large_scale = fig.add_subplot(gs[1:4, 2 * mw / 4:3 * mw / 4])
    small_scale = fig.add_subplot(gs[1:4, 3 * mw / 4:4 * mw / 4])
    feat_header = fig.add_subplot(gs[0:1, -5:])
    ns_label = fig.add_subplot(gs[-1:, -5:])

    # featured = fig.add_subplot(gs[0:-1, 10:])

    elevation = fig.add_subplot(gs[5, :mw])
    xsection = fig.add_subplot(gs[6:h - 1, :mw])
    log = fig.add_subplot(gs[6: h - 1, :mw])
    em = fig.add_subplot(gs[h - 1:, :mw])

    # ------------------------------------------------------------ #
    # Adjust white space between subplots (maybe put in function?)
    # ------------------------------------------------------------ #
    left = 0.05     # left side of the subplots of the figure
    right = 0.95    # right side of the subplots of the figure
    bottom = 0.1   # bottom of the subplots of the figure
    top = 0.9      # top of the subplots of the figure
    wspace = 0.1   # amount of width reserved for blank space between subplots
    hspace = 0.05   # amount of height reserved for white space between subplots
    # adjust white space between subplots

    fig.subplots_adjust(left, bottom, right, top, wspace, hspace)

    # ------------------------------------------------------------#
    # Plot the seismic cross section
    # ------------------------------------------------------------#

    # Loop through each seismic line
    for coords, data in zip(seismic_x, seismic_data):

        z0 = 0
        depth = 2500
        xsection.imshow(data,
                        extent=[np.amin(coords) / 1000.0,
                                np.amax(coords) / 1000.0,
                                depth, z0],
                        aspect="auto", cmap="Greys")

        plot_axis = [extents[0] / 1000., extents[1] / 1000.,
                     extents[2], extents[3]]
        xsection.axis(plot_axis)

    xsection.set_ylabel("Depth [m]", fontsize=8)
    xsection.set_xlabel("Transect range [km]", fontsize=8)

    xsection.tick_params(axis='y', which='major', labelsize=8)
    xsection.tick_params(axis='y', which='minor', labelsize=8)

    xsection.tick_params(axis='x', which='major', labelsize=8)
    xsection.tick_params(axis='x', which='minor', labelsize=8)

    xsection.grid(True)

    # --------------------------------------------------------#
    # Plot the log overlays
    # --------------------------------------------------------#
    for las, pos in zip(log_data, log_x):

        data = np.nan_to_num(las.data["GR"])

        # normalize
        data /= np.amax(data)

        # scale position
        lgsc = 0.015  # hack to compress the log width
        data *= lgsc * (extents[1] - extents[0])
        data += pos

        log.plot(data, las.data['DEPT'],
                 'g', lw=0.5, alpha=0.5)

        log.set_xlim((extents[0], extents[1]))
        log.set_ylim((extents[2], extents[3]))

        log.axis("off")

    log.axis('off')

    # TODO WELL LOG LABELS

    # ----------------------------------------------------------#
    # Elevation and bedrock plot
    # -----------------------------------------------------------#
    # Could be done with LineCollection?
    for i, geo in enumerate(bedrock_data[:-1]):

        lim1 = bedrock_x[i]
        lim2 = bedrock_x[i + 1]

        idx = np.where(np.logical_and(elevation_x >= lim1,
                                      elevation_x <= lim2))[0]  # list

        if idx[-1] < elevation_x.size - 1:
            idx = np.append(idx, (idx[-1] + 1))

        hsv = np.array([[geo["AV_HUE"], geo["AV_SAT"],
                         geo["AV_VAL"]]]).reshape(1, 1, 3)

        color = hsv_to_rgb(hsv / 255.)
        elevation.plot(elevation_x[idx], elevation_data[idx],
                       linewidth=3, color=color.flatten())

    bbox = {'fc': 'w', 'pad': 0, 'ec': 'none', 'alpha': 0.5}
    props = {'ha': 'left', 'va': 'center', 'fontsize': 6, 'bbox': bbox}

    elevation.set_ylim((0, max_height))
    elevation.set_xlim(extents[:2])
    elevation.set_yticks([0, int(max_height),
                          int(np.amax(elevation_data))])
    elevation.set_xticks([])

    elevation.tick_params(axis='y', which='major', labelsize=8)

    elevation.patch.set_alpha(0.1)

    elevation.set_ylabel("Elevation [m]", fontsize=8)

    elevation.grid(True)
    elevation.text(0.0, .5 * (max_height), "Surface Geology",
                   props, rotation=0)
    elevation.set_frame_on(False)

    # ---------------------------------------------------------#
    # Header
    # ---------------------------------------------------------#

    header.axis("off")
    props["fontsize"] = 30
    dy = 0.2
    header.text(0.0, 0.5 + dy, "Transect Title",
                props,
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
                ("some meta data \n" +
                 "about this transect"),
                fontsize=14,
                horizontalalignment='right',
                verticalalignment='bottom')

    # -----------------------------------------------------#
    # Paragraph Description
    # -----------------------------------------------------#
    for i in [0.0, 0.5]:
        # put the paragraph in twice (2 columns)
        description.text(i, 1.0,
                         ("Lorem ipsum dolor sit amet, consectetur \n" +
                          "adipiscing elit, sed do eiusmod tempor \n" +
                          "incididunt ut labore et dolore magna aliqua. \n" +
                          "Ut enim ad minim veniam, quis nostrud \n" +
                          "exercitation ullamco laboris nisi ut aliquip \n" +
                          "ex ea commodo consequat. Duis aute irure \n"),
                         horizontalalignment='left',
                         verticalalignment='top',
                         fontsize=12
                         )

    description.axis('off')
    # description.set_xticks([])
    # description.set_yticks([])

    # TODO figure out how to get a description

    # --------------------------------------------------------#
    # Large scale map
    # ---------------------------------------------------------#
    large_scale.plot(transectx, transecty)
    large_scale.patch.set_facecolor("0.95")
    large_scale.patch.set_edgecolor("0.90")
    large_scale.set_xticks([])
    large_scale.set_yticks([])

    # large_scale.axis("off")

    # --------------------------------------------------------#
    # Small scale map
    # ---------------------------------------------------------#
    small_scale.plot(transectx, transecty)
    small_scale.patch.set_facecolor("0.95")
    small_scale.patch.set_edgecolor("0.90")
    small_scale.set_xticks([])
    small_scale.set_yticks([])

    # -----------------------------------------------------------#
    # Em data
    # -----------------------------------------------------------#

    em.plot(em_x, em_data)
    em.axis("off")

    # TODO

    # -----------------------------------------------------#
    # Feature plot header
    # -----------------------------------------------------#
    feat_header.text(0.0, 1.0,
                     ("Well " + "P-129"),
                     verticalalignment='top',
                     horizontalalignment='left',
                     fontsize=14,
                     fontweight='bold'
                     )
    feat_header.axis("off")
    # -----------------------------------------------------#
    #      Feature Plot - Oh Baby
    # -----------------------------------------------------#
    fname = 'templates/Petrophysics_display_template.csv'
    params = get_curve_params('RT_HRLT', fname)

    logs = log_container.get("P-129")

    Z = logs.data['DEPT']
    GR = logs.data['GR']
    DT = logs.data['DT']
    DPHISS = logs.data['DPHI_SAN']
    NPHISS = logs.data['NPHI_SAN']
    DTS = logs.data['DTS']
    RT_HRLT = logs.data['RT_HRLT']
    RHOB = logs.data['RHOB']
    DRHO = logs.data['DRHO']

    log_dict = {'GR': GR,
                'DT': DT,
                'DPHI_SAN': DPHISS,
                'NPHI_SAN': NPHISS,
                'DTS': DTS,
                'RT_HRLT': RT_HRLT,
                'RHOB': RHOB,
                'DRHO': DRHO
                }

    left = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.1   # the amount of width reserved for blank space between subplots
    hspace = 0.5   # the amount of height reserved for white space between subplots

    window = 51  # window length for smoothing must be an odd integer
    frac = 0.05
    ntracks = 5
    lw = 1.0
    smooth = True
    has_striplog = True
    height = 2.5 * ntracks  # in inches
    width = 1.5 * ntracks  # in inches
    fs = 12  # font size for curve labels
    naxes = 0
    ncurv_per_track = np.zeros(ntracks)

    if has_striplog:
        ncurv_per_track[0] = 1

    for curve, values in log_dict.iteritems():
        naxes += 1
        params = get_curve_params(curve, fname)
        ncurv_per_track[params['track']] += 1

    axss = plt.subplot(gs[2:-1, -5])
    axs0 = [axss, axss.twiny()]
    axs1 = [plt.subplot(gs[2:-1, -4])]
    axs2 = [plt.subplot(gs[2:-1, -3])]
    axs3 = [plt.subplot(gs[2:-1, -2])]
    axs4 = [plt.subplot(gs[2:-1, -1])]

    axs = [axs0, axs1, axs2, axs3, axs4]

    plot_striplog(axs0[0], striplog, width=5, alpha=0.75,
                  ladder=True)

    axs0[0].set_ylim([striplog['bases'][-1], 0])

    # Plot each curve with a white fill to fake the curve fill.

    label_shift = np.zeros(len(axs))

    for curve, values in log_dict.iteritems():

        params = get_curve_params(curve, fname)

        i = params['track']

        if i == 0:
            j = 1

        j = 0  # default number of tracks to index into

        ncurves = ncurv_per_track[i]

        label_shift[i] += 1

        units = '$%s$' % params['units']

        linOrlog = params['logarithmic']

        sxticks = np.array(params['xticks'])
        xticks = np.array(sxticks, dtype=float)
        whichticks = 'major'

        if linOrlog == 'log':
            midline = np.log(np.mean(xticks))
            xpos = midline
            whichticks = 'minor'
        else:
            midline = np.mean(xticks)
            xpos = midline

        if smooth:
            values = rolling_window(values, window)

        if params['acronymn'] == 'GR':
            j = 1  # second axis in first track
            label_shift[i] = 1

        # fill_left
        if params['fill_left_cond']:
            print params['fill_left_cond']
            if params['acronymn'] == 'GR':
                # do the fill for the lithology track
                axs[i][j].fill_betweenx(Z, params['xleft'], values,
                                        facecolor=params['fill_left'],
                                        alpha=1.0, zorder=11)

            if params['acronymn'] == 'DPHI_SAN':
                # do the fill for the neutron porosity track
                nphi = rolling_window(log_dict['NPHI_SAN'], window)
                axs[i][j].fill_betweenx(Z,
                                        nphi,
                                        values,
                                        where=nphi >= values,
                                        facecolor=params['fill_left'],
                                        alpha=1.0,
                                        zorder=11)

                axs[i][j].fill_betweenx(Z,
                                        nphi,
                                        values,
                                        where=nphi <= values,
                                        facecolor='#8C1717',
                                        alpha=0.5,
                                        zorder=12)

        if params['acronymn'] == 'DRHO':
            blk_drho = 3.2
            values += blk_drho   # this is a hack to get DRHO on RHOB scale
            axs[i][j].fill_betweenx(Z,
                                    blk_drho,
                                    values,
                                    where=nphi <= values,
                                    facecolor='#CCCCCC',
                                    alpha=0.5,
                                    zorder=12)

        # fill right
        if params['fill_right_cond']:

            axs[i][j].fill_betweenx(Z, values, params['xright'],
                                    facecolor=params['fill_right'],
                                    alpha=1.0, zorder=12)

        # plot curve
        axs[i][j].plot(values, Z, color=params['hexcolor'],
                       lw=lw, zorder=13)

        # set scale of curve
        axs[i][j].set_xlim([params['xleft'], params['xright']])

        # ------------------------------------------------- #
        # curve label
        # ------------------------------------------------- #

        trans = transforms.blended_transform_factory(axs[i][j].transData,
                                                     axs[i][j].transData)

        axs[i][j].text(xpos, -130 - 40 * (label_shift[i] - 1),
                       params['acronymn'],
                       horizontalalignment='center',
                       verticalalignment='bottom',
                       fontsize=12, color=params['hexcolor'],
                       transform=trans)
        # curve units
        if label_shift[i] <= 1:
            axs[i][j].text(xpos, -110, units,
                           horizontalalignment='center',
                           verticalalignment='top',
                           fontsize=12, color='k',
                           transform=trans)

        axs[i][j].set_xscale(linOrlog)

        axs[i][j].set_ylim([striplog['bases'][-1], 0])

        # set the x-tick positions and tick labels
        axs[i][j].axes.xaxis.set_ticks(xticks)
        axs[i][j].axes.xaxis.set_ticklabels(sxticks, fontsize=8)

        for label in axs[i][j].axes.xaxis.get_ticklabels():
            label.set_rotation(90)

        axs[i][j].tick_params(axis='x', direction='out')

        axs[i][j].xaxis.tick_top()
        axs[i][j].xaxis.set_label_position('top')

        axs[i][j].xaxis.grid(True, which=whichticks,
                             linewidth=0.25, linestyle='-',
                             color='0.75', zorder=100)

        axs[i][j].yaxis.grid(True, which=whichticks,
                             linewidth=0.25, linestyle='-',
                             color='0.75', zorder=100)

        axs[i][j].yaxis.set_ticks(np.arange(0, max(Z), 100))

        if i != 0:
                axs[i][j].set_yticklabels("")

    # add Depth label

    axs[0][0].text(-0.25, 50, 'MD \n $m$ ', fontsize='10',
                   horizontalalignment='left',
                   verticalalignment='center')

    axs[0][0].axes.yaxis.get_ticklabels()
    axs[0][0].axes.xaxis.set_ticklabels('')

    for label in axs[0][0].axes.yaxis.get_ticklabels():
            label.set_rotation(90)
            label.set_fontsize(10)

    for label in axs[1][0].axes.xaxis.get_ticklabels():
        label.set_rotation(90)
        label.set_fontsize(10)

    tops_fname = 'data/wells/P-129/tops/P-129_tops.txt'
    tops = get_tops(tops_fname)

    topx = get_curve_params('DT', fname)

    topmidpt = np.amax((topx)['xright'])

    # plot tops
    for i in range(ntracks):

        for mkr, depth in tops.iteritems():

            # draw horizontal bars at the top position
            axs[i][-1].axhline(y=depth, color='b', lw=2,
                               alpha=0.5, xmin=0.01, xmax=.99,
                               zorder=10000)

            # draw text box at the right edge of the last track
            axs[-1][-1].text(x=topmidpt, y=depth, s=mkr,
                             alpha=0.5, color='k',
                             fontsize='8',
                             horizontalalignment='center',
                             verticalalignment='center',
                             zorder=10000,
                             bbox=dict(facecolor='white', edgecolor='k',
                                       alpha=0.25, lw=0.25),
                             weight='light')

    # --------------------------------------------------------------#
    # LOGO
    # --------------------------------------------------------------#

    img = Image.open("logo.png")
    arr = np.array(img)

    # Broken
    # ns_label.imshow(arr, aspect="auto")
    # ns_label.axis('off')

    ns_label.text(0.9, 0.7, "L O G O",
                  verticalalignment='top',
                  horizontalalignment='right')

    ns_label.text(0.1, 0.7,
                  ("Department of Energy \n" +
                   "Nova Scotia, Canada"),
                  verticalalignment='top',
                  horizontalalignment='left')

    # horizontal line
    ns_label.axhline(y=0.7,
                     xmin=0.1,
                     xmax=0.9,
                     linewidth=1,
                     color='k')
    ns_label.axis("off")