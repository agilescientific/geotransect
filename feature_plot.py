#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for plotting.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import csv
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms

import utils


def get_curve_params(abbrev, fname):
    """
    Builds and returns a dictionary of petrophysical parameters for
    plotting purposes.

    Args:
        abbrev (str): A curve mnemonic or other abbreviation.
        fname (str): The path to a file with the curve configuration.

    Returns:
        dict: A mapping of parameter:value for the curve in question.
    """
    params = {'acronym': abbrev}
    with open(fname, 'rU') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['acronymn'] == abbrev:
                params['track'] = int(row['track'])
                params['units'] = row['units']
                params['xleft'] = float(row['xleft'])
                params['xright'] = float(row['xright'])
                params['logarithmic'] = row['logarithmic']
                params['hexcolor'] = row['hexcolor']
                params['fill_left_cond'] = bool(row['fill_left_cond'])
                params['fill_left'] = row['fill_left']
                params['fill_right_cond'] = bool(row['fill_right_cond'])
                params['fill_right'] = row['fill_right']
                params['xticks'] = row['xticks'].split(',')

    return params


def get_tops(fname):
    """
    Takes a tops_dictionary for plotting in the logs tracks.

    Args:
        fname (str): The path to a file containing the tops.
    """
    tops = {}
    with open(fname) as f:
        for line in f.readlines():
            if not line.startswith('#'):
                temp = line.strip().split(',')
                tops[temp[0]] = float(temp[1])
    return tops


def plot_striplog(ax, striplog, width=1,
                  ladder=False, minthick=1,
                  alpha=0.75):

    pass


def plot_feature_well(tc, gs):
    """
    Plotting function for the feature well.

    Args:
        tc (TransectContainer): The container for the main plot.
        log (axis): A matplotlib axis.
        gs (GridSpec): A matplotlib gridspec.
    Note:
        Use with care! This is not a standalone plotting function. It was
        removed from main plotting method for readability and ease of use.
    """
    fname = tc.settings['curve_display']

    logs = tc.log.get(tc.feature_well)

    Z = logs.data['DEPT']
    curves = ['GR', 'DT',
              'DPHI_SAN',
              'NPHI_SAN',
              'DTS',
              'RT_HRLT',
              'RHOB',
              'DRHO']

    window = 51    # window length for smoothing must be an odd integer
    ntracks = 5
    lw = 1.0
    smooth = True
    has_striplog = True
    naxes = 0
    ncurv_per_track = np.zeros(ntracks)

    if has_striplog:
        ncurv_per_track[0] = 1

    for curve in curves:
        naxes += 1
        params = get_curve_params(curve, fname)
        ncurv_per_track[params['track']] += 1

    gs.update(wspace=0)

    axss = plt.subplot(gs[2:-1, -5])
    axs0 = [axss, axss.twiny()]
    axs1 = [plt.subplot(gs[2:-1, -4])]
    axs2 = [plt.subplot(gs[2:-1, -3])]
    axs3 = [plt.subplot(gs[2:-1, -2])]
    axs4 = [plt.subplot(gs[2:-1, -1])]

    axs = [axs0, axs1, axs2, axs3, axs4]

    #striplog = tc.striplog.get(tc.feature_well)
    #plot_striplog(axs0[0], striplog, width=5, alpha=0.75,
    #             ladder=True)

    axs0[0].set_ylim([Z[-1], 0])
    label_shift = np.zeros(len(axs))

    for curve in curves:
        values = logs.data[curve]
        params = get_curve_params(curve, fname)
        i = params['track']

        j = 0

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
            values = utils.rolling_median(values, window)

        if curve == 'GR':
            j = 1  # second axis in first track
            label_shift[i] = 1
            if params['fill_left_cond']:
                # do the fill for the lithology track
                axs[i][j].fill_betweenx(Z, params['xleft'], values,
                                        facecolor=params['fill_left'],
                                        alpha=1.0, zorder=11)

        if (curve == 'DPHI_SAN') and params['fill_left_cond']:
            # do the fill for the neutron porosity track
            nphi = utils.rolling_median(logs.data['NPHI_SAN'], window)
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

        if curve == 'DRHO':
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
        # curve labels
        # ------------------------------------------------- #

        trans = transforms.blended_transform_factory(axs[i][j].transData,
                                                     axs[i][j].transData)

        axs[i][j].text(xpos, -130 - 40*(label_shift[i]-1),
                       curve,
                       horizontalalignment='center',
                       verticalalignment='bottom',
                       fontsize=12, color=params['hexcolor'],
                       transform=trans)
        # curve units
        if label_shift[i] <= 1:
            axs[i][j].text(xpos, -110,
                           units,
                           horizontalalignment='center',
                           verticalalignment='top',
                           fontsize=12, color='k',
                           transform=trans)

        # ------------------------------------------------- #
        # scales and tickmarks
        # ------------------------------------------------- #

        axs[i][j].set_xscale(linOrlog)
        axs[i][j].set_ylim([Z[-1], 0])
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

    # ------------------------------------------------- #
    # End of curve loop
    # ------------------------------------------------- #

    # Add Depth label
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

    tops_fname = os.path.join(tc.data_dir,
                              tc.tops_file)
    tops = get_tops(tops_fname)
    topx = get_curve_params('DT', fname)
    topmidpt = np.amax((topx)['xright'])

    # plot tops
    for i in range(ntracks):

        for mkr, depth in tops.iteritems():

            # draw horizontal bars at the top position
            axs[i][-1].axhline(y=depth,
                               xmin=0.01, xmax=.99,
                               color='b', lw=2,
                               alpha=0.5,
                               zorder=100)

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
