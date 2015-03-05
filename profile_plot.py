#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Defines various data containers for plotting a transect.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import sys

import numpy as np
from scipy import fft
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

from obspy.segy.core import readSEGY

filename = sys.argv[1]

# np.percentile and np.clip for best imshow display

# Read all traces
section = readSEGY(filename, unpack_headers=True)

print(section)

r_elevs = []
s_elevs = []
esp = []  # energy source point number
ens = []  # ensemble number

for i, trace in enumerate(section.traces):
    nsamples = trace.header.number_of_samples_in_this_trace
    dt = trace.header.sample_interval_in_ms_for_this_trace
    r_elevs.append(trace.header.datum_elevation_at_receiver_group)
    s_elevs.append(trace.header.receiver_group_elevation)
    esp.append(trace.header.energy_source_point_number)
    ens.append(trace.header.ensemble_number)

ntraces = len(section.traces)
tbase = np.arange(0, nsamples * dt / 1000.0, dt)
tstart = 0
tend = np.amax(tbase)
aspect = float(ntraces) / float(nsamples)
nf = 1.0

print ('ntraces', ntraces)
print ('nsamples', nsamples)
print ('dt', dt/nf)

data = np.zeros((nsamples, ntraces))
for i, trace in enumerate(section.traces):
    data[:, i] = trace.data

line_extents = {'first_trace': 1,
                'last_trace': ntraces,
                'start_time': tstart,
                'end_time': tend
                }

clip_val = np.percentile(data, 99.0)

print "clip_val", clip_val
print "max_val", np.amax(data)
print "min_val", np.amin(data)
print "tstart", tstart
print "tend", tend
largest = max(np.amax(data), abs(np.amin(data)))

# MAIN PLOT
h = (tend-tstart) / 250.0
w = ntraces / 250.0


def figsize(scale):
    fig_width_pt = 469.755                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = (fig_width, fig_height)
    return fig_size

sc = 2
page_dims = figsize(scale=sc)

# Allocate areas of the plot in square inches.
nrows = int(np.floor(page_dims[0]))   # number of complete inches in x
ncols = int(np.floor(page_dims[-1]))  # number of complete inches in y
print nrows, ncols


def make_grid_layout(nrows, ncols, data, s_elevs, clip_val, line_extents):
    gs = gridspec.GridSpec(ncols, nrows)
    fig = plt.figure(figsize=page_dims,
                     facecolor='w',
                     edgecolor='k',
                     dpi=None,
                     frameon=True)

    x = np.arange(0, ncols, 1)
    y = np.arange(0, nrows, 1)

    index = [(a, b) for a in x for b in y]
    # Make Background Grid
    for idx, g in zip(index, gs):
        ax = fig.add_subplot(g)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.patch.set_facecolor('grey')
        ax.patch.set_alpha(0.25)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.text(0.5, 0.5,
                str(idx[0])+', '+str(idx[1]),
                ha='center',
                va='center',
                size=sc*5,
                alpha=.5)

    # Make Overlay Grid
    # 'name', (start from top, span_y, start_from_left, end)
    my_axes = {'main': (2, 7, 2, -4),
               'left': (2, 7, 0, 2),
               'right': (3, 8, -4, 0),
               'top line plot': (0, 2, 2, 9),
               'base map': (0, 3, 9, 13),
               'logo': (0, 2, 0, 2),
               'bottom line plot': (-1, 0, 2, -4),
               'lower left info': (-1, 0, 0, 2),
               'lower right info': (7, 8, 9, 13),
               # 'inset a':(2, 7, 4, 5),
               # 'scale':(7, 9, 7, 8),
               # 'inset elevation':(1, 2, 2, 10)
               }

    for key, val in my_axes.iteritems():
        my_ax = fig.add_subplot(gs[val[0]:val[1], val[2]:val[3]])
        my_ax.set_xticks([])
        my_ax.set_yticks([])
        my_ax.patch.set_facecolor('red')
        my_ax.patch.set_alpha(0.2)
        my_ax.text(0.5, 0.5,
                   key,
                   ha='center',
                   va='center',
                   size=sc*8,
                   alpha=1.0)

    # Plot top line profile
    axtop = fig.add_subplot(gs[1:3, 2:9])
    axtop.plot(np.arange(len(s_elevs)), s_elevs, 'k')
    axtop.patch.set_alpha(0.1)
    axtop.set_ylabel('Elevation [m]')
    axtop.set_xticks([])
    axtop.set_xlim([0, len(s_elevs)])
    axtop.set_ylim([-np.amax(s_elevs), 3*np.amax(s_elevs)])
    # axtop.set_tight_layout()

    # Plot main
    xsec_ax = fig.add_subplot(gs[2:7, 2:-4])
    xsec_ax.imshow(data, cmap=cm.gray, origin='upper',
                   vmin=-clip_val,
                   vmax=clip_val,
                   # extent =  (line_extents['first_trace'],
                   #           line_extents['last_trace'],
                   #           line_extents['end_time'],
                   #           line_extents['start_time']),
                   # aspect = aspect*0.75
                   )

    # fig.subplots_adjust(wspace=0.02, hspace=0.02)
    # fig.tight_layout()
    fig.show()

    return fig, gs

fig, gs = make_grid_layout(nrows, ncols, data, s_elevs, clip_val, line_extents)

fig = plt.figure(figsize=(1.5*w, h), facecolor='w')
ax = fig.add_axes([0.1, 0.1, 0.8, 0.85])
im = ax.imshow(data, cmap=cm.gray, origin='upper',
               vmin=-clip_val,
               vmax=clip_val,
               extent=(line_extents['first_trace'],
                       line_extents['last_trace'],
                       line_extents['end_time'],
                       line_extents['start_time']),
               aspect = aspect*0.5
               )
ax.set_ylabel('Two-way time [ms]')
ax.set_xlabel('Trace no.')
ax.grid()
# ax.set_title(segyfile)

# TOP LINE PLOT
axtop = fig.add_axes([0.1, 0.8, 0.8, 0.025])
axtop.plot(np.arange(len(s_elevs)), s_elevs, 'k')
axtop.patch.set_alpha(0.1)
axtop.set_ylabel('Elevation [m]')
axtop.set_xticks([])
axtop.set_xlim([0, len(s_elevs)])
# axtop.set_tight_layout()

# put a colorbar
extreme = max(np.amax(data), abs(np.amin(data)))
colorbar_ax = fig.add_axes([0.8, 0.8, 0.01, 0.07])
fig.colorbar(im, cax=colorbar_ax)
colorbar_ax.text(0.5, -0.1, '%3.0f' % -extreme,
                 transform=colorbar_ax.transAxes,
                 ha='center',
                 verticalalignment='top')
colorbar_ax.text(0.5, 1.1, '%3.0f' % extreme,
                 transform=colorbar_ax.transAxes,
                 ha='center')
colorbar_ax.set_axis_off()

# Do a histogram
hist_ax = fig.add_axes([0.85, 0.35, 0.10, 0.10], axisbg='w')
hist_ax.patch.set_alpha(0.5)
# centering the plot to the same x-range as above the plot
hist_line = hist_ax.hist(np.ravel(data),
                         bins=int(100.0 / (clip_val/largest)),
                         alpha=0.5)
hist_ax.set_xlim(-clip_val, clip_val)
# hist_ax.set_xticklabels([])

hist_ax.set_xticklabels(hist_ax.get_xticks(), fontsize=6)
hist_ax.set_xlabel('amplitude', fontsize=6)
hist_ax.set_ylim(hist_ax.get_ylim()[0], hist_ax.get_ylim()[1]),
hist_ax.set_yticks([])
hist_ax.set_title('Histogram', fontsize=8)


# Do a powerspectrum
S = abs(fft(data[:, 1]))
faxis = np.fft.fftfreq(len(data[:, 1]), d=(1/nf)*dt*1e-6)
# Ss = np.zeros_like(S)
# for signal in np.rollaxis(data,0,1):
#    Ss += abs(fft(signal))
# Ss /= ntraces
spec_ax = fig.add_axes([0.85, 0.15, 0.1, 0.1], axisbg='w')
spec_ax.patch.set_alpha(0.5)
spec_ax.plot(faxis[:len(faxis)//4],
             np.log10(S[0:len(faxis)//4]),
             'b', lw=2,
             alpha=0.5)
spec_ax.set_xlabel('frequency [Hz]', fontsize=6)
spec_ax.set_xticklabels(spec_ax.get_xticks(), fontsize=6)
spec_ax.set_yticklabels(spec_ax.get_yticks(), fontsize=6)
spec_ax.set_ylabel('power [dB]', fontsize=8)
spec_ax.set_title('Power spectrum', fontsize=8)
spec_ax.grid('on')

fig.show()
