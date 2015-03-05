#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Defines various data containers for plotting a transect.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import sys
import os

import numpy as np
from scipy import fft
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import obspy

filename = sys.argv[1]
segyfile = os.path.basename(filename)

# np.percentile and np.clip for best imshow display

# Read all traces
section = obspy.read(filename)

r_elevs = []
s_elevs = []
esp = []     # energy source point number
ens = []     # ensemble number

for t in section.traces:
    nsamples = t.stats.segy.trace_header.number_of_samples_in_this_trace
    dt = t.stats.segy.trace_header.sample_interval_in_ms_for_this_trace
    if dt > 100:
        dt /= 1000.
    r_elevs.append(t.stats.segy.trace_header.datum_elevation_at_receiver_group)
    s_elevs.append(t.stats.segy.trace_header.receiver_group_elevation)
    esp.append(t.stats.segy.trace_header.energy_source_point_number)
    ens.append(t.stats.segy.trace_header.ensemble_number)

ntraces = len(section.traces)
tbase = np.arange(0, nsamples * dt, dt)
tstart = 0
tend = tbase[-1]
aspect = float(ntraces) / float(nsamples)
nf = 1.0

print 'ntraces', ntraces
print 'nsamples', nsamples
print 'dt', dt/nf

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


def add_subplot_axes(ax, rect, axisbg='w'):
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


# MAIN PLOT
h = (tend-tstart) / 250.0
w = ntraces / 250.0

fig = plt.figure(figsize=(10, 10), facecolor='w')

# Seismic data
ax = fig.add_axes([0.05, 0.05, 0.9, 0.95])
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
ax.set_title(segyfile)

# Colourbar
extreme = max(np.amax(data), abs(np.amin(data)))
colorbar_ax = add_subplot_axes(ax, [0.075, 0.075, 0.025, 0.15])
fig.colorbar(im, cax=colorbar_ax)
colorbar_ax.text(1.15, 1.1, '%3.0f' % -extreme,
                 transform=colorbar_ax.transAxes,
                 ha='left',
                 va='top')
colorbar_ax.text(1.15, -0.05, '%3.0f' % extreme,
                 transform=colorbar_ax.transAxes,
                 ha='left', fontsize=10)
colorbar_ax.set_axis_off()

# Power spectrum
S = abs(fft(data[:, 1]))
faxis = np.fft.fftfreq(len(data[:, 1]), d=(1/nf)*dt*1e-6)
spec_ax = add_subplot_axes(ax, [0.50, 0.075, 0.2, 0.15])
spec_ax.plot(faxis[:len(faxis)//4],
             np.log10(S[0:len(faxis)//4]),
             'b', lw=2)
spec_ax.set_xlabel('frequency [Hz]', fontsize=10)
spec_ax.set_xticklabels([0, 100, 200, 300], fontsize=10)
# spec_ax.set_xticklabels(spec_ax.get_xticks(), fontsize=10)
spec_ax.set_yticklabels(spec_ax.get_yticks(), fontsize=10)
spec_ax.set_yticks([])
spec_ax.set_yticklabels([])
spec_ax.text(.95, .9, 'Power spectrum',
             horizontalalignment='right',
             transform=spec_ax.transAxes, fontsize=10
             )
spec_ax.grid('on')

# Histogram
hist_ax = add_subplot_axes(ax, [0.75, 0.075, 0.2, 0.15])
hist_line = hist_ax.hist(np.ravel(data),
                         bins=int(100.0 / (clip_val/largest)))
hist_ax.set_xlim(-clip_val, clip_val)
# hist_ax.set_xticklabels([])
hist_ax.set_yticks([])
hist_ax.set_xticklabels([])
hist_ax.set_ylim(hist_ax.get_ylim()[0], hist_ax.get_ylim()[1]),
hist_ax.set_yticks([])
hist_ax.text(.95, .9, 'Histogram',
             horizontalalignment='right',
             transform=hist_ax.transAxes, fontsize=10
             )

plt.show()
