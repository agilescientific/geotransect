#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Defines various data containers for plotting a transect.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import os
import re

import numpy as np


def all_files(directory, match=None):
    """
    Find files whose names match some regex.
    """
    for path, dirs, files in os.walk(directory):
        for f in files:
            if match:
                if not re.search(match, f, flags=re.IGNORECASE):
                    continue
            yield os.path.join(path, f)


def rolling_median(a, window):
    """
    Apply a moving mdeian filter.
    """
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    rolled = np.lib.stride_tricks.as_strided(a,
                                             shape=shape,
                                             strides=strides)
    rolled = np.median(rolled, -1)
    rolled = np.pad(rolled, window / 2, mode='edge')
    return rolled


def despike(curve, curve_sm, max_clip):
    """
    Remove spikes from a curve.
    """
    spikes = np.where(curve - curve_sm > max_clip)[0]
    spukes = np.where(curve_sm - curve > max_clip)[0]
    out = np.copy(curve)
    out[spikes] = curve_sm[spikes] + max_clip  # Clip at the max allowed diff
    out[spukes] = curve_sm[spukes] - max_clip  # Clip at the min allowed diff
    return out
