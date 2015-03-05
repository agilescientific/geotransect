#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Defines various data containers for plotting a transect.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import os
import re


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
