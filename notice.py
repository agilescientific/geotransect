#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Defines various data containers for plotting a transect.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""


class Notice(object):
    """
    Helper class to make printout more readable.
    """
    styles = {'HEADER': '\033[95m',
              'INFO': '\033[94m',      # blue
              'OK': '\033[92m',        # green
              'WARNING': '\033[93m',   # red
              'FAIL': '\033[91m',
              'BOLD': '\033[1m'
              }
    ENDC = '\033[0m'

    def __init__(self, string, style):
        print self.styles[style.upper()] + string + self.ENDC

    @classmethod
    def title(cls):
        logo = """
                  _                                 _
                 | |                               | |
  __ _  ___  ___ | |_ _ __ __ _ _ __  ___  ___  ___| |
 / _` |/ _ \/ _ \| __| '__/ _` | '_ \/ __|/ _ \/ __| __|
| (_| |  __/ (_) | |_| | | (_| | | | \__ \  __/ (__| |
 \__, |\___|\___/ \__|_|  \__,_|_| |_|___/\___|\___|\__|
  __/ |
 |___/
 """
        return cls(logo, 'WARNING')

    @classmethod
    def warning(cls, string):
        return cls(string, 'WARNING')

    @classmethod
    def header(cls, string):
        return cls('\n'+string+'\n', 'HEADER')

    @classmethod
    def hr_header(cls, string):
        hr = "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
        return cls(hr+string.upper(), 'HEADER')

    @classmethod
    def info(cls, string):
        return cls('\n'+string, 'INFO')

    @classmethod
    def ok(cls, string):
        return cls(string, 'OK')
