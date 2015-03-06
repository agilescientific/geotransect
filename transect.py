#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Master script to build transects.

Put configuration in config.yaml or pass the config file on
the command line.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""

import argparse
import os

import yaml

from containers import TransectContainer


def main(cfg):

    print "Initializing"

    params = cfg['params']
    data = cfg['data']

    root = cfg['params']['data_dir']
    for k, v in data.items():
        data[k] = os.path.join(root, v)

    tc = TransectContainer(params, data)

    # Make the plots
    print "Plotting"
    tc.plot()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--config",
                        type=argparse.FileType('r'),
                        help="The name of a YAML config file.")
    c = parser.parse_args(['--config', 'config.yaml'])

    with c.config as ymlfile:
        config = yaml.load(ymlfile)

    main(config)
