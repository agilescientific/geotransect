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
from collections import OrderedDict

import yaml

from containers import TransectContainer


def ordered_load(stream, Loader=yaml.Loader, object_pairs_hook=OrderedDict):
    """
    Trick to load YAML as OrderedDict instead of ordinary dict.

    From http://stackoverflow.com/questions/5121931/

    Copyright: https://github.com/coldfix
    License: CC-BY-SA
    """
    class OrderedLoader(Loader):
        pass

    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))

    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)

    return yaml.load(stream, OrderedLoader)


def complete_paths(dictionary, root):
    """
    Takes a dictionary whose values are relative file names and
    adds the root directory to them.

    Args:
        dictionary (dict): A dict of relative file names or directories.
        root (str): The absolute path to the relative paths in `dictionary`.
    """
    for k, v in dictionary.items():
        # This feels unpythonic, but it's the best I can come up with.
        # See discussion http://stackoverflow.com/questions/9168904/
        if isinstance(v, basestring):
            dictionary[k] = os.path.join(root, v)
        else:
            # v is iterable. If it isn't an error will be thrown.
            dictionary[k] = []
            for i in v:
                dictionary[k].append(os.path.join(root, i))
    return dictionary


def main(cfg):

    print "Initializing"
    params = cfg['params']
    root = cfg['params']['data_dir']
    data = complete_paths(cfg['data'], root)
    layers = complete_paths(cfg['map'], root)
    tc = TransectContainer(params, layers, data)

    print "Plotting"
    tc.plot()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--config",
                        metavar="file",
                        type=argparse.FileType('r'),
                        default="config.yaml",
                        help="The name of a YAML config file.",
                        nargs="?")
    c = parser.parse_args()

    print "config file", c.config

    with c.config as ymlfile:
        config = ordered_load(ymlfile, yaml.SafeLoader)

    main(config)
