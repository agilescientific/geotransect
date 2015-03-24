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


class Loader(yaml.Loader):
    """
    Overloads pyYAML Loader with ability to include files, and to load
    everything as an OrderedDict, so the order of params can matter, e.g.
    for specifying map layer data.

    From: http://stackoverflow.com/questions/528281
    and http://stackoverflow.com/questions/5121931

    Args:
        stream (YAML stream): A stream from yaml.load()
    """
    def __init__(self, stream):
        super(Loader, self).__init__(stream)
        self._root = os.path.split(stream.name)[0]
        self.add_constructor('!include', Loader.include)
        self.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                             Loader.mapping)

    def include(self, node):
        filename = os.path.join(self._root, self.construct_scalar(node))
        with open(filename, 'r') as f:
            return yaml.load(f, Loader)

    def mapping(self, node):
        self.flatten_mapping(node)
        return OrderedDict(self.construct_pairs(node))


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


def complete_vel_paths(dictionary, root):

    if dictionary["type"] != 'constant':
        dictionary["data"] = os.path.join(root, dictionary["data"])

    return dictionary


def complete_map_paths(dictionary, root):
    """
    The map dictionary is a bit more complicated.

    Args:
        dictionary (dict): A dict of relative file names or directories.
        root (str): The absolute path to the relative paths in `dictionary`.
    """
    for k, v in dictionary.items():
        dictionary[k]['file'] = os.path.join(root, v['file'])
    return dictionary


def complete_pot_paths(dictionary, root):
    """
    The potfield dictionary is more complicated still.

    Args:
        dictionary (dict): A dict of relative file names or directories.
        root (str): The absolute path to the relative paths in `dictionary`.
    """
    for k, v in dictionary.items():
        dictionary[k]['file'] = os.path.join(root, v['file'])
        if dictionary[k].get('colour'):
            fname = os.path.join(root, v['colour'])
            if os.path.isfile(fname):
                dictionary[k]['colour'] = fname
                dictionary[k]['colour_is_file'] = True
            else:
                dictionary[k]['colour_is_file'] = False
    return dictionary


def main(ymlfile):
    """
    Load a config file, instantiate a TransectContainer, and plot it.

    Note:
        Does not return; the plot is a side-effect.

    Args:
        ymlfile (file): A configuration file in YAML format.
    """

    with ymlfile as f:
        cfg = yaml.load(f, Loader)

    kwargs = {}
    cfg['params']['config_file'] = ymlfile.name
    kwargs['params'] = cfg['params']
    root = cfg['params']['data_dir']
    kwargs['data'] = complete_paths(cfg['data'], root)
    kwargs['velocity'] = complete_vel_paths(cfg['velocity'], root)
    kwargs['layers'] = complete_map_paths(cfg['map_layers'], root)
    kwargs['potfields'] = complete_pot_paths(cfg['potfields'], root)

    tc = TransectContainer(**kwargs)

    tc.plot()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config",
                        metavar="file",
                        type=argparse.FileType('r'),
                        default="config.yaml",
                        help="The name of a YAML config file.",
                        nargs="?")
    c = parser.parse_args()
    main(c.config)
