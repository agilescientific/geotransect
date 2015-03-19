#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Crawl a tree for SEGY files, read their trace locations, and
make a shapefile containing the surface traces of the lines.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import os
import sys
import logging

import obspy
import fiona
from fiona import crs
from shapely.geometry import Point, LineString, mapping
import pyproj as pp

import utils

# Set up logging.
log = logging.getLogger('lithlog')
fstring = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

fh = logging.FileHandler('/tmp/lithology.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(logging.Formatter(fstring))
log.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
ch.setFormatter(logging.Formatter(fstring))
log.addHandler(ch)

utm_nad27 = pp.Proj("+init=EPSG:26720")
utm_nad83 = pp.Proj("+init=EPSG:26920")


class ShapeFileExists(Exception):

    pass

def sgy2shp(input_dir, output_dir, convert=False):
    """
    Extracts trace location from SEGY files and saves it in a
    shape file. A shape file is generated for each SEGY file.

    Returns nothing, side effect: writes the shape files.

    Args:
        input_dir (str): Directory containing SEGY files.
        output_dir (str): Directory to save shape files.
    """

    line_out_file = os.path.join(output_dir, "seismic_lines.shp")

    
    if os.path.exists(line_out_file):
        raise ShapeFileExists
    
    # Set up the shapefile schema.
    line_schema = {'geometry': 'LineString',
                   'properties': {'segyfile': 'str',
                                  'line': 'str'
                                  }
                   }

    with fiona.open(line_out_file, "w",
                    driver="ESRI Shapefile",
                    crs=crs.from_epsg(26920),
                    schema=line_schema) as line_out:

        for path in utils.walk(input_dir, "\\.se?gy$"):

            filebase = os.path.splitext(os.path.basename(path))[0]

            # Read in the headers.
            segy = obspy.read(path,
                              headonly=True,
                              unpack_trace_header=True)

            points = []

            point_out_file = os.path.join(output_dir, "." +
                                          filebase + '.shp')

            # Set up the shapefile schema.
            point_schema = {'geometry': 'Point',
                            'properties': {'line': 'str',
                                           'segyfile': 'str',
                                           'trace': 'int'
                                           }
                            }

            with fiona.open(point_out_file, "w",
                            driver="ESRI Shapefile",
                            crs=crs.from_epsg(26920),
                            schema=point_schema) as trace_out:

                for i, trace in enumerate(segy):

                    header = trace.stats.segy.trace_header
                    scalar = header.scalar_to_be_applied_to_all_coordinates
                    if scalar == -100:
                        gain = 0.01
                    elif scalar == -10:
                        gain = 0.1
                    else:
                        gain = 1.0

                    x = float(header.source_coordinate_x) * gain
                    y = float(header.source_coordinate_y) * gain

                    # Sanity check for geometry order of magnitude.
                    if x > 9e5 or y > 55e6:
                        if x > 9e6 or y > 55e7:
                            log.info('Found weird coords: dividing by 100')
                            x = x / 100.0
                            y = y / 100.0
                        else:
                            log.info('Found weird coords: dividing by 10')
                            x = x / 10.0
                            y = y / 10.0
                    else:
                        pass

                    if convert:
                        log.info("Converting from NAD27 to NAD83")
                        x, y = pp.transform(utm_nad27, utm_nad83, x, y)

                    p = Point(x, y)
                    points.append(p)
                    trace_out.write({'geometry': mapping(p),
                                     'properties': {'line': filebase,
                                                    'segyfile': path,
                                                    'trace': i}
                                     })

            #linestring = LineString(points)
            #line_out.write({'geometry': mapping(linestring),
            #                'properties': {'segyfile': path,
            #                               'line': filebase}
            #                })


def main():

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]

    sgy2shp(input_dir, output_dir)

if __name__ == "__main__":
    main()
