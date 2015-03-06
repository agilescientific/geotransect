#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Defines various data containers for plotting a transect.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import os
import fnmatch
from functools import partial

import matplotlib.pyplot as plt
import numpy as np

import pyproj as pp
import rasterio
import fiona
from shapely.geometry import shape, Point, LineString
from shapely.ops import transform
from shapely.prepared import prep
from obspy.segy.core import readSEGY

from las import LASReader
from uberplot import uber_plot
from lithology.lithology import intervals_from_las3_string
from sgy2shp import sgy2shp
from utils import all_files


class BaseContainer(object):
    """
    Holds some basic information that we want in every object. Does not
    contain any data or pointers, only the transect parameters. Maybe
    eventually we can abstract some of the methods here too.
    """
    def __init__(self, params):

        for k, v in params.items():
            setattr(self, k, v)

        # The x extent will be updated at plot time.
        # TODO Make this less cryptic, or use two objects.
        self.extents = [0, 0, self.depth[1], self.depth[0]]


class TransectContainer(BaseContainer):
    """
    Main driver class for generating transect plots. Builds and
    accesses all other plot containers

    usage:
    tc = TransectContainer(params, data)

    @param params: Directory containing shape files of the transects.
    @param data: Contains the shape files of the SEGY trace headers.

    @returns a transectContainer object.
    """

    def __init__(self, params, data):

        super(TransectContainer, self).__init__(params)

        self.tops_file = data['tops_file']
        self.seismic = SeismicContainer(data['seismic_dir'], params)
        self.elevation = ElevationContainer(data['elevation_file'], params)
        self.log = LogContainer(data['well_dir'], params)
        self.bedrock = BedrockContainer(data['bedrock_dir'], params)
        self.striplog = StriplogContainer(data['striplog_dir'], params)
        self.potfield = PotfieldContainer(data['potfield_dir'], params)
        self.locmap = LocmapContainer(data['locmap_contents'], params)

        # Set up 'data' — the transect line — from shapefile.
        self.data = None
        with fiona.open(data['transect_file']) as c:
            for line in c:
                if line['properties']['id'] == self.id:
                    self.data = shape(line['geometry'])
        if not self.data:
            print "No transect with ID", self.id

    def plot(self):
        """
        Generates plot for the transect.
        """

        # Set the extent to the length? Or keep them all the same?
        self.extents[0] = 0
        self.extents[1] = self.data.length

        # update the containers
        self.seismic.update(self.data)
        self.log.update(self.data)
        self.elevation.update(self.data)
        self.bedrock.update(self.data)
        self.striplog.update(self.data)
        self.potfield.update(self.data)
        self.locmap.update(self.data)

        uber_plot(self)

        plt.show()


class LocmapContainer(BaseContainer):
    """
    Class for building and plotting a location map.

    USAGE:
    seis = seismicContainer(seis_dir)

    @param seis_dir: Input directory containing seismic shapefiles.
                     The shapefiles contain points corresponding to
                     UTM trace locations with fields segyfile,
                     trace.

    @returns a seismicContainer object
    """
    def __init__(self, contents, params):

        # First generate the parent object.
        super(LocmapContainer, self).__init__(params)

        # contents is a list of vector or raster layers to plot.
        # The list is in order, the first item will be the top layer.

        self.contents = contents
        self.mid = None

    def update(self, transect):

        pad = 5000  # m ... interior padding for map box

        ll_nad83 = pp.Proj("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs")
        utm_nad83 = pp.Proj("+init=EPSG:26920")
        project = partial(pp.transform, utm_nad83, ll_nad83)

        bounds = transect.bounds
        llx, lly, urx, ury = bounds
        w, h = urx - llx, ury - lly

        # Guarantees the map with have 8:3 aspect
        x_adj, y_adj = 0, 0
        if h > (w*3/8):
            print "adjusting x"
            x_adj = ((8*h/3.) - w) / 2.  # Aspect is hard-coded in uberplot
        else:                            # TODO Fix that!
            print "adjusting y"
            y_adj = ((3*w/8.) - h) / 2.

        ll = transform(project, Point(llx-pad-x_adj, lly-pad-y_adj))
        ur = transform(project, Point(urx+pad+x_adj, ury+pad+y_adj))
        self.ll, self.ur = ll, ur
        self.mid = Point(ll.x + 0.5*(ur.x-ll.x), ll.y + 0.5*(ur.y - ll.y))

        # the goal is to send data to uberplot to plot
        # we just want to find / compute data here
        # we need to distinguish between shp and raster

        for layer in self.contents:
            pass


class SeismicContainer(BaseContainer):
    """
    Class for reading and plotting seismic data.

    USAGE:
    seis = seismicContainer(seis_dir)

    @param seis_dir: Input directory containing seismic shapefiles.
                     The shapefiles contain points corresponding to
                     UTM trace locations with fields segyfile,
                     trace.

    @returns a seismicContainer object
    """
    def __init__(self, seis_dir, params):

        # First generate the parent object.
        super(SeismicContainer, self).__init__(params)

        self.lookup = {}     # Look up table for points: segy/trace
        self.data = []       # plotting data
        self.coords = []     # transect coords of plot data

        # This creates a (hidden) shapefile for each seismic line,
        # then steps over them to read their positions and meta.
        # TODO Simplify this... maybe don't even write the files.
        sgy2shp(seis_dir, seis_dir)

        for f in fnmatch.filter(os.listdir(seis_dir), '.*.shp'):
            shapefile = os.path.join(seis_dir, f)
            with fiona.open(shapefile, "r") as traces:
                for trace in traces:
                    self.lookup[shape(trace["geometry"])] = trace["properties"]

    def update(self, transect):
        """
        Updates the container data to traces that intersect
        the transect line.

        @param transect: A transect line as a shapely LineString
                         object.
        """
        # Preprocessing
        prepared = prep(transect.buffer(self.buffer))

        # Get intersecting points
        points = filter(prepared.contains, self.lookup.keys())

        # Reset data
        self.data = []
        self.coords = []

        # Lookup for grouping traces into segy files
        file_lookup = {}
        for point in points:
            meta = self.lookup[point]
            if(meta["segyfile"] in file_lookup):
                # project onto the transect line
                proj_d = transect.project(point)
                if proj_d:
                    file_lookup[meta["segyfile"]]["pos"].append(proj_d)
                    m = meta["trace"]
                    file_lookup[meta["segyfile"]]["trace"].append(m)
            else:  # Make a new dict entry
                file_lookup[meta["segyfile"]] = {}
                file_lookup[meta["segyfile"]]["trace"] = [meta["trace"]]
                proj_d = transect.project(point)
                file_lookup[meta["segyfile"]]["pos"] = [proj_d]

        # Read in the chunks from the segy file
        for segyfile in file_lookup.keys():
            segy = readSEGY(segyfile)
            traces = file_lookup[segyfile]["trace"]
            coords = file_lookup[segyfile]["pos"]

            # sort the traces to be in order
            idx = sorted(range(len(traces)), key=lambda k: traces[k])

            self.data.append(np.transpose(np.array(
                [segy.traces[traces[i]].data for i in idx])))

            self.coords.append(np.array([coords[i] for i in idx]))


class LogContainer(BaseContainer):
    """
    Container for managing and plotting LAS data.

    usage:

    lc = LogContainer(las_dir)

    @param las_dir: Directory shape files of LAS headers.

    @returns an LogContainer object.
    """
    def __init__(self, las_dir, params):

        # First generate the parent object.
        super(LogContainer, self).__init__(params)

        self.lookup = {}      # maps points to LAS filenames
        self.data = []        # plot data
        self.coords = []      # transect coords of plot data
        self.log_lookup = {}  # look up for log ids

        for shp in all_files(las_dir, '\\.shp$'):
            with fiona.open(shp, "r") as logs:
                for log in logs:
                    self.lookup[shape(log["geometry"])] = log["properties"]

    def update(self, transect):
        """
        Updates the container data to wells that intersect the
        transect line.

        @param transect: A transect line as a shapely LineString
                         object.
        """

        print "Updating log container",

        # Preprocess
        prepared = prep(transect.buffer(self.buffer))

        # Get the intersecting points
        points = filter(prepared.contains, self.lookup.keys())

        # reset the data
        self.data = []
        self.coords = []
        self.log_lookup = {}

        for point in points:
            meta = self.lookup[point]
            name = meta["name"]

            print name,

            # TODO: This has to be more dynamic
            filename = os.path.join('data', 'wells', name,
                                    'wireline_log', name +
                                    '_out.LAS')

            # Write out an error log?
            if not os.path.exists(filename):
                continue

            # store the transecting data
            self.data.append(LASReader(filename, null_subs=np.nan))
            self.coords.append(transect.project(point))
            self.log_lookup[name] = self.data[-1]

    def get(self, log_id):
        """
        Returns data corresponding to log_id
        """

        return self.log_lookup.get(log_id)


class ElevationContainer(BaseContainer):
    """
    Container for managing and plotting elevation data.

    usage:
    ec = ElevationContainer(elevation_raster)

    @param elevation_raster: Raster file containing the elevation
                             profile of the area of interest.

    @returns an ElevationContainer object.
    """
    def __init__(self, elevation_file, params):

        super(ElevationContainer, self).__init__(params)

        self.elevation_profile = []
        self.elevation_grid = []

        self.data = []
        self.coords = []

        decimate = 1

        with rasterio.drivers(CPL_DEBUG=True):
            with rasterio.open(elevation_file) as src:
                self.elevation_profile = src.read()[0, 0:-1:decimate,
                                                    0:-1:decimate]

                # Get as lat/lon using the affine coordinate transform
                # TODO - do this with proj? Why go via lo/la?
                la = np.arange(self.elevation_profile.shape[1])
                lo = np.arange(self.elevation_profile.shape[0])
                lat = la * src.affine[0]*decimate + src.affine[2]
                lon = lo * src.affine[4]*decimate + src.affine[5]
                wgs_grid = np.meshgrid(lat, lon)

                # TODO move to config
                ll_wgs84 = pp.Proj("+init=EPSG:4269")
                utm_nad83 = pp.Proj("+init=EPSG:26920")

                self.elevation_grid = pp.transform(ll_wgs84,
                                                   utm_nad83,
                                                   wgs_grid[0],
                                                   wgs_grid[1])

    def update(self, transect):
        """
        Updates the container data to a profile that intersect the
        transect line.

        @param transect: A transect line as a shapely LineString
                         object.
        """
        # transect coords need to be upsampled
        nsamples = 100

        self.coords = np.zeros(nsamples)
        self.data = np.zeros(nsamples)

        space = np.linspace(0, transect.length, nsamples)
        for i, n in enumerate(space):

            # interpolate along the transect
            x, y = transect.interpolate(n).xy

            # Get the closest elevation points
            xi = np.abs(self.elevation_grid[0][0, :] - x).argmin()
            yi = np.abs(self.elevation_grid[1][:, 0] - y).argmin()

            self.data[i] = self.elevation_profile[yi, xi]

            # add the distance to the coordinates
            self.coords[i] = n


class PotfieldContainer(BaseContainer):
    """
    Contains random data for placeholders.
    """
    def __init__(self, potfield_dir, params):
        super(PotfieldContainer, self).__init__(params)

        self.data = []
        self.coords = []

    def update(self, transect):
        self.coords = np.linspace(0, transect.length, 1000)
        self.data = np.random.randn(self.coords.size)

    def plot(self, extents, xticks):
        plt.plot(self.coords, self.data)
        plt.yticks([-4, 0, 4])
        plt.xticks(xticks, [])

        plt.ylabel("Anomaly [mGal]", fontsize=8)
        plt.tick_params(axis='y', which="major", labelsize=8)

        plt.grid(True)
        plt.xlim((extents[0], extents[1]))
        plt.ylim((-4, 4))


class BedrockContainer(BaseContainer):
    """
    Contains a geological map or similar basic geology shapes.
    """
    def __init__(self, bedrock_dir, params):

        super(BedrockContainer, self).__init__(params)

        self.lookup = {}
        self.data = []
        self.coords = []
        self.buffer = 1000  # [m]

        # Read in all shape files
        for f in os.listdir(bedrock_dir):
            if not f.endswith(".shp"):
                continue
            for line in fiona.open(os.path.join(bedrock_dir, f)):
                self.lookup[(shape(line['geometry']))] = line['properties']

    def update(self, transect):
        points = [Point(xy) for xy in transect.coords]
        self.data = []
        self.coords = []

        for polygon, properties in self.lookup.items():

            for p in filter(polygon.contains, points):

                self.data.append(properties)
                self.coords.append(transect.project(p))

        # Sort the data to be in order
        idx = sorted(range(len(self.coords)),
                     key=lambda k: self.coords[k])

        self.data = [self.data[i] for i in idx]
        self.coords = np.array([self.coords[i] for i in idx])


class StriplogContainer(LogContainer):

    def update(self, transect):
        # Preprocess
        prepared = prep(transect.buffer(self.buffer))

        # Get the intersecting points
        points = filter(prepared.contains, self.lookup.keys())

        # reset the data
        self.data = []
        self.coords = []
        self.log_lookup = {}

        for point in points:

            meta = self.lookup[point]

            name = meta["name"]
            filename = os.path.join('data', 'wells', name,
                                    'lithology_log', name +
                                    '_striplog.las')

            # Write out an error log?
            if not os.path.exists(filename):
                continue

            # store the transecting data
            data = LASReader(filename, null_subs=np.nan,
                             unknown_as_other=True)

            self.data.append(intervals_from_las3_string(data.other))
            self.coords.append(transect.project(point))
            self.log_lookup[name] = self.data[-1]

    def plot(self):

        self.data.plot()
