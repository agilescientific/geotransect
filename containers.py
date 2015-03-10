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

# Import required to avoid bug in Basemap
from mpl_toolkits.basemap import Basemap

import matplotlib.pyplot as plt
import numpy as np
import pyproj as pp
import rasterio
import fiona
from shapely.geometry import shape, Point
from shapely.ops import transform
from shapely.prepared import prep
from obspy.segy.core import readSEGY
import xlrd

from las import LASReader
from plot import plot
from lithology.lithology import intervals_from_las3_string
from sgy2shp import sgy2shp
import utils


class ContainerError(Exception):
    pass


class BaseContainer(object):
    """
    Holds some basic information that we want in every object. Does not
    contain any data or pointers, only the transect parameters. Maybe
    eventually we can abstract some of the methods here too.
    """
    def __init__(self, params):

        for k, v in params.items():
            if k in self.__dict__:
                msg = "Params must have unique names."
                raise ContainerError(msg)
            setattr(self, k, v)

        # The x extent will be updated at plot time.
        # TODO Make this less cryptic, or use two objects.
        self.extents = [0, 0, self.depth[1], self.depth[0]]

    def reset_data(self):
        self.data = []
        self.coords = []
        if isinstance(self, LogContainer):
            self.log_lookup = {}

    def reset_all(self):
        self.reset_data()
        self.lookup = {}


class TransectContainer(BaseContainer):
    """
    Main driver class for generating transect plots. Builds and
    accesses all other plot containers

    Args:
        params (dict): Directory containing shape files of the transects.
        data (dict): Contains the shape files for the various sub-containers.
        layers (dict): Contains the shape files for the map.
    """
    def __init__(self, params, layers, data):

        print "\n+++++++++++++++++++++++++++++++++\nInitializing"

        super(TransectContainer, self).__init__(params)

        print "Starting {0}, id {1}".format(self.title, self.id)

        # Set up 'data' — the transect line — from shapefile.
        self.data = None
        with fiona.open(data['transect_file']) as c:
            for line in c:
                if line['properties']['id'] == self.id:
                    self.data = shape(line['geometry'])
        if not self.data:
            print "No transect with ID", self.id

        # self.data.length holds the length of the transect in metres
        # But often we want ints, and sometimes the number of samples.
        # This will give us a nice linspace. Put them in params.
        # TODO: Allow decimation?
        params['length'] = self.length = int(np.floor(self.data.length))
        params['nsamples'] = self.nsamples = self.length + 1
        params['linspace'] = self.linspace = np.linspace(0, self.length, self.nsamples)

        self.tops_file = data['tops_file']
        self.log = LogContainer(data['well_dir'], params)
        velocity_model = VelocityContainer(data['velocity_file'], self.log, params)
        self.seismic = SeismicContainer(data['seismic_dir'],
                                        velocity_model,
                                        params)
        self.elevation = ElevationContainer(data['elevation_file'], params)
        self.bedrock = BedrockContainer(data['bedrock_dir'], params)
        self.striplog = StriplogContainer(data['striplog_dir'], params)
        self.potfield = PotfieldContainer(data['potfield_dir'], params)
        self.locmap = LocmapContainer(layers, params)

    def plot(self):
        """
        Generates plot for the transect.
        """

        # Set the extent to the length? Or keep them all the same?
        self.extents[0] = 0
        self.extents[1] = self.data.length

        print "\n+++++++++++++++++++++++++++++++++\nUpdating"
        # update the containers
        self.seismic.update(self.data)
        self.log.update(self.data)
        self.elevation.update(self.data)
        self.bedrock.update(self.data)
        self.striplog.update(self.data)
        self.potfield.update(self.data)
        self.locmap.update(self.data)

        print "\n+++++++++++++++++++++++++++++++++\nPlotting"
        plot(self)
        plt.show()


class LocmapContainer(BaseContainer):
    """
    Class for building and plotting a location map.

    Args:
        layers (dict): The relative paths to files to use as layers.
        params (dict): The parameters, as specified in the config.

    Example:
        >>> layers = ['geology.tif', 'roads.shp']
        >>> params = {'domain': 'depth', 'buffer': 300, ... }
        >>> mc = LocmapContainer(layers, params)
    """
    def __init__(self, layers, params):

        # First generate the parent object.
        super(LocmapContainer, self).__init__(params)

        # contents is an OrderedDict of vector or raster layers to plot.
        # The first item will be the top layer.
        self.layers = layers
        self.mid = None

    def update(self, transect):
        """
        Updates the container data to a profile that intersect the
        transect line.

        Returns nothing. Sets attributes as a side effect.

        Args:
            transect (LineString): A transect line.
        """
        print "Updating map container"

        pad = 5000      # m ... interior padding for map box
        aspect = 12/3.  # I was expecting it to be 8:3.

        bounds = transect.bounds
        llx, lly, urx, ury = bounds
        w, h = urx - llx, ury - lly

        # Guarantees the map will have correct aspect
        x_adj, y_adj = 0, 0
        if h > (w/aspect):
            x_adj = ((aspect*h) - w) / 2.  # Aspect is hard-coded in uberplot
        else:                              # TODO Fix that!
            y_adj = ((w/aspect) - h) / 2.

        utm_nad83 = pp.Proj("+init=EPSG:26920")
        ll_nad83 = pp.Proj("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs")
        utm2lola = partial(pp.transform, utm_nad83, ll_nad83)

        ll = transform(utm2lola, Point(llx-pad-x_adj, lly-pad-y_adj))
        ur = transform(utm2lola, Point(urx+pad+x_adj, ury+pad+y_adj))
        self.ll, self.ur = ll, ur
        self.mid = Point(ll.x + 0.5*(ur.x-ll.x), ll.y + 0.5*(ur.y - ll.y))

        # the goal is to send data to uberplot to plot
        # we just want to find / compute data here
        # we need to distinguish between shp and raster

        for layer, fname in self.layers.items():
            print layer, fname
            items = []
            if fname.endswith(".shp"):
                with fiona.open(fname) as c:  # collection
                    for s in c:
                        items.append(shape(s['geometry']))
            setattr(self, layer, items)


class ElevationContainer(BaseContainer):
    """
    Container for managing and plotting elevation data.

    Args:
        elevation_file: Raster file containing the elevation
            profile of the area of interest.
    """
    def __init__(self, elevation_file, params):

        super(ElevationContainer, self).__init__(params)

        self.all_data = []
        self.all_coords = []

        # transect coords need to be upsampled
        # self.nsamples = 100

        self.data = np.zeros(self.nsamples)
        self.coords = np.zeros(self.nsamples)

        decimate = 1

        with rasterio.drivers(CPL_DEBUG=True):
            with rasterio.open(elevation_file) as src:
                self.all_data = src.read()[0,
                                           0:-1:decimate,
                                           0:-1:decimate
                                           ]

                # Get as lat/lon using the affine coordinate transform
                # TODO - do this with proj? Why go via lo/la?
                la = np.arange(self.all_data.shape[1])
                lo = np.arange(self.all_data.shape[0])
                lat = la * src.affine[0]*decimate + src.affine[2]
                lon = lo * src.affine[4]*decimate + src.affine[5]
                wgs_grid = np.meshgrid(lat, lon)

                # TODO move to config
                ll_wgs84 = pp.Proj("+init=EPSG:4269")
                utm_nad83 = pp.Proj("+init=EPSG:26920")

                self.all_coords = pp.transform(ll_wgs84,
                                               utm_nad83,
                                               wgs_grid[0],
                                               wgs_grid[1]
                                               )

    def update(self, transect):
        """
        Updates the container data to a profile that intersect the
        transect line.

        Returns nothing. Sets attributes as a side effect.

        Args:
            transect (LineString): A transect line.
        """
        print "Updating", type(self)

        # space = np.linspace(0, transect.length, self.nsamples)
        for i, n in enumerate(self.linspace):

            # interpolate along the transect
            x, y = transect.interpolate(n).xy

            # Get the closest elevation points
            xi = np.abs(self.all_coords[0][0, :] - x).argmin()
            yi = np.abs(self.all_coords[1][:, 0] - y).argmin()

            self.data[i] = self.all_data[yi, xi]

            # add the distance to the coordinates
            self.coords[i] = n


class BedrockContainer(BaseContainer):
    """
    Contains a geological map or similar basic geology shapes.
    """
    def __init__(self, bedrock_dir, params):

        super(BedrockContainer, self).__init__(params)

        self.reset_all()

        # Read in all shape files
        for f in utils.listdir(bedrock_dir, '\\.shp$'):
            for line in fiona.open(f):
                self.lookup[(shape(line['geometry']))] = line['properties']

    def update(self, transect):
        """
        Updates the container data to a profile that intersect the
        transect line.

        Returns nothing. Sets attributes as a side effect.

        Args:
            transect (LineString): A transect line.
        """
        print "Updating bedrock container"

        self.reset_data()

        points = [Point(xy) for xy in transect.coords]
        for polygon, properties in self.lookup.items():
            for p in filter(polygon.contains, points):
                self.data.append(properties)
                self.coords.append(transect.project(p))

        # Sort the data to be in order
        idx = sorted(range(len(self.coords)),
                     key=lambda k: self.coords[k])

        self.data = [self.data[i] for i in idx]
        self.coords = np.array([self.coords[i] for i in idx])


class SeismicContainer(BaseContainer):
    """
    Class for reading and plotting seismic data.

    Args:
        seis_dir (str): Input directory containing seismic SEGY files.
        params (dict): The parameters, as specified in the config.

    Example:
        >>> seis_dir = 'seismic/segy_files/'
        >>> params = {'domain': 'depth', 'buffer': 300, ... }
        >>> seis = seismicContainer(seis_dir, params)
    """
    def __init__(self, seis_dir, velocity, params):

        # First generate the parent object.
        super(SeismicContainer, self).__init__(params)

        # This creates a (hidden) shapefile for each seismic line,
        # then steps over them to read their positions and meta.
        # TODO Simplify this... maybe don't even write the files.
        sgy2shp(seis_dir, seis_dir)

        self.reset_all()
        self.velocity = velocity

        for f in utils.listdir(seis_dir, '\\..+\\.shp$'):
            with fiona.open(f, "r") as traces:
                for trace in traces:
                    self.lookup[shape(trace["geometry"])] = trace["properties"]

    def update(self, transect):
        """
        Updates the container data to a profile that intersect the
        transect line.

        Returns nothing. Sets attributes as a side effect.

        Args:
            transect (LineString): A transect line.
        """
        print "Updating", type(self)

        self.z = np.linspace(self.depth[0], self.depth[1], 1000)

        self.reset_data()

        # Preprocessing
        prepared = prep(transect.buffer(self.settings['buffer']))

        # Get intersecting points
        points = filter(prepared.contains, self.lookup.keys())

        # Lookup for grouping traces into segy files
        file_lookup = {}
        for point in points:
            meta = self.lookup[point]
            print ".",
            f = meta["segyfile"]
            if f in file_lookup:
                proj_d = transect.project(point)
                if proj_d:
                    file_lookup[f]["pos"].append(proj_d)
                    file_lookup[f]["trace"].append(meta["trace"])
                    file_lookup[f]["point"].append(point)
            else:
                file_lookup[f] = {}
                file_lookup[f]["trace"] = [meta["trace"]]
                file_lookup[f]["pos"] = [transect.project(point)]
                file_lookup[f]["point"] = [point]

        # Read in the chunks from the segy file
        for segyfile in file_lookup.keys():
            print segyfile
            segy = readSEGY(segyfile, unpack_trace_headers=True)
            traces = file_lookup[segyfile]["trace"]
            coords = file_lookup[segyfile]["pos"]
            points = file_lookup[segyfile]["point"]

            # Get the sort order.
            idx = sorted(range(len(traces)), key=lambda k: traces[k])
            idx = filter(None, idx)

            coords = np.array([coords[i] for i in idx])
            data = np.array([self.velocity.time2depth(
                             segy.traces[traces[i]].data,
                             segy.traces[traces[i]].stats["sampling_rate"],
                             points[i], self.z)
                             for i in idx])

            self.data.append(np.transpose(data))
            self.coords.append(coords)

        print file_lookup.keys()


class PotfieldContainer(BaseContainer):
    """
    Contains random data for placeholders.
    """
    def __init__(self, potfield_dir, params):
        super(PotfieldContainer, self).__init__(params)

        self.all_data = {}
        self.all_coords = {}

        self.data = {}
        self.coords = {}

        decimate = 1

        exts = self.settings['raster_extensions']
        regex = ".+\\." + "$|.+\\.".join(exts) + "$"

        for f in utils.listdir(potfield_dir, regex):
            name, ext = os.path.splitext(os.path.basename(f))
            with rasterio.drivers(CPL_DEBUG=True):
                with rasterio.open(f) as src:
                    self.all_data[name] = src.read()[0,
                                                     0:-1:decimate,
                                                     0:-1:decimate
                                                     ]

                    # Get as lat/lon using the affine coordinate transform
                    # TODO - do this with proj? Why go via lo/la?
                    la = np.arange(self.all_data[name].shape[1])
                    lo = np.arange(self.all_data[name].shape[0])
                    lat = la * src.affine[0]*decimate + src.affine[2]
                    lon = lo * src.affine[4]*decimate + src.affine[5]
                    wgs_grid = np.meshgrid(lat, lon)

                    # TODO move to config
                    ll_wgs84 = pp.Proj("+init=EPSG:4269")
                    utm_nad83 = pp.Proj("+init=EPSG:26920")

                    self.all_coords[name] = pp.transform(ll_wgs84,
                                                         utm_nad83,
                                                         wgs_grid[0],
                                                         wgs_grid[1]
                                                         )

        print "Found rasters", self.all_data.keys()

    def update(self, transect):
        """
        Updates the container data to a profile that intersect the
        transect line.

        Returns nothing. Sets attributes as a side effect.

        Args:
            transect (LineString): A transect line.
        """
        print "Updating potfield container"

        for k, v in self.all_data.items():
            self.coords[k] = self.linspace
            self.data[k] = np.zeros_like(self.linspace)
            for i in self.linspace:
                i = int(i)
                x, y = transect.interpolate(i).xy
                xi = np.abs(self.all_coords[k][0][0, :] - x).argmin()
                yi = np.abs(self.all_coords[k][1][:, 0] - y).argmin()

                self.data[k][i] = v[yi, xi]


class LogContainer(BaseContainer):
    """
    Container for managing and plotting LAS data.

    Args:
        las_dir (str): Directory shape files of LAS headers.

    Example:
        >>> las_dir = 'wells/logs/'
        >>> params = {'domain': 'depth', 'buffer': 300, ... }
        >>> lc = LogContainer(las_dir, params)
    """
    def __init__(self, las_dir, params):

        # First generate the parent object.
        super(LogContainer, self).__init__(params)

        self.reset_all()
        self.log_lookup = {}  # look up for log ids

        for shp in utils.walk(las_dir, '\\.shp$'):
            with fiona.open(shp, "r") as logs:
                for log in logs:
                    self.lookup[shape(log["geometry"])] = log["properties"]['name']

    def update(self, transect):
        """
        Updates the container data to a profile that intersect the
        transect line.

        Returns nothing. Sets attributes as a side effect.

        Args:
            transect (LineString): A transect line.
        """
        print "Updating log container",

        # Preprocess
        prepared = prep(transect.buffer(self.settings['buffer']))

        # Get the intersecting points
        points = filter(prepared.contains, self.lookup.keys())

        self.reset_data()

        for point in points:
            name = self.lookup[point]
            # name = meta["name"]

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

        print "\n",

    def get(self, log_id):
        """
        Returns data corresponding to log_id
        """
        return self.log_lookup.get(log_id)

    def get_point(self, log_id):

        ids = self.lookup.values()
        points = self.lookup.keys()
        index = ids.index(log_id)

        if index:
            return points[index]
        else:
            return None


class StriplogContainer(LogContainer):

    def __init__(self, striplog_dir, params):

        # First generate the parent object.
        super(LogContainer, self).__init__(params)

        self.reset_all()

        pass

    def update(self, transect):
        # Preprocess
        prepared = prep(transect.buffer(self.settings['buffer']))

        # Get the intersecting points
        points = filter(prepared.contains, self.lookup.keys())

        self.reset_data()

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


class VelocityContainer(BaseContainer):
    """
    Container for managing velocity profiles.

    Args:
    vel_xls (str): Velocity depth lookup tables in xls format.
        Each sheet name should correspond to an existing
        log, which will contain the position information.
    las_container (LogContainer): Log container used to get position
        information.
    """

    def __init__(self, vel_xls, las_container, params):
        super(VelocityContainer, self).__init__(params)
        book = xlrd.open_workbook(vel_xls)
        self.lookup = {}
        for sheet in book.sheets():
            name = sheet.name
            point = las_container.get_point(str(name))
            if point:
                twt = np.array([i.value for i in sheet.col(1)[1:]])
                depth = np.array([i.value for i in sheet.col(2)[1:]])
                self.lookup[point] = [twt, depth]

    def get_profile(self, point):
        """
        Returns the velocity profile closest to the point

        Args:
        point (Point): Point object of the location of the desired
            velocity profile.

        Returns:
            (array, 1d) of the velocity profile closest to input point.
        """
        vel_points = self.lookup.keys()
        distance = [(np.array(p)[0] - np.array(point)[0])**2.0 +
                    (np.array(p)[0] - np.array(point)[1])**2.0 for
                    p in vel_points]
        profile = self.lookup[vel_points[np.argmin(distance)]]
        return profile

    def time2depth(self, trace, sample_rate, point, z):
        """
        Converts a seismic trace from time to depth.

        Args:
            trace (array, 1d): A 1D numpy array.
            sample_rate (float): The sample rate of the of the data [samp/sec].
            point (Point): Point object corresponding to the trace location.
        """
        profile = self.get_profile(point)
        trace_time = np.arange(0, trace.size)/sample_rate
        # Interpolate lookup table to be more samples.
        time = profile[0]
        depth = profile[1]
        vrms = 2 * depth / (time+0.000001)
        # Fix the divide by zero.
        vrms[0] = vrms[1]
        vrms = np.interp(trace_time, time, vrms)
        depth = vrms * trace_time
        return np.interp(z, depth, trace, left=0, right=0)
