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
import time

# Optional backend specification
# import matplotlib; matplotlib.use("WX")

# Import required to avoid bug in Basemap
from mpl_toolkits.basemap import Basemap

# 3rd party
import matplotlib.pyplot as plt
import numpy as np
import pyproj as pp
import rasterio
import fiona
from shapely.geometry import shape, Point
from shapely.ops import transform
from shapely.prepared import prep
from obspy.segy.core import readSEGY

# Other Agile libraries
from agilegeo.avo import time_to_depth, depth_to_time
from striplog import Well, Lexicon

from plot import plot
from sgy2shp import sgy2shp, ShapeFileExists
from notice import Notice
import utils


class ContainerError(Exception):
    pass


class BaseContainer(object):
    """
    Holds some basic information that we want in every object. Does not
    contain any data or pointers, only the transect parameters. Maybe
    eventually we can abstract some of the methods here too.

    Args:
        params (dict): The params you want to store. None are compulsory.
    """
    def __init__(self, params):

        for k, v in params.items():
            if k in self.__dict__:
                msg = "Params must have unique names."
                raise ContainerError(msg)
            setattr(self, k, v)

        # The x extent will be updated at plot time.
        # TODO Make this less cryptic, or use two objects.
        rng = getattr(self, 'range', None)
        if rng:
            self.extents = [0, 0, self.range[1], self.range[0]]

    def reset_data(self):
        """
        Set some basic objects to empty placeholders.

        No args, no return. Side-effect: sets attributes.
        """
        self.data = []
        self.coords = []
        if isinstance(self, LogContainer):
            self.log_lookup = {}

    def reset_all(self):
        """
        Called at the start of an update. Just to make
        sure everything is empty.

        No args, no return. Side-effect: sets attributes.
        """
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
    def __init__(self, **kwargs):

        Notice.title()
        Notice.hr_header("Initializing")

        params = kwargs.get('params')
        layers = kwargs.get('layers')
        potfields = kwargs.get('potfields')
        data = kwargs.get('data')
        velocity = kwargs.get('velocity')

        super(TransectContainer, self).__init__(params)

        print "Starting {}, id {}, file {}".format(self.title,
                                                   self.id,
                                                   self.config_file)

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
        params['length'] = self.length = int(np.floor(self.data.length))
        params['nsamples'] = self.nsamples = self.length + 1
        params['linspace'] = self.linspace = np.linspace(0, self.length,
                                                         self.nsamples)

        self.time = time.strftime("%Y/%m/%d %H:%M", time.localtime())
        self.tops_file = data['tops_file']
        self.log = LogContainer(data['well_dir'], params)
        self.velocity = self.__velocity_factory(velocity, params)
        self.seismic = SeismicContainer(data['seismic_dir'],
                                        self.velocity,
                                        params)
        self.horizons = HorizonContainer(data['horizon_dir'],
                                         self.velocity,
                                         params)
        self.elevation = ElevationContainer(data['elevation_file'], params)
        self.bedrock = BedrockContainer(data['bedrock_dir'], params)
        self.potfield = PotfieldContainer(potfields, params)
        self.locmap = LocmapContainer(layers, params)

    def plot(self):
        """
        Generates plot for the transect.
        """

        # Set the extent to the length? Or keep them all the same?
        self.extents[0] = 0
        self.extents[1] = self.data.length

        Notice.hr_header("Updating")
        # update the containers
        self.velocity.update(self.data)
        self.log.update(self.data)
        self.seismic.update(self.data)
        self.horizons.update(self.data)
        self.elevation.update(self.data)
        self.bedrock.update(self.data)
        self.potfield.update(self.data)
        self.locmap.update(self.data)

        Notice.hr_header("Plotting")
        plot(self)
        plt.show()

    def __velocity_factory(self, model_params, params):
        """
        Factory function for returning a VelocityContainer matching the
        given params
        """
        if model_params["type"] == "constant":
            return ConstantVelocityContainer(model_params["data"], params)
        elif model_params["type"] == "segy":
            return SegyVelocityContainer(model_params["data"], params)
        elif model_params["type"] == "simple":
            return SimpleVelocityContainer(model_params["data"], params)
        else:
            raise VelocityError("Invalid velocity type")


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
        Notice.info("Updating " + self.__class__.__name__)

        pad = self.settings['map_padding']
        aspect = 12/3.  # I was expecting it to be 8:3.

        # Calculate the map bounds and centre.
        bounds = transect.bounds
        llx, lly, urx, ury = bounds
        w, h = urx - llx, ury - lly

        x_adj, y_adj = 0, 0
        if h > (w/aspect):
            x_adj = ((aspect*h) - w) / 2.  # Aspect is hard-coded in uberplot
        else:
            y_adj = ((w/aspect) - h) / 2.

        utm_nad83 = pp.Proj("+init=EPSG:26920")
        ll_nad83 = pp.Proj("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs")
        utm2lola = partial(pp.transform, utm_nad83, ll_nad83)

        ll = transform(utm2lola, Point(llx-pad-x_adj, lly-pad-y_adj))
        ur = transform(utm2lola, Point(urx+pad+x_adj, ury+pad+y_adj))
        self.ll, self.ur = ll, ur
        self.mid = Point(ll.x + 0.5*(ur.x-ll.x), ll.y + 0.5*(ur.y - ll.y))

        # Go over the layers and collect data.
        for layer, details in self.layers.items():
            path = details['file']
            print layer, path

            # Set up convenient params dictionary for plotting function.
            params = {k: v for k, v in details.items() if k != 'file'}
            self.layers[layer]['params'] = params

            # Get a list of shapes from the file.
            shapes = []
            fname, ext = os.path.splitext(os.path.basename(path))
            if ext.strip('.').lower() in self.settings['raster_extensions']:
                # TODO: Deal with rasters.
                pass
            elif ext.strip('.').lower() == 'shp':
                with fiona.open(path) as c:
                    for s in c:
                        shapes.append(shape(s['geometry']))
                        # name = s.get('name') or s.get('id') or None
                        # data = {name: shape(s['geometry'])}
                        # setattr(self, 'data', data)
                setattr(self, layer, shapes)
            else:
                pass


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
        self.coords = self.linspace

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
        Notice.info("Updating " + self.__class__.__name__)

        # space = np.linspace(0, transect.length, self.nsamples)
        for i, n in enumerate(self.linspace):

            # interpolate along the transect
            x, y = transect.interpolate(n).xy

            # Get the closest elevation points
            xi = np.abs(self.all_coords[0][0, :] - x).argmin()
            yi = np.abs(self.all_coords[1][:, 0] - y).argmin()

            self.data[i] = self.all_data[yi, xi]


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
        Notice.info("Updating " + self.__class__.__name__)

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


class SegyContainer(BaseContainer):
    """
    Superclass for handling Segy data.

    Args:
        segy_dir (str): Input directory containing Segy files.
        params (dict): The parameters, as specified in the config

        Example:
        >>> seis_dir = 'seismic/segy_files/'
        >>> params = {'domain': 'depth', 'buffer': 300, ... }
        >>> segy = SegyContainer(seis_dir, params)
    """

    def __init__(self, segy_dir, params):

        super(SegyContainer, self).__init__(params)

        try:
            # This creates a (hidden) shapefile for each seismic line,
            # then steps over them to read their positions and meta.
            # TODO Simplify this... maybe don't even write the files.
            sgy2shp(segy_dir, segy_dir)
        except ShapeFileExists:
            pass

        self.reset_all()

        for f in utils.listdir(segy_dir, '\\..+\\.shp$'):
            with fiona.open(f, "r") as traces:
                for trace in traces:
                    self.lookup[shape(trace["geometry"])] = trace["properties"]

    def update(self, transect, flat=False):
        """
        Updates the container data to a profile that intersect the
        transect line.

        Returns nothing. Sets attributes as a side effect.

        Args:
            transect (LineString): A transect line.
            flat (Bool): Reads data into a flat list instead of
                         sorting by files.
        """
        Notice.info("Updating " + self.__class__.__name__)

        self.reset_data()

        # Preprocessing
        prepared = prep(transect.buffer(self.settings['buffer']))

        # Get intersecting points
        points = filter(prepared.contains, self.lookup.keys())

        # Lookup for grouping traces into segy files
        count = 0
        file_lookup = {}
        for point in points:
            meta = self.lookup[point]
            count += 1
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
            print os.path.basename(segyfile)
            segy = readSEGY(segyfile, unpack_trace_headers=True)
            traces = file_lookup[segyfile]["trace"]
            coords = file_lookup[segyfile]["pos"]
            points = file_lookup[segyfile]["point"]

            # Get the sort order.
            idx = sorted(range(len(traces)), key=lambda k: traces[k])
            idx = filter(None, idx)

            coords = [coords[i] for i in idx]
            data = [segy.traces[traces[i]] for i in idx]

            if flat:
                self.data += data
                self.coords += coords
            else:
                self.data.append(data)
                self.coords.append(coords)


class SeismicContainer(SegyContainer):
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
        super(SeismicContainer, self).__init__(seis_dir, params)
        self.velocity = velocity

    def update(self, transect):
        """
        Gather the data near the transect. Depth convert if necessary.

        Returns nothing. Side effect: sets attributes.

        Args:
            transect (LineString): A shapely LineString object.
        """
        # Do the super class
        super(SeismicContainer, self).update(transect)

        data = []
        # Loop through files
        for segydata, coords in zip(self.data, self.coords):
            # Through traces
            traces = []
            for trace, coord in zip(segydata, coords):
                # Get the sample rate, in Hz.
                samp = trace.stats["sampling_rate"]

                if self.domain.lower() in ['depth', 'd', 'z']:
                    trz = self.velocity.time2depth(trace.data, coord, 1/samp,
                                                   self.dz)
                    traces.append(trz)
                    samp = 1.0/self.dz

                else:
                    traces.append(trace.data)
                    samp /= 1000.0

            struct = {"sample_interval": 1.0/samp,
                      "traces": np.transpose(np.array(traces))}
            data.append(struct)

        self.data = data


class HorizonContainer(BaseContainer):
    """
    Class for reading and plotting seismic horizons.

    Args:
        hor_dir (str): Input directory containing seismic SEGY files.
        velocity (Object): A velocity container for time-depth
                           conversion.
        params (dict): The parameters, as specified in the config.
    """

    def __init__(self, hor_dir, velocity, params):

        # First generate the parent object.
        super(HorizonContainer, self).__init__(params)

        self.velocity = velocity

        self.data = {}
        self.coords = {}
        self.lookup = {}

        for fname in utils.listdir(hor_dir):
            with open(fname) as f:
                samples = f.readlines()
            name = samples.pop(0).strip().strip('#')
            points = []
            for s in samples:
                line, cdp, x, y, t, surv = s.split()
                x, y = int(float(x)), int(float(y))
                t = float(t)
                points.append(Point(x, y, t))
            if self.lookup.get(name):
                self.lookup[name] += points
            else:
                self.lookup[name] = points

    def update(self, transect):
        """
        Gather the data near the transect. Depth convert if necessary.

        Returns nothing. Side effect: sets attributes.

        Args:
            transect (LineString): A shapely LineString object.
        """
        Notice.info("Updating " + self.__class__.__name__)

        b = self.settings['fine_buffer']
        prepared = prep(transect.buffer(b))

        for horizon, points in self.lookup.items():
            l = len(points)
            points = filter(prepared.contains, points)
            print horizon, len(points), "of", l, "points"
            data, coords = [], []
            for p in points:
                coords.append(transect.project(p))
                if self.domain.lower() in ['depth', 'd', 'z']:
                    "Depth converting horizon"
                    zpt = self.velocity.time2depthpt(p.z/1000, coords[-1])
                    data.append(zpt)
                else:
                    data.append(p.z)

            self.data[horizon] = np.array(data)
            self.coords[horizon] = np.array(coords)


class PotfieldContainer(BaseContainer):
    """
    Contains potential field data.
    """
    def __init__(self, potfields, params):
        super(PotfieldContainer, self).__init__(params)

        self.data = {}

        for name, pf_params in potfields.items():

            # Populate these now.
            all_data = {}
            all_coords = {}

            all_data, all_coords = self.__get_all_data(pf_params)

            c_def = self.settings['default_colour']
            cif = pf_params['colour_is_file']
            if cif:
                all_colour = self.__get_all_data(pf_params, colour=True)[0]
            else:
                all_colour = pf_params.get('colour', c_def)

            cmap = pf_params.get('cmap')
            if not cmap:
                cmap = self.settings.get('default_cmap')

            scale = pf_params.get('scale')

            payload = {'all_data': all_data,
                       'all_coords': all_coords,
                       'all_colour': all_colour,
                       'colour_is_file': cif,
                       'cmap': cmap,
                       'scale': scale}

            self.data[name] = payload

    def update(self, transect):
        """
        Updates the container data to a profile that intersect the
        transect line.

        Returns nothing. Sets attributes as a side effect.

        Args:
            transect (LineString): A transect line.
        """
        Notice.info("Updating " + self.__class__.__name__)

        for name, payload in self.data.items():
            print name
            payload['coords'] = self.linspace
            payload['data'] = np.zeros_like(self.linspace)
            if payload['colour_is_file']:
                payload['colour'] = np.zeros_like(self.linspace)
            else:
                payload['colour'] = payload['all_colour']

            for i in self.linspace:
                i = int(i)
                x, y = transect.interpolate(i).xy
                xi = np.abs(payload['all_coords'][0][0, :] - x).argmin()
                yi = np.abs(payload['all_coords'][1][:, 0] - y).argmin()

                payload['data'][i] = payload['all_data'][yi, xi]
                if payload['colour_is_file']:
                    try:
                        payload['colour'][i] = payload['all_colour'][yi, xi]
                    except IndexError:
                        payload['colour'][i] = 0

            self.data[name] = payload

    def __get_all_data(self, params, colour=False):
        """
        Just a raster data reader.

        TODO: Needs refactoring.
        """
        with rasterio.drivers(CPL_DEBUG=True):

            if colour:
                f = params['colour']
            else:
                f = params['file']

            with rasterio.open(f) as src:

                # Read and decimate the data.
                decimate = params.get('decimate', 1)
                all_data = src.read()[0,
                                      0:-1:decimate,
                                      0:-1:decimate
                                      ]

                # Clip the data to deal with outliers
                perc = params.get('clip', 100)
                if colour:
                    perc = 95
                vmin = np.percentile(all_data, 100-perc)
                vmax = np.percentile(all_data, perc)
                all_data[all_data > vmax] = vmax
                all_data[all_data < vmin] = vmin

                # Get as lat/lon using the affine coordinate transform
                # TODO - do this with proj?
                y = np.arange(all_data.shape[1])
                x = np.arange(all_data.shape[0])
                y = y * src.affine[0]*decimate + src.affine[2]
                x = x * src.affine[4]*decimate + src.affine[5]
                wgs_grid = np.meshgrid(y, x)

                utm_nad83 = pp.Proj("+init=EPSG:26920")

                all_coords = pp.transform(utm_nad83,
                                          utm_nad83,
                                          wgs_grid[0],
                                          wgs_grid[1]
                                          )
        return all_data, all_coords


class LogContainer(BaseContainer):
    """
    Container for managing and plotting LAS data.

    Args:
        well_dir (str): Directory shape files of LAS headers.

    Example:
        >>> well_dir = 'wells/logs/'
        >>> params = {'domain': 'depth', 'buffer': 300, ... }
        >>> lc = LogContainer(well_dir, params)
    """
    def __init__(self, well_dir, params):

        # First generate the parent object.
        super(LogContainer, self).__init__(params)

        self.well_dir = well_dir
        self.reset_all()

        for shp in utils.walk(well_dir, '\\.shp$'):
            with fiona.open(shp, "r") as wells:
                for well in wells:
                    shp = shape(well['geometry'])
                    self.lookup[shp] = well["properties"]['name']

    def update(self, transect):
        """
        Updates the container data to a profile that intersect the
        transect line.

        Returns nothing. Sets attributes as a side effect.

        Args:
            transect (LineString): A transect line.
        """
        Notice.info("Updating " + self.__class__.__name__)

        # Preprocess
        prepared = prep(transect.buffer(self.settings['buffer']))

        # Get the intersecting points
        points = filter(prepared.contains, self.lookup.keys())

        self.reset_data()
        self.names = []

        for point in points:
            name = self.lookup[point]
            self.names.append(name)
            print name,

            pattern = "^" + name + ".*out.las"
            for fname in utils.walk(self.well_dir, pattern):
                # This is a loop but there should only be one matching file.
                well = Well(fname, null_subs=np.nan)
                print well.curves.names
                self.data.append(well)
                self.log_lookup[name] = self.data[-1]

            if not self.log_lookup.get(name):
                print
                self.data.append(None)

            sl_name = getattr(self, 'striplog', None)
            sl = None
            if sl_name and (name == self.feature_well):
                lexicon = Lexicon.default()
                pattern = "^" + name + ".*striplog.las"
                for fname in utils.walk(self.well_dir, pattern):
                    # Load the striplog.
                    sl = Well(fname, lexicon=lexicon, null_subs=np.nan)

                    # Add it to the well
                    self.log_lookup[name].add_striplog(sl.striplog[sl_name],
                                                       sl_name)

            self.coords.append(transect.project(point))

    def get(self, log_id):
        """
        Returns data corresponding to log_id.
        """
        return self.log_lookup.get(log_id)

    # def get_point(self, log_id):
    #     ids = self.lookup.values()
    #     points = self.lookup.keys()
    #     index = ids.index(log_id)

    #     if index:
    #         return points[index]
    #     else:
    #         return None


class VelocityError(Exception):
    """
    Gerneric error class.
    """
    pass


class ConstantVelocityContainer(BaseContainer):
    """
    Class for handling constant velocity models
    """

    def __init__(self, velocity, params):
        super(ConstantVelocityContainer, self).__init__(params)
        self.velocity = velocity

    def update(self, transect):
        """
        Does nothing, velocity is constant
        """
        pass

    def depth2timept(self, d, point):
        """
        Converts a single sample from depth to time.

        Args:
            d (float): Depth in metres.
            point (float): Range along transect in m. Not used
                here, only for consistent API.

        Returns:
            Float: TWT [s]
        """
        return 2*d/self.velocity

    def time2depthpt(self, t, point):
        """
        Converts a single sample from time to depth.

        Args:
            t (float): TWT in seconds.
            point (float): Range along transect. Not used
                here, only for consistent API.

        Returns:
            Float: Depth [m].
        """
        return t*self.velocity/2.0

    def time2depth(self, data, point, dt, dz):
        """
        Converts a data array from time to depth

        Args:
            data (array, 1d): A 1D numpy array.
            point (float):  distance along transect corresponding to
            the trace location. Not actually used, as this model
            assumes a constant velocity

            dt (float): The sample interval of the input data.
            dz (float): The sample interval of the depth converted
                        data.
        """
        velocity = self.velocity
        return time_to_depth(data, np.ones(data.size)*velocity, dt, dz)

    def depth2time(self, data, point, dz, dt):
        """
        Converts a data array in the depth domain to the time
        domain.

        Args:
            data (array, 1d): A 1D numpy array.
            point (float):  distance along transect corresponding to
                the trace location. Not actually used, as this model
                assumes a constant velocity. Kept for consistency with other
                velocity methods.
            dz (float): The sample interval of the input data.
            dt (float): The sample interval of the converted
                        data.
        """
        velocity = self.velocity
        return depth_to_time(data, np.ones(data.size)*velocity, dz, dt)


class SimpleVelocityContainer(BaseContainer):
    """
    Class for handling simple velocity profiles read from text
    files. Simple example below, # are not parsed and are used as
    comments.

    coordinates, 150,20.25

    # time [s], depth [m]
    0, 0
    1000, 1000
    2000, 2000
    3000, 3000

    """
    def __init__(self, vel_dir, params):

        # Initialize the base class
        super(SimpleVelocityContainer, self).__init__(params)

        self.profiles = {}
        for profile in fnmatch.filter(os.listdir(vel_dir), '*.txt'):

            with open(os.path.join(vel_dir, profile), 'r') as f:

                # Read the location out of the file
                header = False
                loc = None
                while not header:
                    l = f.readline().lstrip()
                    if l.startswith('#'):
                        header = True
                    elif l.startswith('coordinates'):

                        x, y = l.split(',')[1:]
                        x = float(x)
                        y = float(y)

                        loc = Point((x, y))
                        header = True

                if not loc:
                    # There were no coordinates.
                    loc = Point((0, 0))

                # Build up the time depth profile
                self.profiles[loc] = []
                for line in f.readlines():

                    # Ignore comments
                    if line.strip().startswith('#'):
                        continue

                    if len(line.strip().split(',')) == 2:
                        time, depth = line.strip().split(',')
                    elif len(line.strip().split()) == 2:
                        time, depth = line.strip().split()
                    else:
                        Notice.warning("Could not read velocity file.")
                        continue

                    self.profiles[loc].append([float(time), float(depth)])

    def update(self, transect):
        """
        Gather the data near the transect. Depth convert if necessary.

        Returns nothing. Side effect: sets attributes.

        Args:
            transect (LineString): A shapely LineString object.
        """
        self.data = []
        self.coords = []

        # Grab the closest velocity profile for each point
        for transect_point in transect.coords:
            min_dist = np.Inf
            for point in self.profiles.keys():

                if((point.coords[0][0]-transect_point[0])**2 +
                   (point.coords[0][1] - transect_point[1])**2) < min_dist:

                    # update the data
                    self.data.append(np.array(self.profiles[point]))
                    self.coords.append(transect.project(point))

    def time2depthpt(self, t, point):
        """
        Converts a single sample from time to depth.

        Args:
            t (float): TWT in seconds.
            point (float): Range along transect [m].

        Returns:
            Float: Depth in m.
        """
        # Find the nearest profile
        coords = np.array(self.coords)
        closest = np.abs(coords - point).argmin()

        time = self.data[closest][:, 0]
        depth = self.data[closest][:, 1]

        return np.interp(t, time, depth)

    def depth2timept(self, d, point):
        """
        Converts a single sample from depth to time.

        args:
            d (float): depth in meters
            point (float): range along transect

        returns: twt [s]
        """
        # Find the nearest profile
        coords = np.array(self.coords)
        closest = np.abs(coords - point).argmin()

        time = self.data[closest][:, 0]
        depth = self.data[closest][:, 1]

        return np.interp(d, depth, time)

    def time2depth(self, data, point, dt, dz):
        """
        Converts a data array from time to depth

        Args:
            data (array, 1d): A 1D numpy array.
            point (float):  distance along transect corresponding to
            the trace location. Not actually used, as this model
            assumes a constant velocity

            dt (float): The sample interval of the input data.
            dz (float): The sample interval of the depth converted
                        data.
        """
        # Find the nearest profile
        coords = np.array(self.coords)
        closest = np.abs(coords - point).argmin()

        time = self.data[closest][:, 0]
        depth = self.data[closest][:, 1]

        time[0] = time[1]
        vavg = depth / time

        # resample vavg to match the data
        input_time = np.arange(data.size) * dt

        vavg = np.interp(input_time, time, vavg, vavg[0], vavg[1])

        # calculate the z axis
        z_in = vavg * input_time
        z_out = np.arange(z_in[0], z_in[-1], dz)

        # do the conversion
        return np.interp(z_out, z_in, data, data[0], data[1])

    def depth2time(self, data, point, dz, dt):
        """
        Converts a data array in the depth domain to the time
        domain.

        Args:
            data (array, 1d): A 1D numpy array.
            point (float):  distance along transect corresponding to
                the trace location. Not actually used, as this model
                assumes a constant velocity. Kept for consistency with other
                velocity methods.
            dz (float): The sample interval of the input data.
            dt (float): The sample interval of the converted data.
        """
        # Find the nearest profile
        coords = np.array(self.coords)
        closest = np.abs(coords - point).argmin()
        time = self.data[closest][:, 0]
        depth = self.data[closest][:, 1]

        time[0] = time[1]  # divide by zero hack
        vavg = depth / time
        print "vavg", vavg[:5], vavg[-5:]

        # Get an actual interval
        if np.size(dz) > 1:
            dz = dz[2] - dz[1]

        # resample vavg to match the data
        input_depth = np.arange(data.size) * dz
        print "dz...", dz
        print "input_depth", input_depth[:5], input_depth[-5:]

        vavg = np.interp(input_depth, depth, vavg, vavg[0], vavg[1])
        print "vavg now", vavg[:5], vavg[-5:]

        # calculate the linear z axis
        t_in = input_depth / vavg
        print "len tin", len(t_in)
        print "t_in...", t_in[0], t_in[1:5], t_in[-1]
        print "dt.....", dt
        t_out = np.arange(t_in[0], t_in[-1], dt)

        # do the conversion
        return np.interp(t_out, t_in, data, data[0], data[1])


class SegyVelocityContainer(SegyContainer):
    """
    Container for managing velocity profiles.

    See SegyContainer.
    """

    def __init__(self, vel_dir, params):
        super(SegyVelocityContainer, self).__init__(vel_dir, params)

    def update(self, transect):
        super(SegyVelocityContainer, self).update(transect, flat=True)

    def depth2timept(self, d, point):
        """
        Converts a single sample from depth to time.

        args:
            d (float): depth in meters
            point (float): range along transect

        returns: twt [s]
        """
        # Find the nearest profile
        coords = np.array(self.coords)
        closest = np.abs(coords - point).argmin()
        profile = self.data[closest]

        # HACK TO MAKE FAKE VELOCITIES
        velocity = np.ones(profile.data.size) * 2000.0

        samp = profile.stats["sampling_rate"]
        z = np.arange(velocity.size) * 1.0 / samp
        vrms = np.cumsum(velocity) / (np.arange(velocity.size)+1)
        vrms = np.interp(d, z, vrms)
        return 2.0*d/vrms

    def time2depthpt(self, t, point):
        """
        Converts a single sample from depth to time.

        args:
            t (float): time in [s]
            point (float): range along transect

        returns: depth [m]
        """
        # Find the nearest profile
        coords = np.array(self.coords)
        closest = np.abs(coords - point).argmin()
        profile = self.data[closest]

        # HACK TO MAKE FAKE VELOCITIES
        velocity = np.ones(profile.data.size) * 2000.0

        samp = profile.stats["sampling_rate"]
        t_profile = np.arange(velocity.size) * 1.0 / samp
        vrms = np.cumsum(velocity) / (np.arange(velocity.size)+1)
        vrms = np.interp(t, t_profile, vrms)
        return t*vrms/2.

    def depth2time(self, trace, point, dz, dt):
        """
        Converts a data array in the depth domain to the time
        domain.

        Args:
            data (array, 1d): A 1D numpy array.
            point (float):  distance along transect corresponding to
                the trace location. Not actually used, as this model
                assumes a constant velocity. Kept for consistency with other
                velocity methods.
            dz (float): The sample interval of the input data.
            dt (float): The sample interval of the converted
                        data.
        """
        distance = [(p - point)**2.0 for p in self.coords]
        idx = np.argmin(distance)
        profile = self.data[idx]
        seismic = np.array(trace)

        # HACK TO MAKE FAKE VELOCITIES
        velocity = np.ones(profile.data.size) * 2000.0

        print idx+1, 'of', len(self.coords)
        samp = profile.stats["sampling_rate"]
        vel_z = np.arange(velocity.size) * 1.0 / samp
        z = np.arange(seismic.size)*dz
        velocity = np.interp(z, vel_z, velocity, velocity[0], velocity[-1])
        t = np.arange(self.range[0], self.range[1], dt)
        data = depth_to_time(seismic, velocity, z, t, mode="cubic")
        return data

    def time2depth(self, trace, point, dt, dz):
        """
        Converts a seismic trace from time to depth.

        Args:
            trace (array, 1d): A 1D numpy array.
            point (float):  distance along transect corresponding to
                the trace location.
        """
        distance = [(p - point)**2.0 for p in self.coords]
        idx = np.argmin(distance)
        profile = self.data[idx]
        seismic = np.array(trace)

        # HACK TO MAKE FAKE VELOCITIES
        velocity = np.ones(profile.data.size) * 2000.0

        print idx+1, 'of', len(self.coords)
        samp = profile.stats["sampling_rate"]
        vel_t = np.arange(velocity.size) * 1.0 / samp
        t = np.arange(seismic.size)*dt
        velocity = np.interp(t, vel_t, velocity, velocity[0], velocity[-1])
        z = np.arange(self.range[0], self.range[1], dz)
        data = time_to_depth(seismic, velocity, t, z, mode="cubic")
        return data
