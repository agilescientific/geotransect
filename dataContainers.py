from obspy.segy.core import readSEGY
from las import LASReader

import fnmatch

import os
from itertools import product

from matplotlib import pyplot as plt, gridspec
import numpy as np

import rasterio
from fiona import collection
from shapely.geometry import shape, mapping
from shapely.ops import unary_union
from shapely.prepared import prep


import pyproj as pp

class transectContainer():

    def __init__(self, transect_dir, seismic_shape,
                 elevation_raster):

        self.seismic = seismicContainer(seismic_shape)
        self.elevation = elevationContainer(elevation_raster)
        
        self.data = []
        self.buffer = 10000000000.0

        for f in os.listdir(transect_dir):

            if not f.endswith(".shp"):
                continue
            
            for line in collection(os.path.join(transect_dir,f)):
            
                self.data.append(shape(line['geometry']))

            
        

    def plot(self):

        for transect in self.data:

            fig = plt.figure()
            gs = gridspec.GridSpec(12, 12)
            
            self.seismic.update(transect.buffer(self.buffer))
                                
            self.seismic.plot(fig, gs[6:,6:])

            plt.show()
 

    
class seismicContainer():
    """
    Class for reading and plotting seismic data.

    USAGE:

    seis = seismicContainer(seis_dir)

    @param seis_dir: Input directory containing seismic shapefiles.
    The shapefiles contain points corresponding to UTM trace locations
    with properties segyfile, trace.
    """
    def __init__(self, seis_dir):

        self.positions = []
        self.metadata = []

        self.data = []
        
        for root, dirs, files in os.walk(seis_dir):

            try:
                shapefile = \
                  os.path.join(root,fnmatch.filter(files,'*.shp')[0])
            except:
                continue
            
            with collection(shapefile, "r") as traces:

                for trace in traces:

                    self.positions.append(shape(trace["geometry"]))
                    self.metadata.append(trace["properties"])
                    
                    
    def update(self, transect):
        """
        Updates the data to data that intersects the transect line.

        @param transect: A transect line as a buffered shapely object.
        """

        lookup_dict = {}
        prepared = prep(transect)

        for point, meta in zip(self.positions, self.metadata):

            if prepared.contains(point):

                if meta["segyfile"] in lookup_dict:
                    lookup_dict[meta["segyfile"]]["trace"].append(meta["trace"])

                    ## TODO Project onto the transect ##
                    lookup_dict[meta["segyfile"]]["pos"].append(point.coords[0])
                    
                else:
                    lookup_dict[meta["segyfile"]] = {}
                    lookup_dict[meta["segyfile"]]["trace"] = \
                      [meta["trace"]]
                       
                    lookup_dict[meta["segyfile"]]["pos"] = \
                      [point.coords[:]]
                    
        

        # Read in the chunks from the segy file
        for segyfile in lookup_dict.keys():

            segy = readSEGY(segyfile)

            self.data.append([trace.data for trace in segy.traces])

        
    def plot(self, fig, axis):

        fig.add_subplot(axis)
        plt.imshow(np.array(self.data)[0,:,:], aspect='auto')
        plt.axis('off')

class lasContainer():

    def __init__(self, las_file):

        self.data = LASReader(las_file, null_subs=np.nan)

    def plot(self, fig, axis, log):

        fig.add_subplot(axis)

        plt.plot(self.data.data[log], self.data.data['DEPT'])
        plt.axis('off')

    def feature_plot(self, fig, axis, log):


        fig.add_subplot(axis)

        plt.plot(self.data.data[log], self.data.data['DEPT'])
        
        

class elevationContainer():

    def __init__(self, elevation_file):

        with rasterio.drivers(CPL_DEBUG=True):
            with rasterio.open(elevation_file) as src:
                self.data = src.read()[0,:,:]

                # Get as lat/lon
                xlat = np.arange(self.data.shape[0]) * src.affine[0]\
                  + src.affine[2]
                ylat = np.arange(self.data.shape[1]) * src.affine[1]\
                  + src.affine[2]

                wgs_grid = np.meshgrid(xlat,ylat)

                ll_wgs84 = pp.Proj("+init=EPSG:4326")
                utm_nad83 = pp.Proj("+proj=utm +zone=20T,"+
                                    "+north +datum=NAD83 +units=m +"+
                                    "no_defs")

                self.coords = pp.transform(ll_wgs84, utm_nad83,
                                           wgs_grid[0], wgs_grid[1])
                
                                
                
    def plot(self, fig, axis):

        return

