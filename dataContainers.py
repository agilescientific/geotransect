from obspy.segy.core import readSEGY
from las import LASReader

import fnmatch

import os
from itertools import product

from matplotlib import pyplot as plt, gridspec
import numpy as np
from scipy.interpolate import griddata

import rasterio
from fiona import collection
from shapely.geometry import shape, mapping, LineString
from shapely.ops import unary_union
from shapely.prepared import prep


import pyproj as pp

class transectContainer():

    def __init__(self, transect_dir, seismic_shape,
                 elevation_raster, las_file, extents):

        self.extents = extents
        self.seismic = seismicContainer(seismic_shape)
        self.elevation = elevationContainer(elevation_raster)
        self.log = lasContainer(las_file)
        
        self.data = []
   

        for f in os.listdir(transect_dir):

            if not f.endswith(".shp"):
                continue
            
            for line in collection(os.path.join(transect_dir,f)):
            
                self.data.append(shape(line['geometry']))

            
        

    def plot(self):

        for transect in self.data:

            self.extents[1] = transect.length

            fig = plt.figure()
            gs = gridspec.GridSpec(12, 12)
            
            self.seismic.update(transect)
            fig.add_subplot(gs[2:10,4:])
            self.seismic.plot(self.extents)


            self.log.update(transect)
            fig.add_subplot(gs[2:10,6:])
            self.log.plot(self.extents,"GR")
            fig.add_subplot(gs[2:12,0:3])
            self.log.feature_plot("GR")
           
            self.elevation.update(transect)
            fig.add_subplot(gs[0:2,4:])
            self.elevation.plot(self.extents)

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
    def plot(self, extents):


        # Loop through each seismic line
        for coords, data in zip(self.coords, self.data):

            z0 = 0
            depth = 5000
            plt.imshow(data,
                       extent=[np.amin(coords),
                               np.amax(coords),
                               depth, z0],
                               aspect="auto",
                       cmap="Greys")
            plt.axis(extents)

    

    def __init__(self, seis_dir):

        self.lookup = {}
        self.data = []
        self.coords = []
        self.buffer = 300.0
        
        for f in fnmatch.filter(os.listdir(seis_dir), '.*.shp'):

            shapefile = os.path.join(seis_dir, f)
        
            with collection(shapefile, "r") as traces:

                for trace in traces:

                    self.lookup[shape(trace["geometry"])] =\
                                trace["properties"]
                    
                    
    def update(self, transect):
        """
        Updates the data to data that intersects the transect line.

        @param transect: A transect line as a buffered shapely object.
        """
        
        prepared = prep(transect.buffer(self.buffer))

        points = filter(prepared.contains, self.lookup.keys())

        self.data = []
        self.coords = []

        file_lookup = {}
        for point in points:

            meta = self.lookup[point]

      
            if(meta["segyfile"] in file_lookup):

                ## project onto the transect
                proj_d = transect.project(point)

                if proj_d:
                    file_lookup[meta["segyfile"]]["pos"].append(proj_d)
                    file_lookup[meta["segyfile"]]["trace"].append(meta["trace"])               
                    
            else:
                
                file_lookup[meta["segyfile"]] = {}
                file_lookup[meta["segyfile"]]["trace"] = \
                  [meta["trace"]]

                proj_d = transect.project(point)
                file_lookup[meta["segyfile"]]["pos"] = [proj_d]
                    
    
        # Read in the chunks from the segy file
        for segyfile in file_lookup.keys():

            segy = readSEGY(segyfile)
            traces = file_lookup[segyfile]["trace"]
            coords = file_lookup[segyfile]["pos"]

            # sort the traces
            idx = sorted(range(len(traces)), key=lambda k: traces[k])
            
            self.data.append(np.transpose(np.array(
                [segy.traces[traces[i]].data for i in idx])))
            
            self.coords.append(np.array([coords[i] for i in idx]))
            
        
class lasContainer():

    def __init__(self, las_dir):

        # Read in the shape file
        self.lookup = {}
        self.data = []
        self.coords = []
        self.buffer = 300 # m

        for root, dirs, files in os.walk(las_dir):

            try:
                shapefile = \
                  os.path.join(root,fnmatch.filter(files,'*.shp')[0])
            except:
                continue


            with collection(shapefile, "r") as logs:
               
                for log in logs:

                    # add to the lookup table
                    self.lookup[shape(log["geometry"])] =\
                        log["properties"]


    def update(self, transect):

        prepared = prep(transect.buffer(self.buffer))

        points = filter(prepared.contains, self.lookup.keys())

        self.data = []
        self.coords = []

        for point in points:

            meta = self.lookup[point]

            name = meta["name"]
            filename = os.path.join('data', 'wells', name,
                                    'wireline_log', name +
                                    '_out.LAS')

            if not os.path.exists(filename): continue
            self.data.append(LASReader(filename, null_subs=np.nan))
            self.coords.append(transect.project(point))

            
    def plot(self, extents, log):

        for las, pos in zip(self.data, self.coords):
            data = np.nan_to_num(las.data[log])
            data /= np.amax(data)
            data *= .1*(extents[1] - extents[0])
            data += pos
        
            plt.plot(data, las.data['DEPT'])
        
            plt.xlim((extents[0], extents[1]))
            plt.ylim((extents[2], extents[3]))
        
            plt.axis("off")
        plt.axis('off')
            

    def feature_plot(self, log):

        if self.data:
            plt.plot(self.data[0].data[log],
                     self.data[0].data['DEPT'])
        
        

class elevationContainer():

    def __init__(self, elevation_file):

        self.elevation_profile = []
        self.elevation_grid = []

        self.data = []
        self.coords = []

        decimate = 1
        
        with rasterio.drivers(CPL_DEBUG=True):
            with rasterio.open(elevation_file) as src:
                self.elevation_profile = src.read()[0,0:-1:decimate,
                                                    0:-1:decimate]

                # Get as lat/lon using the affine coordinate
                # transform
                lat = np.arange(self.elevation_profile.shape[1]) \
                  *src.affine[0]*decimate + src.affine[2]
                  
                lon = np.arange(self.elevation_profile.shape[0])\
                   * src.affine[4]*decimate + src.affine[5]

                wgs_grid = np.meshgrid(lat,lon)

                # TODO make sure these projections are legitimate
                ll_wgs84 = pp.Proj("+init=EPSG:4269")
                utm_nad83 = pp.Proj("+init=EPSG:26920")

                self.elevation_grid = pp.transform(ll_wgs84,
                                                   utm_nad83,
                                                   wgs_grid[0],
                                                   wgs_grid[1])
                
                                

    def update(self, transect):
        
        # transect coords need to be upsampled
        nsamples = 100

        self.coords = np.zeros(nsamples)
        self.data = np.zeros(nsamples)

        for i,n in enumerate(np.linspace(0,transect.length,
                                         nsamples)):

            # interpolate along the transect
            x,y = transect.interpolate(n).xy


            # Get the closest elevation points
            xi = np.abs(self.elevation_grid[0][0,:] - x).argmin()
            yi = np.abs(self.elevation_grid[1][:,0] - y).argmin()

         
            self.data[i] = self.elevation_profile[yi,xi]
            
            # add the distance to the coordinates
            self.coords[i] = n


    def plot(self, extents):

        plt.plot(self.coords, self.data)
        plt.xlim((extents[0], extents[1]))
        plt.axis('off')

    

        

