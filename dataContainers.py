from obspy.segy.core import readSEGY
from las import LASReader
from matplotlib import pyplot as plt
import numpy as np

import rasterio
from fiona import collection
from shapely.geometry import shape, mapping
from shapely.ops import unary_union
from shapely.prepared import prep

class transectContainer():

    def __init__(self, transect_file, seismic_shape):

        self.seismic = seismicContainer(seismic_shape)
        self.data = []
        self.buffer = 1000.0
        
        for line in collection(transect_file):
            
            self.data.append(shape(line['geometry']))

            
        

    def plot(self):

        for transect in self.data:
            print len(self.seismic.update(transect.buffer(self.buffer)))
 

    
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

        for point, meta  in zip(self.positions, self.metadata):

            if prepared.contains(point):

                if meta["segyfile"] in lookup_dict:
                    lookup_dict[meta["segyfile"]]["trace"].append(meta["trace"])

                    ## TODO Project onto the transect ##
                    lookup_dict[meta["segyfile"]]["pos"].append(point.coords[:])
                    
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
        plt.imshow(self.data.traces, aspect='auto')
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
                self.res = src.res
                
    def plot(self, fig, axis, x1,x2,y1,y2):

        fig.add_subplot(axis)
        length = int(np.hypot(x2-x1, y2-y1))
        x, y = np.linspace(x1, x2, length),np.linspace(y1, y2, length)

        zi = self.data[x.astype(np.int), y.astype(np.int)]

        plt.plot(zi)
        plt.axis('off')
