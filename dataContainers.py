from obspy.segy.core import readSEGY
from las import LASReader
from matplotlib import pyplot as plt
import numpy as np

import rasterio
import fiona


class transectContainer():

    def __init__(self, transect_file):

        with fiona.drivers(CPL_DEBUG=True):
            with fiona.open(transect_file) as src:

                self.bounds = src.bounds

                # Assumes we have only one shape in the filter
                for i in src.filter():

                    self.coords = i["geometry"]["coordinates"]


class seismicContainer():

    def __init__(self, segy_file):

        self.data = readSEGY(segy_file, unpack_headers=True)

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
