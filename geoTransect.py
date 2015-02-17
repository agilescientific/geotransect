# import modules
from osgeo import gdal
from osgeo import ogr

from matplotlib import pyplot as plt, gridspec

from obspy.segy.core import readSEGY

from dataContainers import *

## These will be read from a config file .... ##
transect_file = 'data/transect_shp/transect_test.shp'
segy_file = 'data/seismic/NS00_Wind-01_PSTM.sgy'
las_file = 'data/wells/P-111/wireline_log/P-111_out.LAS'
elevation_file = 'data/elevation_raster/NS_DEM_ascii.asc'

# make the data containers
trans = transectContainer(transect_file)
seis = seismicContainer(segy_file)
las = lasContainer(las_file)
elevation = elevationContainer(elevation_file)

fig = plt.figure()
# initialize the plot
gs = gridspec.GridSpec(12, 12)

# plot the elements
seis.plot(fig, gs[2:9,4:])
las.plot(fig, gs[2:9,4:], "GR")
las.feature_plot(fig, gs[:,0:3], "GR")

elevation.plot(fig, gs[0:2,4:], 1,4000,1,10000)


plt.show()
