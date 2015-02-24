from dataContainers import transectContainer
import sys

## These will be read from a config file .... ##
transect_file = 'data/transect_shp'
seismic_dir = 'data/seismic/'
las_shape = 'data/wells'
elevation_file = 'data/elevation_raster/NS_DEM_ascii.asc'
bedrock_dir = 'data/bedrock_shapefile'
striplog_dir = 'data/wells'


extents = [0,40000, 2500,0]

# Initial the main object
tc = transectContainer(transect_file, seismic_dir, elevation_file,
                       las_shape, bedrock_dir,
                       striplog_dir,extents)



# Make the plots
tc.plot()


