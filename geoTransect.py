from dataContainers import transectContainer

## These will be read from a config file .... ##
transect_file = 'data/transect_shp'
seismic_dir = 'data/seismic/NS00_Wind-01_PSTM'
las_shape = 'data/wells/'
elevation_file = 'data/elevation_raster/NS_DEM_ascii.asc'
extents = [0,50000, -10000,0]

# Initial the main object
tc = transectContainer(transect_file, seismic_dir, elevation_file,
                       las_shape, extents)


tc.plot()

