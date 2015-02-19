from dataContainers import transectContainer

## These will be read from a config file .... ##
transect_file = 'data/transectTest/all'
seismic_dir = 'data/seismic/'
las_file = 'data/wells/P-111/wireline_log/P-111_out.LAS'
elevation_file = 'data/elevation_raster/NS_DEM_ascii.asc'
extents = [0,100000, -10000,0]

# Initial the main object
tc = transectContainer(transect_file, seismic_dir, elevation_file,
                       las_file, extents)

tc.plot()

