# geotransect
Tool for visualizing subsurface data along 2D transects.

[![Codacy Badge](https://www.codacy.com/project/badge/83d94f2a212c46b291c9e8c72ba6ed3f)](https://www.codacy.com/app/matt/geotransect_2)

## Creating an envrironment
Best thing is to install Continuum's Anaconda. Then, from this directory:

    conda create -n geotransect --file conda_requirements.txt
    source activate geotransect
    pip install -r pip_requirements.txt

    python transect.py -c my_config.yaml

## Dependencies
Hopefully the evironment 'just worked' but if not, you might have to install dependencies by hand...

Some of these dependencies have dependencies of their own, some of which are substantial and need compiling. If you're on a Mac, I recommend using Homebrew for the tricky things like gdal and PROJ4.

- [pyproj](https://pypi.python.org/pypi/pyproj)
- [fiona](https://github.com/sgillies/fiona): high-level shapefile read/write, wraps GDAL
- [rasterio](https://github.com/sgillies/rasterio): high-level raster read/write, wraps GDAL
- [shapely](https://pypi.python.org/pypi/Shapely) probably
- [pyYAML](https://pypi.python.org/pypi/PyYAML)
- [NumPy](https://github.com/numpy/numpy) of course
- [SciPy](http://scipy.org/) for some interpolation functions
- [matplotlib](http://matplotlib.org/) for plotting
- [matplotlib Basemap](https://github.com/matplotlib/basemap) for maps
- [Pillow](http://pillow.readthedocs.org/installation.html)
- [obspy](https://github.com/obspy/obspy)
- [agilegeo](https://github.com/agile-geoscience/agilegeo/tree/develop), for depth conversion (need **develop** branch at the moment)
- [striplog](https://github.com/agile-geoscience/striplog)
