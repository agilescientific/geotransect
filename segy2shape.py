
import fnmatch
from obspy.segy.core import readSEGY
import fiona
from fiona import collection, crs
from shapely.geometry import Point, mapping
import os
import sys

def segy2shape(input_dir, output_dir):
    """
    Extracts trace location from SEGY files and saves it in a 
    shape file. A shape file is generated for each SEGY file.

    @param input_dir: Directory containing SEGY files
    @param output_dir: Directory to save shape files
    """

    for segyfile in os.listdir(input_dir):

        # setup the output file
        output_filebase,ext = os.path.splitext(os.path.basename(segyfile))

        if ext not in ['.SEGY', '.segy', '.SGY', '.sgy']:
            print ext
            continue
        
        outfile = os.path.join(output_dir, output_filebase)
        
        # Read in the headers
        segy = readSEGY(os.path.join(input_dir,segyfile),
                        headonly=True,
                        unpack_trace_header=True)
        schema = { 'geometry': 'Point',
                   'properties': {'segyfile': 'str', 'trace':'int' }}


        with fiona.open(outfile, 'w', driver='ESRI Shapefile',
                        crs=crs.from_epsg(26920),
                        schema=schema) as output:

            for i, trace in enumerate(segy):

                header = trace.stats.segy.trace_header
                coord_scalar = header.scalar_to_be_applied_to_all_coordinates
                if coord_scalar == -100:
                    coord_gain = 0.01
                elif coord_scalar == -10:
                    coord_gain = 0.1
                else: 
                    coord_gain = 1.0

                # UTM coordinates
                print header.source_coordinate_x, header.source_coordinate_y
                x = float(header.source_coordinate_x)*coord_gain
                y = float(header.source_coordinate_y)*coord_gain

                if x > 9e5 or y > 55e6:
                    if x > 9e6 or y > 55e7:
                        print ('found weird geometries: dividing by 100!')
                        x = x / 100.0
                        y = y / 100.0 
                    else: 
                        print ('found weird geometries: dividing by 10!')
                        x = x / 10.0
                        y = y / 10.0 
                else:
                    pass
                p = Point(x,y)
                output.write({
                'properties': {
                    'segyfile': os.path.join(input_dir,segyfile),
                    'trace': i},
                'geometry': mapping(p)})

        
def main():
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]

    segy2shape(input_dir, output_dir)

if __name__ == "__main__":
    main()

