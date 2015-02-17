#!/usr/bin/python	

from __future__ import print_function
import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys,os 
from obspy.segy.core import readSEGY as streamReadSEGY
from obspy.segy.segy import readSEGY
from obspy.core.utcdatetime import UTCDateTime
from osgeo import gdal
from osgeo import ogr
import string
        
# find SEGY files
def file_is_segy(path):
    segy_extension = ('.SEGY','.segy','.SGY','.sgy')
    found = False
    if path.endswith(segy_extension):
        if path.find("3D") == -1: # omitting for 3D surveys
            found = True
    else:
        found = False
    return found

def remove_vintage(name):
    #find index of last underscore in name
    idx = name.rfind("_")
    return name[:idx] 

def get_vintage(name):
    idx = name.rfind("_")
    fdx = name.rfind(".")
    return name[idx+1:fdx].upper() 
    
# find unique lines
def path_to_lines(start_path):
    line_list = []
    line_paths = []
    for root, dirs, files in os.walk(start_path):
        for f in files:
            #print (f)
            if file_is_segy(f):
                segy_path = root + '/' + f
                linename = remove_vintage(f)
                if linename in line_list:
                    pass
                else:
                    line_list.append(linename)
                    line_paths.append(segy_path)
    return line_list, line_paths

def make_vintage_dict(line_list, line_path):
    vintages = {}  # line_name (key): ['PSTM','UPSTM','MIG'] (value
    for root, dirs, files in os.walk(start_path):
        for f in files:
            #print (f)
            if file_is_segy(f):
                linename = remove_vintage(f)
                vintage = get_vintage(f)
                survey = root
                if linename in vintages.keys():
                    vintages[linename].append(vintage)
                else:
                    vintages[linename]= [vintage]
    return vintages

def make_survey_dict(line_list, line_path):
    survey = {}  # line_name (key): ['PSTM','UPSTM','MIG'] (value
    for root, dirs, files in os.walk(start_path):
        for f in files:
            #print (f)
            if file_is_segy(f):
                linename = remove_vintage(f)
                trimpath = root[2:]
                survey[linename] = trimpath[:trimpath.find('/')]
                print (survey[linename])
    return survey
        
# Creates a dictionary of seismic surveys and attributes 
# from NS_seismic_surveys.csv
# In order to create attribute fields for seismic shapefiles.
# Open file and read rows into rowdata
with open('NS_seismic_surveys.csv', 'r') as f:
     rowdata = [row for row in csv.reader(f.read().splitlines())]
     f.close()

keys = rowdata[0]

surveys = { }

for row in rowdata:
    survey_fields = row[1:]
    surveys[row[0]] = survey_fields

start_path = "."

line_list, segylist = path_to_lines(".")

print (segylist)

vintages = make_vintage_dict(line_list, segylist)
print ('Done making Vintage Dictionary')

# get survey dictionary read from file directory
survey_dict = make_survey_dict(line_list, segylist)
print ('Done linking line to survey')

driverName = "ESRI Shapefile"
drv = ogr.GetDriverByName( driverName )
if drv is None:
    print("{0} driver not available.\n".format(driverName))
    sys.exit( 1 )

output = "Antigonish_Cape_Breton_Seismic.shp"
overwrite = True
if overwrite:
  drv.DeleteDataSource(output)

ds = drv.CreateDataSource(output)
if ds is None:
    print("Creation of output file failed.\n")
    sys.exit( 1 )

lyr = ds.CreateLayer( "line_out", None, ogr.wkbLineString )
if lyr is None:
    print("Layer creation failed.\n")
    sys.exit( 1 )

# ***append any attributes here
fieldnames = ["SURVEY","REGION","YEAR", "COMPANY", "ONSHORE","TYPE" # 1-6 USER To get from survey dictionary
               ,"LINE_NAME"           # 7 i.e. NS03-WIND-04
               ,"VINTAGES"            # 8 list of vintages: i.e. [PSTM, UPSTM, STK] 
               ,"PREF_VIN"
               ,"NSAMPS" 
               ,"NTRACES" 
               ,"SAMPRATE" 
               ,"RECLENGTH"
                #"PREF_VINT ",         # 9 Leave Blank, enter manually afterwards [PSTM]
                #"CONTRACTOR",         # 10 available from textual header
                #"PROCESSOR",          # 11 available from textual header?    
                #"NSAMPLES",           # 12 number of samples per trace
                #"SAMPLE_RAsuveTEMS",      # 13 sample rate in milliseconds
                #"NTRACES",            # 14 number of traces / cdps / shotpoints
                #"BYTE_LOC_SP",        # 15 Byte locations of SP number: e.g. 17-20
                #"BYTE_LOC_CDP",       # 16 Byte locations of CDP number; e.g. 21-24
                #"BYTE_LOC_XCOORD",    # 17 Byte locations of X-coordinate; e.g. 73-76
                #"BYTE_LOC_YCOORD",    # 18 Byte locations of Y-coordinate 77-80
                #"COORD_SC",           # 19 Coordinate scaling factor
                #"TRACE_START_TIME"    # 20 Trace start time [ms]
                #"TRACE_END_TIME"      # 21 Trace end time [ms]
                ]
fieldtypes = [ogr.OFTString, ogr.OFTString, ogr.OFTInteger, ogr.OFTString, ogr.OFTString, ogr.OFTString 
                ,ogr.OFTString   # LINE_NAME 
                ,ogr.OFTString    # VINTAGES
                ,ogr.OFTString    # PREF_VIN
                ,ogr.OFTString    # NSAMPS
                ,ogr.OFTString    # NTRACES
                ,ogr.OFTString    # SAMPRATE
                ,ogr.OFTString    # RECLENGTH
                ]
fieldwidths = [32,32,4,32,10,3 
                ,50              # LINE_NAME 
                ,100             # VINTAGES
                ,12              # PREFFERED VINTAGE
                ,10              # NSAMPS
                ,10              # NTRACES
                ,10              # SAMPRATE
                ,10              # RECLENGTH
                ]
fieldprecisions = [0,0,0,0,0,0
                    ,0           # LINE_NAME
                    ,0           # VINTAGES
                    ,0           # PREFFERED VINTAGE
                    ,0           # NSAMPS
                    ,0           # NTRACES
                    ,0           # SAMPRATE
                    ,2           # RECLENGTH
                    ]

# Set Up
for fieldname, fieldtype, fieldwidth, fieldprecision in zip(fieldnames, fieldtypes, fieldwidths, fieldprecisions):
  field_defn = ogr.FieldDefn( fieldname, fieldtype )
  field_defn.SetWidth( fieldwidth )
  field_defn.SetPrecision( fieldprecision )

  if lyr.CreateField ( field_defn ) != 0:
    print("Creating {0} field failed.\n".format(fieldname))
    sys.exit( 1 )

for segyfile in segylist:
  
  # Read Source Data (Read Method 2 from obspy)
  #filename = datadir + "/" + segyfile
  stream = streamReadSEGY(segyfile, unpack_headers=True)

  print("START")
  #print("SEGY File: {0}".format(segyfile))
  print("SEGY File: {0}".format(segyfile),file=sys.stderr)
  print("")

  # Construct Feature
  feat = ogr.Feature( lyr.GetLayerDefn())
  #feat.SetField( "LINE", segyfile[segyfile.rfind('/'):] )
  #feat.SetField( "SAMPLEMS", stream.stats.binary_file_header.sample_interval_in_microseconds )
  seismicline = ogr.Geometry(ogr.wkbLineString)

  # Read Source Data (Read Method 3 from obspy)
  segy = readSEGY(segyfile, unpack_headers=True)
  dtstart = UTCDateTime("1970-01-01T00:00:00.000000Z")
  dtend = UTCDateTime("1970-01-01T00:00:01.600000Z")
  line = segyfile[segyfile.rfind('/')+1:segyfile.rfind('_')]
  print (line)

  # Load Trace Data Coordinates
  for i,trace in enumerate(segy.traces):
    coord_scalar = trace.header.scalar_to_be_applied_to_all_coordinates
    if coord_scalar == -100:
        coord_gain = 0.01
    elif coord_scalar == -10:
        coord_gain = 0.1
    else: 
        coord_gain = 1.0
    
    x = float(trace.header.source_coordinate_x)*coord_gain
    y = float(trace.header.source_coordinate_y)*coord_gain
    
    nsamples = trace.header.number_of_samples_in_this_trace
    sample_rate = trace.header.sample_interval_in_ms_for_this_trace 
    ntraces = len(stream.traces)

    # Sanity check for geometry having correct order of magnitude
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
    
    # Fix NS02-CUMB-55X lines
    #if line in ['NS02-CUMB-551','NS02-CUMB-552','NS02-CUMB-553','NS02-CUMB-554']:
    #    x *= 100.0
    #    y *= 100.0
            
    seismicline.AddPoint(x,y)
    
    survey = survey_dict[line]
    region = surveys[survey][0]
    year =  surveys[survey][1]
    company =  surveys[survey][2]
    onshore = surveys[survey][3]
    stype = surveys[survey][4]
    format = surveys[survey][5]
    status = surveys[survey][6]
    pref_vin = surveys[survey][7]
    
    feat = ogr.Feature( lyr.GetLayerDefn())
    feat.SetField( "SURVEY", survey )
    feat.SetField( "REGION", region )
    feat.SetField( "YEAR", year )
    feat.SetField( "COMPANY", company )
    feat.SetField( "ONSHORE", onshore )
    feat.SetField( "TYPE", stype )
    feat.SetField( "FORMAT", format )
    feat.SetField( "STATUS", status )
    feat.SetField( "LINE_NAME", line )
    feat.SetField( "VINTAGES", str(vintages[line]).strip('[]') )
    feat.SetField( "PREF_VIN", pref_vin )
    feat.SetField( "NSAMPS", nsamples  )
    feat.SetField( "NTRACES", ntraces  )
    feat.SetField( "SAMPRATE", sample_rate )
    feat.SetField( "RECLENGTH", nsamples*sample_rate/1000.0 )

  # Add Feature to Layer
  feat.SetGeometry(seismicline)
  if lyr.CreateFeature(feat) != 0:
    print("Failed to create feature in shapefile.\n")
    sys.exit( 1 )

  # Start Again
  feat.Destroy()

  print("END")
  print("")
  print("")
    
ds = None

