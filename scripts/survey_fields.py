# Creates a dictionary of seismic surveys and attributes 
# from NS_seismic_surveys.csv
# In order to create attribute fields for seismic shapefiles.

import csv

# Open file and read rows into rowdata
with open('NS_seismic_surveys.csv') as f:
     rowdata = []
     reader = csv.reader(f)
     reader.next()
     for row in reader:
        rowdata.append(row)
     f.close()

# Make the first row the keys of the dictionary
keys = rowdata[0]

# Make an empty dictionaru
surveys = { }

# For each row, make a dictionary entry for each survey.
for row in rowdata:
    survey_fields = row[1:]
    surveys[row[0]] = survey_fields
