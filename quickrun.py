import os
import os.path
import argparse
import datetime
import sqlite3

# Appends text to file
# Parse args
parser = argparse.ArgumentParser()
parser.add_argument("-basedir", help="dirwithall")
parser.add_argument("-basefastafile", help="dirwithall")
parser.add_argument("-queryfastafile", help="dirwithall")
parser.add_argument("-maxeval", help="dirwithall")
args=parser.parse_args()

basedir=args.basedir




# Make dirs
if not os.path.isdir(args.basedir):
     os.makedirs(args.basedir)

environmentaldata="/Users/security/science/RNA/mimebs_metadata/metadata.2009.csv"


if True:
    prog="/Users/security/science/software/ISsuite/"+"wrangle_results_sigA.py"
    runstring="python "+prog

    runstring+=" -in_db="+basedir+"appended_large.db"
    runstring+=" -out="+basedir+"bajs"
    runstring+=" -outcsv="+basedir+"summary.csv"
    runstring+=" -environmentaldata="+environmentaldata
    runstring+=" -basepicout="+basedir

    os.system(runstring)
# wrangle_results.py uses this
#    csvh=open("/Users/security/Desktop/headers2.csv","r")

print "finished program"




