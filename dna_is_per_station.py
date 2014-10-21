import os
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import xlrd
import sqlite3
import rpy2.robjects as robjects
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.packages import importr
import pandas as pd
import numpy as np
import rpy2.robjects.pandas2ri
from rpy2.robjects.vectors import DataFrame
import math
import datetime

parser = argparse.ArgumentParser()


parser.add_argument("-in_csv", help="")
parser.add_argument("-out", help="")

args=parser.parse_args()
folders=getdir(indir)

stations={}

for folder in folders:
    if "done.tar.gz" in folder:

        #county+1

        #if county>2:
        #    continue
        print folder
        ungzcom="tar -xzf " +indir+folder+" -C /Users/security/science/"
        print ungzcom
        os.system(ungzcom)

        untarred_folder="/Users/security/science/Users/"

        untarred_folder_fh=open("/Users/security/science/Users/security/science/metagen2009_2/"+folder.replace("_done.tar.gz","")+"/allgenomesinfo.csv","r")
        untarred_read=untarred_folder_fh.read()
        splitt=untarred_read.split("\n")
        splitt=splitt[1:]
        loc=folder.replace("_done.tar.gz","")
        print loc
        raw_input()
        for u in splitt:
            if u=="":
                continue
            uu=u.split(",")
            


            #if not uu[1] in namedic.keys():
        
         
            #outline+=loc+","+uu[0]+","+uu[2]+"\n"
            
        
        untarred_folder="/Users/security/science/Users/"
        rmline="rm -R "+untarred_folder
        os.system(rmline)

