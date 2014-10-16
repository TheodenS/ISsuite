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

infile="/Users/security/science/bigoutput.csv"


dataf = DataFrame.from_csvfile(infile, sep = ",",header=True)

# Get statistics for investigated seqs
#rmean = robjects.r['mean']
#rmed = robjects.r['median']
#rmax = robjects.r['max']
#rsd = robjects.r['sd']
#rsum = robjects.r['sum']
#
#ma=rmax(dataf.rx('Length'))
#
#as_vec = robjects.r['as.vector']
#as_num = robjects.r['as.numeric']
#as_mat = robjects.r['as.matrix']
#
#test22=as_vec(dataf.rx('Length'))
#test23=as_num(test22[0])
#
#r_base = importr('base')
#newsum= r_base.summary(test23)
##for k, v in newsum.items():
##    #print("%s: %.3f\n" %(k, v))
##    print k
##    print v
##print "end summary"
#
#
##print "median"
##text_log+="median: "+str(rmed(test23)[0])+end
##text_log+="average: "+str(rmean(test23)[0])+end
##text_log+="sum: "+str(rsum(test23)[0])+end
#
#roughbin= round(ma[0]/100)
#bins=round(roughbin/100)*100


#ma2=rmax(ed)

#dataf_subset = dataf.rx(dataf.rx2("contig").ro >= 18, true)

scales = importr('scales')

gp = ggplot2.ggplot(dataf)
	#geom_histogram(aes(y = ..density..))
	#   ggplot2.geom_density()+\

	    # pp = gp + ggplot2.aes_string(x='%s(contrrr)') +  ggplot2.geom_histogram()+ggplot2.scale_y_sqrt()
bins=10
teest3=robjects.r('theme(axis.text.x=element_text(angle=90))')

pp = gp + \
ggplot2.aes_string(x='Length') +  \
ggplot2.geom_histogram()+\
ggplot2.ggtitle("Found IS fragment lengths")+ \
ggplot2.scale_x_continuous(name="fragment lengths, bin="+str(bins),breaks=scales.pretty_breaks(20)) +\
ggplot2.scale_y_continuous(labels=scales.comma,name="Count",breaks=scales.pretty_breaks(10))+ \
teest3
pp.plot()
robjects.r.ggsave("/Users/security/science/dna_subj_hist.pdf")
