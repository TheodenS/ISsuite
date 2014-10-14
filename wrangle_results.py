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
import csv
#
parser = argparse.ArgumentParser()


parser.add_argument("-in_csv", help="")
#parser.add_argument("-in_csv_all", help="")
parser.add_argument("-out", help="")

args=parser.parse_args()


db="/Users/security/science/dest3.db"
con = sqlite3.connect(db)
cur = con.cursor()
#
#columns="orf_id,group_,GS670_3p0"
groupthing="Cyanobacteria"
#
#
csvh=open("/Users/security/Desktop/headers2.csv","r")
csread=csvh.read()
css=csread.split("\n")

heads=[]
for u in css:
    if u=="":
        continue
    uu=u.split(",")
    print u
    heads.append(uu[1])


dethg=""
for d in heads:
    if "3p0" in d:
        dethg+="SUM("+d+")+"
        #dethg+=d+","
print dethg


dethg=dethg[:-1]
#columns="SUM("+dethg+")"
columns="SUM("+"GS670_0p8"+")"
#columns=dethg
print columns
selstring="SELECT "+columns+" FROM allist WHERE group_='"+groupthing+"'"


with con:
        cur.execute(selstring)

reply=cur.fetchall()


broup="SELECT DISTINCT group_ FROM allist"
with con:
        cur.execute(broup)
orgroups=cur.fetchall()
#
#
dets=heads[33:87]
#print dets
#
had="mystation,"
for i in orgroups:
    print "text"+i[0]+"tes"
    hed=i[0]
    if len(hed)<5:
        hed="unknown"
    had+=hed.replace(" ","").replace("/","")+","

#raw_input()
#
had=had[:-1]
had+="\n"
outstring=had

#for det in dets:
#    outstring+=det+","
#    for i in orgroups:
#        print i[0]
#        columns="SUM("+det+")"
##columns=dethg
#        selstring="SELECT "+columns+" FROM allist WHERE group_='"+i[0]+"'"
#        with con:
#            cur.execute(selstring)
#        reply=cur.fetchall()
#        for u in reply:
#            print u[0]
#            outstring+=str(u[0])+","
#    outstring=outstring[:-1]
#    outstring+="\n"
#
#for det in dets:
#    outstring+=det+","
#    for i in orgroups:
#        print i[0]
#        columns="SUM("+det+")"
##columns=dethg
#        selstring="SELECT "+columns+" FROM allist WHERE isname!='' AND group_='"+i[0]+"'"
#        print selstring
#        #selstring="SELECT "+columns+" FROM allist WHERE group_='"+i[0]+"'"
#        with con:
#            cur.execute(selstring)
#        reply2=cur.fetchall()
#        for u in reply2:
#            print u[0]
#            outstring+=str(u[0])+","
#    outstring+="\n"

#output="/Users/security/science/output3.csv"
#csvWriter = open(output, "w")
#csvWriter.write(outstring)

output="/Users/security/science/fractions.csv"
#output="/Users/security/science/output4.csv"
#csvWriter = open(output, "w")
#csvWriter.write(outstring)

#
#
robjects.r('library(reshape2)')
robjects.r('library(RColorBrewer)')
robjects.r('library(stats)')
#
##dataf=robjects.r('mydataset <- read.delim("'+args.in_csv+'", sep="," , header=TRUE)')
dataf=robjects.r('mydataset <- read.delim("'+output+'", sep="," , header=TRUE)')
robjects.r('attach(mydataset)')



for depth in ["0p1","0p8","3p0"]:



    o3p0=robjects.r('only_depth<-mydataset[grep("'+depth+'", mydataset$mystation), ]')

# remove landsort
    robjects.r('new_d <- only_depth[! grepl("GS678",only_depth$mystation),]')

# only landsort
    #robjects.r('onlyls <- only_depth[grepl("GS678",only_depth$mystation),]')
    
    robjects.r('testtt<-melt(new_d)')

    robjects.r('cols <- colorRampPalette(brewer.pal(9, "Accent"))')
    robjects.r('print(cols)')

    allISfiltersgraph=robjects.r('allplot<-ggplot(data=testtt, aes(x=testtt$mystation, y=testtt$value,fill=testtt$variable))') 
    allISfiltersgraph=robjects.r('allplot<-allplot+geom_bar(color="black",stat="identity",position="fill",size=0.2)')

    allISfiltersgraph=robjects.r('allplot<-allplot+theme(axis.text.x=element_text(angle=90))')
    #allISfiltersgraph=robjects.r('allplot<-allplot+theme(axis.text.x=element_text(angle=90),axis.text.x = element_blank())')

    allISfiltersgraph=robjects.r('allplot<-allplot+scale_fill_manual(values=cols(19))')
    gols=robjects.r('scale_fill_manual(values=cols(19))')
    allISfiltersgraph=robjects.r('allplot<-allplot+guides(fill = guide_legend(reverse=TRUE,title = "Group", title.position = "top"))')
    allISfiltersgraph=robjects.r('allplot<-allplot+scale_y_continuous(breaks=seq(0, 1, 0.1),name="Fraction of transcripts")')
    allISfiltersgraph=robjects.r('allplot<-allplot+scale_x_discrete(name="Station")')

    allISfiltersgraph=robjects.r('allplot<-allplot+labs(title="'+"fractions_"+depth+'")')


    robjects.r.ggsave("/Users/security/science/"+depth+"_fractions.pdf")
    #robjects.r.ggsave(args.out.replace(".pdf","")+depth+".pdf")



