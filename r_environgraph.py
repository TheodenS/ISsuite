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
parser.add_argument("-in_csv_joined", help="")
parser.add_argument("-out", help="")
parser.add_argument("-plotfraction", help="")

args=parser.parse_args()



# Get statistics for investigated seqs, examples
####################
rmean = robjects.r['mean']
rmed = robjects.r['median']
rmax = robjects.r['max']
rsd = robjects.r['sd']
rsum = robjects.r['sum']
rseq = robjects.r['seq']


r_base = importr('base')



dataf2=robjects.r('joined <- read.delim("'+args.in_csv_joined+'", sep="," , header=TRUE)')

dataf=robjects.r('allfilters <- read.delim("'+args.in_csv+'", sep="," , header=TRUE)')

# Get rid of landsort
####################
smalls=robjects.r('joined_no_ls<-subset(joined,filtersize %in% c("0p1","0p8","3p0"))')
smalls1=robjects.r('allfilters_no_ls<-subset(allfilters,filtersize %in% c("0p1","0p8","3p0"))')


# make environmental graph pt1
ggobj=robjects.r('p<-ggplot(data=joined_no_ls)')
no3=robjects.r('p<-p+geom_line(aes(x=mystation, y=NO3,group=2,color="NO3"))' )
nh4=robjects.r('p<-p+geom_line(aes(x=mystation, y=NH4*10,group=2,color="NH4"))' )
urea=robjects.r('p<-p+geom_line(aes(x=mystation, y=Urea*10,group=2,color="Urea"))' )
po4=robjects.r('p<-p+geom_line(aes(x=mystation, y=PO4*10,group=2,color="PO4"))' )
np=robjects.r('p<-p+geom_line(aes(x=mystation, y=NP,group=2,color="NP"))' )
si=robjects.r('p<-p+geom_line(aes(x=mystation, y=Si,group=2,color="Si"))' )
dop=robjects.r('p<-p+geom_line(aes(x=mystation, y=DOP*10,group=2,color="DOP"))' )

ang=robjects.r('p<-p+theme(axis.text.x=element_text(angle=90))')
teest3=robjects.r('p<-p+coord_fixed(ratio=0.1)')
teest3.plot()
robjects.r.ggsave(args.out)


# make environmental graph pt2
#onlyfilters2=robjects.r('alltra22<-subset(datass,filtersize %in% c("0p1","0p8","3p0"))')
img2=robjects.r('plt2<-ggplot(data=joined_no_ls)')

robjects.r('plt2<-plt2+geom_line(aes(x=mystation, y=Sal,group=1,colour="Sal"))' )
robjects.r('plt2<-plt2+geom_line(aes(x=mystation, y=Temp,group=2,color="Temp"))' )
robjects.r('plt2<-plt2+geom_line(aes(x=mystation, y=pH,group=2,color="pH"))' )
robjects.r('plt2<-plt2+geom_line(aes(x=mystation, y=O,group=2,color="O"))' )
robjects.r('plt2<-plt2+geom_line(aes(x=mystation, y=Chla,group=2,color="Chla"))' )

robjects.r('plt2<-plt2+theme(axis.text.x=element_text(angle=90))')
img2=robjects.r('plt2<-plt2+coord_fixed(ratio=0.1)')

img2.plot()
robjects.r.ggsave((str(args.out).replace(".pdf",""))+"2.pdf")


# Make graph of all IS

if args.plotfraction=="True":
    allISfiltersgraph=robjects.r('allplot<-ggplot(data=allfilters_no_ls, aes(x=mystation, y=ratio,fill="mauve")) + geom_bar(stat="identity",position=position_dodge())+theme(axis.text.x=element_text(angle=90))')

# Plot hits instead of fractions of hits
if args.plotfraction=="False":
    allISfiltersgraph=robjects.r('allplot<-ggplot(data=allfilters_no_ls, aes(x=mystation, y=istranscripts,fill="mauve")) + geom_bar(stat="identity",position=position_dodge())+theme(axis.text.x=element_text(angle=90))')
#teest=robjects.r('p<-p+coord_fixed(ratio=400)')

robjects.r.ggsave((str(args.out).replace(".pdf",""))+"_hits_station.pdf")





#latits=dataf2.rx('Lat')
#maxval=rmax(latits)

# extra functions
#as_vec = robjects.r['as.vector']
#as_num = robjects.r['as.numeric']
#as_mat = robjects.r['as.matrix']

#lat_as_vector=as_vec(latits)

#lat_as_num=as_num(lat_as_vector[0])

#print "summary for lattitude"
#newsum= r_base.summary(lat_as_num)
#for k, v in newsum.items():
#    print k
#    print v
#print "end summary"



#roughbin= round(maxval[0]/100)
#bins=round(roughbin/100)*100


#ma2=rmax(ed)

scales = importr('scales')


graphtitle="Latitudes"
xlabel="xlabel"
ylabel="ylabel"

#names = robjects.StrVector(("blue", "red"))
grdevices = importr('grDevices')
clustsim = importr('clusterSim')


onlyfilts=robjects.r('onlyfilt<-subset(joined,filtersize %in% c("0p1","0p8","3p0"))')
#onlysurf=robjects.r('onlyfilt2<-subset(datass,Sampledepth %in% c("0p1","0.30487806"))')
#onlysurfno=robjects.r('onlyfilt2<-')
onlysurf=robjects.r('onlysurface<-onlyfilt[ which("Sampledepth" =="0.30487806"),]')
#print "onlysurf"
#print onlysurf

#colours2 = grdevices.topo_colors(10)
colours2 = grdevices.cm_colors(10)
#colours2 = grdevices.rainbow(20)
#print colours2
#colours = ggplot2.rainbow(54)
#bins=10
gp = ggplot2.ggplot(onlysurf)
#gp = ggplot2.ggplot(onlyfilts)

gp=gp+ggplot2.aes_string(x="Lon", y="Lat", col="Temp",label="Station")
gp=gp+ggplot2.scale_colour_gradientn(colours=colours2)
gp=gp+ggplot2.geom_text(col="black",offset = 10)
gp=gp+ggplot2.geom_point(position="jitter")
gp=gp+ggplot2.ggtitle(graphtitle)

robjects.r('library(ggmap)')
robjects.r('library(mapproj)')
robjects.r('map <- get_map(location = "Europe", zoom = 4)')
robjects.r('ggmap(map)')

#robjects.r('library(maps)')

#robjects.r('map("world", interior = FALSE)')

#robjects.r('map("state", boundary = FALSE, col="gray", add = TRUE)')
#gp.plot()

'''
pp = gp + \
ggplot2.aes_string(x="Lon", y="Lat", col="Temp",label="Station") + \
ggplot2.scale_colour_gradientn(colours=colours2)+ \
ggplot2.geom_text(col="black")+ \
ggplot2.geom_point()

ggplot2.ggtitle(graphtitle)
pp.plot()
'''
#robjects.r.ggsave((str(args.out).replace(".pdf",""))+"map.pdf")

onlyfiltxxx=robjects.r('onl<-subset(allfilters,filtersize %in% c("0p1","0p8","3p0"))')
#print "here is"
#print robjects.r('print(datass1)')
#print "there was"
'''
ggplot2.scale_x_continuous(name=xlabel,breaks=scales.pretty_breaks(20)) +\
ggplot2.scale_y_continuous(labels=scales.comma,name=ylabel,breaks=scales.pretty_breaks(10))+ ggplot2.theme(title=ggplot2.element_text(colour="blue",face="bold"))''' 

if args.plotfraction=="False":
    barsdata=robjects.r('p<-ggplot(data=onl, aes(x=mystation, y=istranscripts,fill=filtersize)) + geom_bar(stat="identity",position=position_dodge())+theme(axis.text.x=element_text(angle=90))')
if args.plotfraction=="True":
    barsdata=robjects.r('p<-ggplot(data=onl, aes(x=mystation, y=ratio,fill=filtersize)) + geom_bar(stat="identity",position=position_dodge())+theme(axis.text.x=element_text(angle=90))')
#barsdata=robjects.r('p<-ggplot(data=datass1, aes(x=mystation, y=ratio,fill=filtersize)) + geom_bar(stat="identity",position=position_dodge())+theme(axis.text.x=element_text(angle=90))')
#teest=robjects.r('p<-p+coord_fixed(ratio=400)')
#barsg = ggplot2.ggplot(onlysurf)
#barsg=barsg+ggplot2.aes_string(x="mystation", y="ratio")
#barsg=barsg+ggplot2.geom_bar(stat="identity")
#bbb=robjects.r('p<-ggplot(data=smalls11, aes(x=mystation, y=ratio,fill=filtersize)) + geom_bar(stat="identity",position=position_dodge())+theme(axis.text.x=element_text(angle=90))')
#barsg=barsg+ggplot2.geom_bar(stat="identity",position=ggplot2.position_dodge())
#barsg=barsg+ggplot2.theme(axis.text.x=ggplot2.element_text(angle=90))
barsdata.plot()
robjects.r.ggsave((str(args.out).replace(".pdf",""))+"allfilters.pdf")


#dataframe=robjects.r('mydata <- read.delim("/Users/security/science/transcless3ke5/joined_all.csv", sep="," , header=TRUE)')
