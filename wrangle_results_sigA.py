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
import random


def get_map_settings(title,yname):
#graphtitle="TESTTI"

    plottitlesize=str(8)
    axistextsize=str(8)
    axistitlesize=str(8)
    xname="Longditude"

    #robjects.r('m<-m +scale_y_continuous(breaks = (-2:2) * 30)')
    #robjects.r('m<-m +scale_x_continuous(breaks = (-4:4) * 45)')
    map_settings='ggtitle("'+title+'")' 
    map_settings+='+theme('
    map_settings+='plot.title = element_text(size='+plottitlesize+',lineheight=.8)'
    map_settings+=',axis.text.x=element_text(size='+axistextsize+',angle=0)'
    map_settings+=',axis.text.y=element_text(size='+axistextsize+',angle=0)'
    map_settings+=',axis.title.x=element_text(size='+axistitlesize+')'
    map_settings+=',axis.title.y=element_text(size='+axistitlesize+')'
    map_settings+=',legend.position="none"'
    map_settings+=')'
    map_settings+='+scale_x_continuous(name="'+xname+'")'     
    map_settings+='+scale_y_continuous(name="'+yname+'")'     

#teest=robjects.r('p<-p+coord_fixed(ratio=400)')
    return map_settings

def get_bar_settings(title,yname):
#graphtitle="TESTTI"

    mycolorfirst = "%06x" % random.randint(0,0xFFFFFF)
    mycolor="#"+str(mycolorfirst).upper()
    print mycolor
    plottitlesize=str(8)
    axistextsize=str(8)
    axistitlesize=str(8)
    xname="station"

    robjects.r('library(scales)')
    bar_settings='geom_bar(stat="identity",position=position_dodge(),fill="'+mycolor+'")'    

    bar_settings+='+ggtitle("'+title+'")' 
    bar_settings+='+theme('
    bar_settings+='plot.title = element_text(size='+plottitlesize+',lineheight=0.8)'
    bar_settings+=',axis.text.x=element_text(size='+axistextsize+',angle=0)'
    bar_settings+=',axis.text.y=element_text(size='+axistextsize+',angle=0)'
    bar_settings+=',axis.title.x=element_blank()'
    bar_settings+=',axis.title.y=element_text(size='+axistitlesize+')'
    bar_settings+=',legend.position="none"'
    bar_settings+=')'
    bar_settings+='+scale_x_discrete(name="'+xname+'")'     
    bar_settings+='+scale_y_continuous(labels=comma, name="'+yname+'")'     

#teest=robjects.r('p<-p+coord_fixed(ratio=400)')
    return bar_settings


# Trying with maps
##############

#robjects.r('library(maps)')
#robjects.r('map("world", interior = FALSE)')
#robjects.r('map("state", boundary = FALSE, col="gray", add = TRUE)')
#gp.plot()
#mp <- NULL
#mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
#mp <- ggplot() + mapWorld
# 
##Now Layer the cities on top
#mp <- mp+ geom_point(aes(x=visit.x, y=visit.y) ,color="blue", size=3)
#mp

def bar_allsize(incsv,outname,w,h,graphtitle,plotvar):
    #graphtitle="Total number of reads"
#robjects.r('mylocation <- c(lon=18,lat=61)')
    robjects.r('library(ggmap)')
    station_dataset=robjects.r('chosendataset <- read.delim("'+incsv+'", sep="," , header=TRUE)')
    allISfiltersgraph=robjects.r('oio2<-ggplot(data=chosendataset, aes(x=mystation, y='+"n_hit_station/alltranscripts"+',fill="green")) + geom_bar(stat="identity",position=position_dodge())+theme(axis.text.x=element_text(angle=90))')
    #allISfiltersgraph=robjects.r('oio2<-ggplot(data=chosendataset, aes(x=mystation, y='+plotvar+',fill="mauve")) + geom_bar(stat="identity",position=position_dodge())+theme(axis.text.x=element_text(angle=90))')
#teest=robjects.r('p<-p+coord_fixed(ratio=400)')
    robjects.r('oio2 <- oio2+ggtitle("'+graphtitle+'")')
    robjects.r('oio2 <- oio2+theme(plot.title = element_text(lineheight=.8, face="bold"))')
    #robjects.r.ggsave("/Users/security/science/ttest.pdf")
    robjects.r('pdf("'+outname+'", width='+str(w)+',height='+str(h)+')')
    im4=robjects.r('print(oio2)')
    robjects.r('dev.off()')

    #station_dataset=robjects.r('chosendataset <- read.delim("'+incsv+'", sep="," , header=TRUE)')
    ##station_dataset=robjects.r('chosendataset <- read.delim("'+'/Users/security/science/output.csv'+'", sep="," , header=TRUE)')

    #robjects.r('pdf("'+outname+'")')
    #robjects.r('oio <- ggmap(map) + geom_point(data=chosendataset,color="red",aes(x=Lon, y = Lat,size=alltranscripts))')
    #robjects.r('oio <- oio + geom_text(data=chosendataset,size=3,color="red",aes(x=Lon,y=Lat,label=mystation),hjust=0.5, vjust=-0.7)')

def gogmap_google(incsv,outname,w,h):

#robjects.r('mylocation <- c(lon=18,lat=61)')
    robjects.r('library(ggmap)')
    station_dataset=robjects.r('chosendataset <- read.delim("'+incsv+'", sep="," , header=TRUE)')
    #robjects.r('map <- get_googlemap(maptype ="satellite",center = c(lon=18,lat=61.5),zoom=4,color="bw")')
    robjects.r('map <- get_googlemap(maptype ="hybrid",center = c(lon=18,lat=61.5),zoom=4,color="bw")')
    #robjects.r('map <- get_googlemap(maptype ="roadmap",center = c(lon=18,lat=61.5),zoom=4,color="bw")')
    #robjects.r('map <- get_googlemap(maptype ="terrain",center = c(lon=18,lat=61.5),zoom=4,color="bw")')

    station_dataset=robjects.r('chosendataset <- read.delim("'+incsv+'", sep="," , header=TRUE)')
    #station_dataset=robjects.r('chosendataset <- read.delim("'+'/Users/security/science/output.csv'+'", sep="," , header=TRUE)')

    robjects.r('pdf("'+outname+'")')
    robjects.r('oio <- ggmap(map) + geom_point(data=chosendataset,color="red",aes(x=Lon, y = Lat,size=alltranscripts))')
    robjects.r('oio <- oio + geom_text(data=chosendataset,size=3,color="red",aes(x=Lon,y=Lat,label=mystation),hjust=0.5, vjust=-0.7)')
    im4=robjects.r('print(oio)')
    robjects.r('dev.off()')

def with_rworld(incsv,outname,w,h):

    station_dataset=robjects.r('chosendataset <- read.delim("'+incsv+'", sep="," , header=TRUE)')

    robjects.r('library(grid)')
    robjects.r('library(gridExtra)')
    
    robjects.r('library(rworldmap)')
    robjects.r('worldMap <- getMap(resolution="high")')

    robjects.r('world.points <- fortify(worldMap)')
    robjects.r('world.points$region <- world.points$id')
    robjects.r('world.df <- world.points[,c("long","lat","group", "region")]')
    robjects.r('geom_polygon(data = world.df, aes(x = long, y = lat, group = group))')
    robjects.r('m <- ggplot() ')

    robjects.r('m<-m +geom_path(data = world.df, aes(x = long, y = lat, group = group))')
    robjects.r('m<-m + theme(legend.position="none")')
    #robjects.r('m<-m +scale_y_continuous(breaks = (-2:2) * 30)')
    #robjects.r('m<-m +scale_x_continuous(breaks = (-4:4) * 45)')
    robjects.r('m<-m +coord_map("ortho", orientation=c(10,16, 0),xlim=c(8,24),ylim=c(54,70))')
    robjects.r('m<-m +geom_point(data=chosendataset,color="red",aes(x=Lon, y = Lat,size=chosendataset$alltranscripts))')
    robjects.r('m <- m+ geom_text(data=chosendataset,size=3,color="red",aes(x=Lon,y=Lat,label=mystation),hjust=0.5, vjust=-0.7)')
   
    #robjects.r('pdf("'+outname+'")')
    robjects.r('pdf("'+outname+'", width='+str(w)+',height='+str(h)+')')
    im4=robjects.r('print(m)')
    robjects.r('dev.off()')







def gogmap3(incsv_surf,incsv_mid,outname,w,h):
    robjects.r('library(gridExtra)')
    robjects.r('library(grid)')
    robjects.r('library(ggmap)')
    robjects.r('library(rworldmap)')
    robjects.r('library(rworldxtra)')
# Three bar plots for shallow
    robjects.r('chosendataset_surf <- read.delim("'+incsv_surf+'", sep="," , header=TRUE)')



    mycolorfirst = "%06x" % random.randint(0,0xFFFFFF)
    mycolor="#"+str(mycolorfirst).upper()
    mycolor="blue"
    print mycolor

    graphtitle="Total_number_of_transcripts"
    plotline='bargall_shallow<-ggplot(data=chosendataset_surf, aes(x=mystation, y='+"alltranscripts"+'),fill="blue",colour="blue")'
    robjects.r(plotline)
    #robjects.r('bargall_shallow<-ggplot(data=chosendataset_surf, aes(x=mystation, y='+"alltranscripts"+',fill="'+mycolor+'"))')
    lol='bargall_shallow<-bargall_shallow+'+get_bar_settings(graphtitle,"Reads_with_ISs")
    print lol
    robjects.r(lol)

    graphtitle="IS fraction of transcripts"
    station_dataset=robjects.r('chosendataset_surf_shallow <- read.delim("'+incsv_surf+'", sep="," , header=TRUE)')
    allISfiltersgraph=robjects.r('bargall_shallow_frac<-ggplot(data=chosendataset_surf_shallow, aes(x=mystation, y='+"n_hit_station/alltranscripts"+',fill="green")) ')
    lol2='bargall_shallow_frac<-bargall_shallow_frac+'+get_bar_settings(graphtitle,"IS frac. of total")
    robjects.r(lol2)

    graphtitle="Number of IS transcripts"
    station_dataset=robjects.r('chosendataset_surf_shallow <- read.delim("'+incsv_surf+'", sep="," , header=TRUE)')
    allISfiltersgraph=robjects.r('bargall_shallow_is<-ggplot(data=chosendataset_surf_shallow, aes(x=mystation, y='+"n_hit_station"+',fill="deeppink4")) ')
    lol3='bargall_shallow_is<-bargall_shallow_is+'+get_bar_settings(graphtitle,"Total reads")
    robjects.r(lol3)

# Three bar plots for medium 
    station_dataset=robjects.r('chosendataset_mid <- read.delim("'+incsv_mid+'", sep="," , header=TRUE)')


    graphtitle="Total_number_of_transcripts"
    robjects.r('bargall_mid_total_transcripts<-ggplot(data=chosendataset_mid, aes(x=mystation, y='+"alltranscripts"+',fill="'+mycolor+'"))')
    lol='bargall_mid_total_transcripts<-bargall_mid_total_transcripts+'+get_bar_settings(graphtitle,"Total reads")
    robjects.r(lol)

    graphtitle="Total number of hits"
    robjects.r('bargall_mid_total<-ggplot(data=chosendataset_mid, aes(x=mystation, y='+"n_hit_station"+',fill="green"))')
    lol='bargall_mid_total<-bargall_mid_total+'+get_bar_settings(graphtitle,"Reads with ISs")
    robjects.r(lol)

    graphtitle="fraction_of_total"
    robjects.r('bargall_mid_frac<-ggplot(data=chosendataset_mid, aes(x=mystation, y='+"n_hit_station/alltranscripts"+',fill="'+mycolor+'"))')
    lol='bargall_mid_frac<-bargall_mid_frac+'+get_bar_settings(graphtitle,"IS frac. of total")
    robjects.r(lol)


# Make map
#########

    robjects.r('worldMap <- getMap(resolution="high")')
    robjects.r('world.points <- fortify(worldMap)')
    robjects.r('world.points$region <- world.points$id')
    robjects.r('world.df <- world.points[,c("long","lat","group", "region")]')
    robjects.r('geom_polygon(data = world.df, aes(x = long, y = lat, group = group))')


    robjects.r('m <- ggplot() ')
    robjects.r('m<-m +geom_path(data = world.df, size=0.3,aes(x = long, y = lat, group = group))')
    gostr="m<-m + "+ get_map_settings("","Latitude")
    robjects.r(gostr)

# scandic view
    robjects.r('mymap<-m +coord_map("ortho", orientation=c(10,16, 0),xlim=c(8,24),ylim=c(54,70))')
# space view
    #robjects.r('m<-m +coord_map("ortho", orientation=c(10,16, 0),xlim=c(-50,49),ylim=c(10,70))')


# point size by alltranscripts
    #graphtitle="title"



    robjects.r('surf_points<-mymap +geom_point(data=chosendataset_surf,color="red",aes(x=Lon, y = Lat,size=alltranscripts))')

    robjects.r('surf_points <- surf_points+ geom_text(data=chosendataset_surf,size=3,color="red",aes(x=Lon,y=Lat,label=mystation),hjust=0.5, vjust=-0.7)')
    #robjects.r('surf_points<-surf_points +ggtitle("'+graphtitle+'")')



    robjects.r('mid_points<-mymap +geom_point(data=chosendataset_mid,color="#000000",aes(x=Lon, y = Lat,size=alltranscripts))')

    robjects.r('mid_points <- mid_points+ geom_text(data=chosendataset_mid,size=3,color="red",aes(x=Lon,y=Lat,label=mystation),hjust=0.5, vjust=-0.7)')
    robjects.r('mid_points <- mid_points+ggtitle("'+graphtitle+'")')


    robjects.r('pdf("'+outname+'", width='+str(w)+',height='+str(h)+')')

    grid_arrangement='grid.arrange('
    #grid_arrangement='arrangeGrob('
    grid_arrangement+='surf_points,'
    grid_arrangement+='arrangeGrob('
    grid_arrangement+='bargall_shallow_frac,bargall_shallow_is,bargall_shallow'
    grid_arrangement+=',ncol=1'
    grid_arrangement+=',heights=c(2,1,1)),'
    #grid_arrangement+=',ncol=2)'
    grid_arrangement+='mid_points,'
    grid_arrangement+='arrangeGrob(bargall_mid_frac,bargall_mid_total_hits,bargall_mid_total_hits,ncol=1,heights=c(2,1,12,1,1)),'
    #grid_arrangement+='deep_points,'
    #grid_arrangement+='arrangeGrob(bargall_deep,bargall_deep,ncol=1),' 
    grid_arrangement+='ncol=2,'
    grid_arrangement+='main=textGrob("Insertion sequences in the Baltic Sea 2009", gp = gpar(fontsize=18, fontface="bold.italic", fontsize=18))'
    grid_arrangement+=',widths=1:1)'

    hh='grid.arrange( arrangeGrob( surf_points,arrangeGrob(bargall_shallow_frac,bargall_shallow_is,bargall_shallow,ncol=1,heights=c(2,1,1)) ,ncol=2,main="Shallow"),arrangeGrob(mid_points,arrangeGrob(bargall_mid_frac,bargall_mid_total_hits,bargall_mid_total_transcripts,ncol=1,heights=c(2,1,1)) ,ncol=2,main="Deep") ,ncol=1,main=textGrob("\nInsertion sequences in the Baltic Sea 2009\n", gp = gpar(fontsize=18, fontface="bold", fontsize=18)),widths=1:1)'




    print grid_arrangement
    robjects.r(hh)
    #robjects.r(grid_arrangement)
    #robjects.r('grid.arrange( surf_points,mid_points,deep_points, ncol=2)')

    #im4=robjects.r('print(m)')
    robjects.r('dev.off()')


def getfromcc(name):
    print name
    badst= "nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2,nod2"
    acsv=""
    for lin in metalines:
    #for cline in clines:
        if lin=="":
            continue
        #print lin
        #print cline
        lintabs=lin.split("\t")
        #clinetab=cline.split(",")
        me1name=lintabs[0]
        if me1name==name:
            for m in lintabs:
                acsv+=str(m)+","
            #print acsv
            new=acsv[:-1]
            #print new
            return new 

    return badst 

parser = argparse.ArgumentParser()
#parser.add_argument("-in_csv", help="")
parser.add_argument("-in_db", help="")
#parser.add_argument("-in_csv_all", help="")
parser.add_argument("-out", help="")
parser.add_argument("-outcsv", help="")
parser.add_argument("-environmentaldata", help="")
parser.add_argument("-basepicout", help="")

args=parser.parse_args()
basepicdir=args.basepicout
metadata="/Users/security/science/RNA/mimebs_metadata/metadata.2009.csv"

# Read environmental csv
###########
metafh=open(args.environmentaldata,"r")
metaread=metafh.read()
metalines=metaread.split("\n")
metahead=metalines[0]
metalines=metalines[1:]

#db="/Users/security/science/dest3.db"
db=args.in_db
con = sqlite3.connect(db)
cur = con.cursor()

# what columns are there?
# made with pramga t outcsv
##############
csvh=open("/Users/security/Desktop/headers2.csv","r")
csread=csvh.read()
css=csread.split("\n")

heads=[]
for u in css:
    if u=="":
        continue
    uu=u.split(",")
    heads.append(uu[1])

print "heads"
print heads


stations=[]
for i in range(33,87):
   # print heads[i]
    stations.append(heads[i])




dethg=""
for d in heads:
    if "3p0" in d:
        dethg+="SUM("+d+")+"
        #dethg+=d+","


dethg=dethg[:-1]
#columns="SUM("+dethg+")"
#columns=dethg
print "dethg"
print dethg


summarycsv=""


#selstring="SELECT "+"GS670_0p8"+" FROM allist WHERE group_='"+groupthing+"'"
#selstring="SELECT "+columns+" FROM allist WHERE group_='"+groupthing+"'"

# Sample one column
############
columns="SUM("+"GS670_0p8"+")"
selstring="SELECT "+columns+" FROM allist"
with con:
        cur.execute(selstring)
reply=cur.fetchall()
onesampl=reply[0][0]
summarycsv+="Sample: sum of all transcripts at GS670_0p8:,"+str(onesampl)+"\n"


# Sample one column with hits
############
columns="SUM("+"GS670_0p8"+")"
selstring="SELECT "+columns+" FROM allist WHERE hit_abb!='nod'"
#selstring="SELECT "+columns+" FROM allist WHERE isname!=''"
with con:
        cur.execute(selstring)
reply=cur.fetchall()
onesampl=reply[0][0]
summarycsv+="Sample: sum of all hit transcripts at GS670_0p8:,"+str(onesampl)+"\n"


# compile a list with stations
###########
stationcsv="name,mystation,filtersize,alltranscripts,n_hit_station,"+metahead.replace("\t",",")+"\n"
#for mystation in stations:
#    print "getting hit counts for "+mystation
#    stat_loc=mystation.split("_")[0]
#    stat_filt=mystation.split("_")[1]
#
#    columns="SUM("+mystation+")"
#    selstring="SELECT "+columns+" FROM allist"
#    with con:
#            cur.execute(selstring)
#    reply=cur.fetchall()
#    onesampl=reply[0][0]
#
#    columns="SUM("+mystation+")"
#    selstring="SELECT "+columns+" FROM allist WHERE hit_abb!='nod'"
#    with con:
#            cur.execute(selstring)
#    reply=cur.fetchall()
#    hitcount=reply[0][0]
#
#    stationcsv+=mystation+","+stat_loc+","+stat_filt+","+str(onesampl)+","+str(hitcount)+","+getfromcc(mystation.split("_")[0])+"\n"
##
#open(args.outcsv,"w").write(summarycsv)
#fhh=open(args.outcsv+"stats.csv","w")
#fhh.write(stationcsv)
#fhh.close()


station_dataset=robjects.r('statdataset <- read.delim("'+args.outcsv+"stats.csv"+'", sep="," , header=TRUE)')
print station_dataset

robjects.r('attach(statdataset)')

robjects.r('library(reshape2)')
robjects.r('library(RColorBrewer)')
robjects.r('library(stats)')

grdevices = importr('grDevices')




# Make stationsdb
#################
stations_db_place=basepicdir+"stationsdb.db"
con4 = sqlite3.connect(stations_db_place)
cur4 = con4.cursor()

dropstring="DROP TABLE IF EXISTS "+"t"
cur4.execute(dropstring)
cur4.execute("CREATE TABLE t (name,mystation,filtersize,n_hit_station,alltranscripts,Lon,Lat,O,Sampledepth);")

with open(args.outcsv+"stats.csv",'r') as fin: # `with` statement available in 2.5+
    # csv.DictReader uses first line in file for column headings by default
    dr = csv.DictReader(fin) # comma is default delimiter
    to_db = [(i['name'],i['mystation'],i['filtersize'], i['n_hit_station'],i['alltranscripts'],i['Lon'], i['Lat'],i['O'],i['Sampledepth']) for i in dr]

cur4.executemany("INSERT INTO t (name,mystation,filtersize,n_hit_station,alltranscripts, Lon,Lat,O,Sampledepth) VALUES (?,?,?, ?,?,?, ?,?,?);", to_db)
con4.commit()


# Made stations db
############

#make csv for all depths, added filtersizes
have_lats_filtadd_csv=basepicdir+"have_lats_filtersadd.csv"
data = cur4.execute("SELECT name,mystation, sum(n_hit_station),sum(alltranscripts), Lon, Lat, Sampledepth FROM t WHERE Lat!='nod2' GROUP BY mystation")
with open(have_lats_filtadd_csv, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['name','mystation','n_hit_station','alltranscripts','Lon', 'Lat','Sampledepth'])
    writer.writerows(data)

graphtitle="Total number of reads"
#bar_allsize(have_lats_filtadd_csv,basepicdir+"allbars_alltranscripts.pdf",10,10,graphtitle,"alltranscripts")
graphtitle="Total number of rpod"
#bar_allsize(have_lats_filtadd_csv,basepicdir+"allbars_rpods.pdf",10,10,graphtitle,"n_hit_station")
graphtitle="n_hit_station_alltranscripts"
#bar_allsize(have_lats_filtadd_csv,basepicdir+"allbars_rpods_div_alltranscripts.pdf",10,10,graphtitle,"n_hit_station/alltranscripts")
#make csv for all depths, all filtersizes
have_lats_csv=basepicdir+"have_lats.csv"
data = cur4.execute("SELECT name,mystation, filtersize,n_hit_station,alltranscripts, Lon, Lat, Sampledepth FROM t WHERE Lat!='nod2'")
with open(have_lats_csv, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['name','mystation','filtersize','n_hit_station','alltranscripts','Lon', 'Lat','Sampledepth'])
    writer.writerows(data)

#graphtitle="Total number of reads"
#bar_allsize(have_lats_csv,basepicdir+"allbars_alltranscripts.pdf",10,10,graphtitle,"alltranscripts")
#graphtitle="Total number of rpod"
#bar_allsize(have_lats_csv,basepicdir+"allbars_rpods.pdf",10,10,graphtitle,"n_hit_station")
#graphtitle="n_hit_station_alltranscripts"
#bar_allsize(have_lats_csv,basepicdir+"allbars_rpods_div_alltranscripts.pdf",10,10,graphtitle,"n_hit_station/alltranscripts")
#
# Collect data with filter sizes
#data = cur4.execute("SELECT name,mystation, alltranscripts, Lon, Lat, Sampledepth FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)<=0.30487806")

surfLats_csv=basepicdir+"surf_lats.csv"
# Make database with latitudes with sampledepth < 0.3
#######################
#data = cur4.execute("SELECT * FROM t")
#data = cur4.execute("SELECT * FROM t WHERE Lat!='nod2'")
data = cur4.execute("SELECT name,mystation, SUM(n_hit_station),SUM(alltranscripts), Lon, Lat, Sampledepth FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)<=0.30487806 GROUP BY mystation")
with open(surfLats_csv, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['name','mystation','n_hit_station','alltranscripts','Lon', 'Lat','Sampledepth'])
    writer.writerows(data)

    
surfLats_csv2=basepicdir+"surf_lats2.csv"
# Make database with latitudes with sampledepth < 0.3
#######################
#data = cur4.execute("SELECT * FROM t")
#data = cur4.execute("SELECT * FROM t WHERE Lat!='nod2'")
data = cur4.execute("SELECT name,mystation, SUM(n_hit_station),SUM(alltranscripts), Lon, Lat, Sampledepth FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)<=0.30487806 GROUP BY mystation")
with open(surfLats_csv2, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['name','mystation','n_hit_station','alltranscripts','Lon', 'Lat','Sampledepth'])
    writer.writerows(data)


# write surface maps
#gogmap_google(surfLats_csv,basepicdir+"surface_hybrid_map.pdf",10,10)
#with_rworld(surfLats_csv,basepicdir+"surface_rworld.pdf",10,10)
#gogmap_google('/Users/security/science/output.csv',"/Users/security/science/surfacemap.pdf")


midLats_csv=basepicdir+"mid_lats.csv"
# Make database with latitudes with sampledepth < middle
#######################
data = cur4.execute("SELECT name,mystation, sum(n_hit_station),SUM(alltranscripts), Lon, Lat, Sampledepth FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)>0.30487806 AND cast(Sampledepth as float)<20 GROUP BY mystation")
with open(midLats_csv, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['name','mystation','n_hit_station','alltranscripts','Lon', 'Lat','Sampledepth'])
    writer.writerows(data)

#gogmap_google(midLats_csv,basepicdir+"middle_hybrid_map.pdf",10,10)
#with_rworld(midLats_csv,basepicdir+"middlemap_rworld.pdf",10,10)



#deepLats_csv=basepicdir+"deep_lats.csv"
## Deep lvl
#################
#data = cur4.execute("SELECT name,mystation, SUM(alltranscripts), Lon, Lat, Sampledepth FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)>20 AND cast(Sampledepth as float)<200 GROUP BY mystation")
#with open(deepLats_csv, 'w') as f:
#    writer = csv.writer(f)
#    writer.writerow(['name','mystation','alltranscripts','Lon', 'Lat','Sampledepth'])
#    writer.writerows(data)


#gogmap_google(deepLats_csv,basepicdir+"deep_hybrid_map.pdf",10,10)
#with_rworld (deepLats_csv,basepicdir+"deep_rworld_map.pdf",10,10)


gogmap3(surfLats_csv,midLats_csv,basepicdir+"deep_map_multi.pdf",10,10)




#for filtersize in ["0p1","0p8","3p0"]:
#
#    # take only one filtersize
#    o3p0=robjects.r('only_filter<-statdataset[grep("'+filtersize+'", statdataset$filtersize), ]')
#
## remove landsort
#    robjects.r('only_filter_nols <- only_filter[! grepl("GS678",only_filter$mystation),]')
#
## only landsort
#    #robjects.r('onlyls <- only_depth[grepl("GS678",only_depth$mystation),]')
#   
#    # melt the data
#    robjects.r('molten_stations<-melt(only_filter_nols)')
#
#    robjects.r('cols <- colorRampPalette(brewer.pal(9, "Accent"))')
#
#
#    allISfiltersgraph=robjects.r('allplot<-ggplot(data=molten_stations, aes(x=molten_stations$mystation, y=molten_stations$value,fill=molten_stations$variable))') 
#    #allISfiltersgraph=robjects.r('allplot<-ggplot(data=molten_stations, aes(x=molten_stations$mystation, y=molten_stations$value,fill=molten_stations$variable))') 
#    allISfiltersgraph=robjects.r('allplot<-allplot+geom_bar(color="black",stat="identity",position="fill",size=0.2)')
#
#    allISfiltersgraph=robjects.r('allplot<-allplot+theme(axis.text.x=element_text(angle=90))')
#    #allISfiltersgraph=robjects.r('allplot<-allplot+theme(axis.text.x=element_text(angle=90),axis.text.x = element_blank())')
#
#    allISfiltersgraph=robjects.r('allplot<-allplot+scale_fill_manual(values=cols(19))')
#    gols=robjects.r('scale_fill_manual(values=cols(19))')
#    allISfiltersgraph=robjects.r('allplot<-allplot+guides(fill = guide_legend(reverse=TRUE,title = "Group", title.position = "top"))')
#    allISfiltersgraph=robjects.r('allplot<-allplot+scale_y_continuous(breaks=seq(0, 1, 0.1),name="Fraction of transcripts")')
#    allISfiltersgraph=robjects.r('allplot<-allplot+scale_x_discrete(name="Station")')
#
#    allISfiltersgraph=robjects.r('allplot<-allplot+labs(title="'+"fractions_"+filtersize+'")')
#
#
#    robjects.r.ggsave("/Users/security/science/"+filtersize+"siga_fractions2.pdf")
#    #robjects.r.ggsave(args.out.replace(".pdf","")+depth+".pdf")
#    

#broup="SELECT DISTINCT group_ FROM allist"
#with con:
#        cur.execute(broup)
#orgroups=cur.fetchall()
##
##
#dets=heads[33:87]
##print dets
##
#had="mystation,"
#for i in orgroups:
#    print "text"+i[0]+"tes"
#    hed=i[0]
#    if len(hed)<5:
#        hed="unknown"
#    had+=hed.replace(" ","").replace("/","")+","

#raw_input()
#
#had=had[:-1]
#had+="\n"
#outstring=had

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

#output="/Users/security/science/fractions.csv"
#output="/Users/security/science/output4.csv"
#csvWriter = open(output, "w")
#csvWriter.write(outstring)

#
#
#robjects.r('library(reshape2)')
#robjects.r('library(RColorBrewer)')
#robjects.r('library(stats)')
##
###dataf=robjects.r('mydataset <- read.delim("'+args.in_csv+'", sep="," , header=TRUE)')
#dataf=robjects.r('mydataset <- read.delim("'+output+'", sep="," , header=TRUE)')
#robjects.r('attach(mydataset)')
#
#
#
#for depth in ["0p1","0p8","3p0"]:
#
#
#
#    o3p0=robjects.r('only_depth<-mydataset[grep("'+depth+'", mydataset$mystation), ]')
#
## remove landsort
#    robjects.r('new_d <- only_depth[! grepl("GS678",only_depth$mystation),]')
#
## only landsort
#    #robjects.r('onlyls <- only_depth[grepl("GS678",only_depth$mystation),]')
#    
#    robjects.r('testtt<-melt(new_d)')
#
#    robjects.r('cols <- colorRampPalette(brewer.pal(9, "Accent"))')
#    robjects.r('print(cols)')
#
#    allISfiltersgraph=robjects.r('allplot<-ggplot(data=testtt, aes(x=testtt$mystation, y=testtt$value,fill=testtt$variable))') 
#    allISfiltersgraph=robjects.r('allplot<-allplot+geom_bar(color="black",stat="identity",position="fill",size=0.2)')
#
#    allISfiltersgraph=robjects.r('allplot<-allplot+theme(axis.text.x=element_text(angle=90))')
#    #allISfiltersgraph=robjects.r('allplot<-allplot+theme(axis.text.x=element_text(angle=90),axis.text.x = element_blank())')
#
#    allISfiltersgraph=robjects.r('allplot<-allplot+scale_fill_manual(values=cols(19))')
#    gols=robjects.r('scale_fill_manual(values=cols(19))')
#    allISfiltersgraph=robjects.r('allplot<-allplot+guides(fill = guide_legend(reverse=TRUE,title = "Group", title.position = "top"))')
#    allISfiltersgraph=robjects.r('allplot<-allplot+scale_y_continuous(breaks=seq(0, 1, 0.1),name="Fraction of transcripts")')
#    allISfiltersgraph=robjects.r('allplot<-allplot+scale_x_discrete(name="Station")')
#
#    allISfiltersgraph=robjects.r('allplot<-allplot+labs(title="'+"fractions_"+depth+'")')
#
#
#    robjects.r.ggsave("/Users/security/science/"+depth+"siga_fractions.pdf")
#    #robjects.r.ggsave(args.out.replace(".pdf","")+depth+".pdf")



