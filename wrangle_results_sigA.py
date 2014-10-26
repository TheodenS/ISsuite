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
def gogmap_google(incsv,outname):

    robjects.r('library(ggmap)')
    robjects.r('map <- get_googlemap(maptype ="roadmap",center = c(lon=18,lat=61.5),zoom=4,color="bw")')

    station_dataset=robjects.r('chosendataset <- read.delim("'+'/Users/security/science/output.csv'+'", sep="," , header=TRUE)')
    #station_dataset=robjects.r('chosendataset <- read.delim("'+'/Users/security/science/output.csv'+'", sep="," , header=TRUE)')

    robjects.r('pdf("'+outname+'")')
    robjects.r('oio <- ggmap(map) + geom_point(data=chosendataset,color="red",aes(x=Lon, y = Lat,size=alltranscripts))')
    robjects.r('oio <- oio + geom_text(data=chosendataset,size=3,color="red",aes(x=Lon,y=Lat,label=mystation),hjust=0.5, vjust=-0.7)')
    im4=robjects.r('print(oio)')
    robjects.r('dev.off()')

def gogmap_ggmap(incsv,outname):

    robjects.r('library(ggmap)')
    robjects.r('library(rgdal)')
    #robjects.r('map.world <- map_data (map = "world")')
    robjects.r('map.world<-map.world[-which(getMap()$ADMIN== "Antarctica"),]')
    robjects.r('map.world <- spTransform(map.world, CRS=CRS("+proj=robin +ellps= WGS84"))')

    station_dataset=robjects.r('chosendataset <- read.delim("'+'/Users/security/science/output.csv'+'", sep="," , header=TRUE)')
    #station_dataset=robjects.r('chosendataset <- read.delim("'+'/Users/security/science/output.csv'+'", sep="," , header=TRUE)')

    robjects.r('pdf("'+outname+'")')
    robjects.r('p1 <- ggplot(map.world, aes (x = long, y = lat, group = map.world$group,colour="black"))')
    robjects.r('p1 <- p1 + geom_path (colour="black")')
    robjects.r('p1 <- p1 + geom_point(data=chosendataset,color="red",aes(x=Lon, y = Lat,size=alltranscripts))')
    robjects.r('p1 <- p1 + theme (legend.position= "none")')
    robjects.r('p1 <- p1 + xlim(-10,50)')
    robjects.r('p1 <- p1 + ylim(50,70)')
    im4=robjects.r('print(p1)')
    robjects.r('dev.off()')

def gogmap3(incsv,outname):



    robjects.r('library(grid)')
    robjects.r('library(gridExtra)')
    robjects.r('pushViewport(viewport(layout = grid.layout(1, 2)))')
    
    robjects.r('library(rworldmap)')
    robjects.r('worldMap <- getMap(resolution="high")')

    robjects.r('world.points <- fortify(worldMap)')
    robjects.r('world.points$region <- world.points$id')
    robjects.r('world.df <- world.points[,c("long","lat","group", "region")]')
    robjects.r('geom_polygon(data = world.df, aes(x = long, y = lat, group = group))')
    robjects.r('m <- ggplot() ')

    robjects.r('m<-m +geom_path(data = world.df, aes(x = long, y = lat, group = group))')
    #robjects.r('m<-m +scale_y_continuous(breaks = (-2:2) * 30)')
    #robjects.r('m<-m +scale_x_continuous(breaks = (-4:4) * 45)')
    #robjects.r('m<-m +coord_map("ortho", orientation=c(51,40, 0),xlim=c(-10,50),ylim=c(50,70))')
    robjects.r('m<-m +coord_map("ortho", orientation=c(10,16, 0),xlim=c(8,24),ylim=c(54,70))')
    #robjects.r('m<-m +coord_map("ortho", orientation=c(51,40, 0),xlim=10,ylim=40)')
    robjects.r('m<-m +geom_point(data=chosendataset,color="red",aes(x=Lon, y = Lat,size=1))')
    #robjects.r('m<-m +geom_point(data=chosendataset,color="red",aes(x=Lon, y = Lat,size=chosendataset$alltranscripts))')
    #im4=robjects.r('g<-multiplot(m,m,col=2)')
    #robjects.r('m <- m + xlim(-10,50)')
    #robjects.r('m <- m + ylim(50,70)')
    #objects.r('geom_point(data = mypoints, aes(x = long, y = lat), color = "blue3")')

    #robjects.r('m<-m + xlim(0, 37) + ylim(53, 70)')
    #robjects.r('m<-ggplot() + coord_map(xlim = c(-169, -121), ylim = c(45, 70)) + geom_polygon(data = mymap, aes(long, lat, group = group), color = "grey20", fill = "grey15", size = 0.3) + geom_point(data = mypoints, aes(x = long, y = lat), color = "blue3") + xlim(-179, -90) + ylim(30, 80)')

   
    robjects.r('pdf("'+outname+'")')
    robjects.r('grid.arrange( m, m, ncol=2)')
    #robjects.r('m <- m + geom_point(data=chosendataset,color="red",aes(x=Lon, y = Lat,size=alltranscripts))')

    #robjects.r('oio <- oio + geom_text(data=chosendataset,size=3,color="red",aes(x=Lon,y=Lat,label=mystation),hjust=0.5, vjust=-0.7)')
    im4=robjects.r('print(m)')
    #im4=robjects.r('print(m, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))')
    robjects.r('dev.off()')



def gogmap2(incsv,outname):


    robjects.r('library(rworldmap)')


    #newmap <- getMap(resolution = "low")

    #robjects.r('newmap <- getMap(resolution = "low")')


    #robjects.r('library(rgdal)')
    #robjects.r('map.world <- map_data (map = "world")')
    #robjects.r('map.world<-map.world[-which(getMap()$ADMIN== "Antarctica"),]')



    robjects.r('data("countryExData", envir=environment(), package="rworldmap",asp=1.25)')
    robjects.r('mymap <- joinCountryData2Map(countryExData, joinCode = "ISO3", nameJoinColumn = "ISO3V10", mapResolution = "high",projection=1.25)')
    #robjects.r('mymap <- mapCountryData(mymap, aspect=1.25)')
    #robjects.r('mymap <- spTransform(mymap, CRS=CRS("+proj=robin +ellps=WGS84"))')
    robjects.r('mymap <- fortify(mymap)') 
    #robjects.r('mypoints <- data.frame(lat = rep(55, 3), long = c(-145, -147, -149))')

    robjects.r('m<-ggplot()')
    #robjects.r('m<-m+ coord_map(xlim = c(-10, 50), ylim = c(50, 70))')
    robjects.r('m<-m + geom_path(data = mymap, aes(mymap$long, mymap$lat, group = mymap$group), color = "grey20", fill = "grey15", size = 0.3)')

    robjects.r('m<-m +geom_point(data=chosendataset,color="red",aes(x=Lon, y = Lat))')
    #robjects.r('m<-m +geom_point(data=chosendataset,color="red",aes(x=Lon, y = Lat,size=chosendataset$alltranscripts))')
    #objects.r('geom_point(data = mypoints, aes(x = long, y = lat), color = "blue3")')

    robjects.r('m<-m + xlim(0, 37) + ylim(53, 70)')
    #robjects.r('m<-ggplot() + coord_map(xlim = c(-169, -121), ylim = c(45, 70)) + geom_polygon(data = mymap, aes(long, lat, group = group), color = "grey20", fill = "grey15", size = 0.3) + geom_point(data = mypoints, aes(x = long, y = lat), color = "blue3") + xlim(-179, -90) + ylim(30, 80)')

   
    robjects.r('pdf("'+outname+'")')
    #robjects.r('m <- m + geom_point(data=chosendataset,color="red",aes(x=Lon, y = Lat,size=alltranscripts))')

    #robjects.r('oio <- oio + geom_text(data=chosendataset,size=3,color="red",aes(x=Lon,y=Lat,label=mystation),hjust=0.5, vjust=-0.7)')
    im4=robjects.r('print(m)')
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

args=parser.parse_args()

metadata="/Users/security/science/RNA/mimebs_metadata/metadata.2009.csv"

metafh=open(args.environmentaldata,"r")
metaread=metafh.read()
metalines=metaread.split("\n")
metahead=metalines[0]
metalines=metalines[1:]



	
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



#db="/Users/security/science/dest3.db"
db=args.in_db
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
############
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


#robjects.r('pdf("'+"/Users/security/science/testmap2.pdf"+'",width=10,heigh=10)')


# pick only the rows with long values

#ppg=robjects.r('usedataset <- statdataset[! grepl("nod2",statdataset$Lon),]')
ppg=robjects.r('usedataset <- statdataset')
print "dataset"
print ppg

onlysurf=robjects.r('onlysurface<-usedataset[ which("Sampledepth" =="0.30487806"),]')

robjects.r('mylocation <- c(lon=18,lat=61)')
#print"longs py print"
#lo=robjects.r('lol<-usedataset$Lon')
#la=robjects.r('lot<-usedataset$Lat')
#print lo
#print la
#
#print"longs r print"
#robjects.r('print(usedataset$Lon)')
#longlen=robjects.r('xxx<-length(usedataset$Lon)')
#latlen=robjects.r('xxy<-length(usedataset$Lat)')
#framelen=robjects.r('framelen<-length(usedataset)')
#print longlen
#print latlen
#print framelen
#
##robjects.r('print(usedataset$Lon)')
#robjects.r('library(mapproj)')
##robjects.r('map <- get_map(source="google",maptype ="roadmap",location = mylocation, zoom = 4)')
##robjects.r('map <- get_map(source="google",maptype ="satellite",location = mylocation, zoom = 4)')
##robjects.r('map <- get_map(source="google",maptype ="terrain",location = mylocation, zoom = 4)')
#robjects.r('map <- get_googlemap(maptype ="satellite",center = c(lon=18,lat=61.5),zoom=4,color="bw")')






#robjects.r('oio<-ggmap(map)')
#robjects.r('map <- get_map(location = "Europe", zoom = 4)')




con4 = sqlite3.connect("/Users/security/science/stationsdb.db")
cur4 = con4.cursor()

dropstring="DROP TABLE IF EXISTS "+"t"
cur4.execute(dropstring)
cur4.execute("CREATE TABLE t (name,mystation,alltranscripts,Lon,Lat,Sampledepth);")

with open(args.outcsv+"stats.csv",'r') as fin: # `with` statement available in 2.5+
    # csv.DictReader uses first line in file for column headings by default
    dr = csv.DictReader(fin) # comma is default delimiter
    to_db = [(i['name'],i['mystation'],i['alltranscripts'] , i['Lon'] ,i['Lat'],i['Sampledepth']) for i in dr]

cur4.executemany("INSERT INTO t (name,mystation,alltranscripts, Lon,Lat,Sampledepth) VALUES (?,?,?, ?,?,?);", to_db)
con4.commit()

#data = cur4.execute("SELECT * FROM t")
#data = cur4.execute("SELECT * FROM t WHERE Lat!='nod2'")
data = cur4.execute("SELECT name,mystation, SUM(alltranscripts), Lon, Lat, Sampledepth FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)<=0.30487806 GROUP BY mystation")

with open('/Users/security/science/output.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['name','mystation','alltranscripts','Lon', 'Lat','Sampledepth'])
    writer.writerows(data)


gogmap_google(args.outcsv+"stats.csv","/Users/security/science/surfacemap.pdf")



data = cur4.execute("SELECT name,mystation, SUM(alltranscripts), Lon, Lat, Sampledepth FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)>0.30487806 AND cast(Sampledepth as float)<20 GROUP BY mystation")

with open('/Users/security/science/output.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['name','mystation','alltranscripts','Lon', 'Lat','Sampledepth'])
    writer.writerows(data)

gogmap_google(args.outcsv+"stats.csv","/Users/security/science/middlemap.pdf")

data = cur4.execute("SELECT name,mystation, SUM(alltranscripts), Lon, Lat, Sampledepth FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)>20 AND cast(Sampledepth as float)<200 GROUP BY mystation")

with open('/Users/security/science/output.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['name','mystation','alltranscripts','Lon', 'Lat','Sampledepth'])
    writer.writerows(data)

gogmap3(args.outcsv+"stats.csv","/Users/security/science/deepmap.pdf")
#station_dataset=robjects.r('chosendataset <- read.delim("'+'/Users/security/science/output.csv'+'", sep="," , header=TRUE)')
#nonsurf=robjects.r('nosurf<-chosendataset[which(chosendataset$Sampledepth !=0.30487806),]')
#onlysurf=robjects.r('chosendataset<-chosendataset[which(chosendataset$Sampledepth ==0.30487806),]')
##onlysurf0=robjects.r('chosendataset0<-chosendataset[which(chosendataset$Sampledepth >0.30487806),]')
#print onlysurf
#
#robjects.r('s <- map_data("world")')
#
#robjects.r('m <- ggplot(s, aes(x=long, y=lat, group=group)) + geom_polygon(fill="white", colour="black")')
#robjects.r('print(m)')
#
#
#
#robjects.r('pdf("'+"/Users/security/science/testmap3.pdf"+'")')
#
#robjects.r('library(rworldmap)') 
#robjects.r('newmap <- getMap(resolution = "high")') 
#robjects.r('nmap<-plot(newmap, xlim = c(5, 30), ylim = c(53, 69), asp = 1)')
#robjects.r('points(18, 62, col = "red", cex = .6)')	
#robjects.r('dev.off()')
#
#robjects.r('pdf("'+"/Users/security/science/testmap2.pdf"+'")')
#robjects.r('oio <- ggmap(map) + geom_point(data=chosendataset,color="red",aes(x=Lon, y = Lat,size=alltranscripts))')
##robjects.r('oio <- oio + aes_string(x=Lon, y = Lat,label=name)')
##robjects.r('oio <- oio + geom_text(data=chosendataset,aes(x=Lon,y=Lat,label=mystation,size=60),hjust=0, vjust=0)')
#robjects.r('oio <- oio + geom_text(data=chosendataset,size=3,color="red",aes(x=Lon,y=Lat,label=mystation),hjust=0.5, vjust=-0.7)')
#robjects.r('oio <- oio + geom_text(data=nosurf,size=3,color="blue",aes(x=Lon,y=Lat,label=mystation),hjust=0.5, vjust=0.7)')
#
#
##robjects.r('oio <- ggmap(map) + geom_point(aes(x = Lon, y = Lat, size = 20,data=usedataset,color="red"))')
##robjects.r('oio <- ggmap(map) + geom_point(aes(x = 18, y = 62, size = 200), alpha = .5)')
#im4=robjects.r('print(oio)')
#
#
#robjects.r('dev.off()')

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



