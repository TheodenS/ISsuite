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

def singleplot(yvar):
    # NH4
    ####################
    pointsize=str(1)
    axistextsize=str(5)
    titlesize=str(6)
    yvar_in="modeldata$"+yvar
    modelmkstring=yvar+'plot<-ggplot(data=modeldata, aes(x='+yvar_in+', y=fraction,fill="green"))'
    modelmkstring+='+ geom_point(data=modeldata,aes(x='+yvar_in+', y=fraction,color="blue", size='+pointsize+'))'
    modelmkstring+='+stat_smooth(data=modeldata,method = "lm", se=TRUE, color="black", formula = y ~ x)'
    modelmkstring+='+scale_x_continuous("'+yvar+'")'
    modelmkstring+='+ggtitle("'+yvar+'")'
    #robjects.r('corstats<-cor.test('+yvar_in+' , modeldata$fraction,model="person")')

    modelmkstring+='+theme('
    modelmkstring+='axis.text.x=element_text(size='+axistextsize+',angle=0)'
    modelmkstring+=',axis.text.y=element_text(size='+axistextsize+',angle=0)'
    modelmkstring+=',axis.title.y=element_text(size='+axistextsize+',angle=90)'
    modelmkstring+=',axis.title.x=element_text(size='+axistextsize+',angle=0)'
    modelmkstring+=',legend.position="none"'
    modelmkstring+=',plot.title = element_text(size='+titlesize+',lineheight=1, face="bold")'
    modelmkstring+=')'
    #robjects.r(modelmkstring) 

    return modelmkstring

#arrangeGrob(bargall_mid_frac,bargall_mid_total_transcripts,bargall_mid_hits,ncol=1,heights=c(2,1,12,1,1))

    
    #'oio2<-mp <- mp+ geom_point(aes(x=visit.x, y=visit.y) ,color="blue", size=3) 
#teest=robjects.r('p<-p+coord_fixed(ratio=400)')

def make_blank_window(statvar):
    textsize=str(4)
    make_text_plot='+scale_x_continuous(limits=c(1, 100))'
    make_text_plot+='+scale_y_continuous(limits=c(1, 100))'

    make_text_plot+='+geom_text(aes(x = 3, y = 85, size='+textsize+',label = paste("p-value: ",round('+statvar+'$p.value,3))), data = textdf, colour = I("black"),hjust=0)'
    make_text_plot+='+geom_text(aes(x = 3, y = 45, size='+textsize+',label = paste("rho: ",round('+statvar+'$estimate,3))), data = textdf, colour = I("black"),hjust=0)'

    make_text_plot+='+theme('
    make_text_plot+='legend.position="none"'
    make_text_plot+=',axis.title.x=element_blank()'
    make_text_plot+=',axis.title.y=element_blank()'
    make_text_plot+=',axis.text.x=element_blank()'
    make_text_plot+=',axis.text.y=element_blank()'
    make_text_plot+=',axis.ticks = element_blank()'
    make_text_plot+=',panel.grid.minor=element_blank()'
    make_text_plot+=',panel.grid.major=element_blank()'
    make_text_plot+=',panel.border=element_blank()'
    make_text_plot+=',panel.background = element_rect(fill="white", colour="black")'
    make_text_plot+=')'
    return make_text_plot


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

    bar_settings+='+ggtitle("'+title+mycolor+'")' 
    bar_settings+='+theme('
    bar_settings+='plot.title = element_text(size='+plottitlesize+',lineheight=0.8)'
    bar_settings+=',axis.text.x=element_text(size='+axistextsize+',angle=0)'
    bar_settings+=',axis.text.y=element_text(size='+axistextsize+',angle=0)'
    bar_settings+=',axis.title.x=element_blank()'
    bar_settings+=',axis.title.y=element_text(size='+axistitlesize+')'
    bar_settings+=',legend.position="none"'
    bar_settings+=')'

    #bar_settings+='+geom_text(size=3,color="red",aes(x=0,y=0,label="'+mycolor+'"))'
    bar_settings+='+scale_x_discrete(name="'+xname+'")'     
    bar_settings+=''     
    bar_settings+='+scale_y_continuous(labels=comma, name="'+yname+'")'     

#teest=robjects.r('p<-p+coord_fixed(ratio=400)')
    return bar_settings

def get_line_settings(dataset,title,yname):
#graphtitle="TESTTI"

    mycolorfirst = "%06x" % random.randint(0,0xFFFFFF)
    mycolor="#"+str(mycolorfirst).upper()
    print mycolor
    plottitlesize=str(6)
    axistextsize=str(8)
    axistitlesize=str(4)
    legendsize=str(4)
    xname="station"

    robjects.r('library(scales)')
    env_vars=['pH','NH4','Chla','P']
    line_settings='geom_line(data='+dataset+',aes(x=mystation, y=scale(Temp),group=2,color="Temp"))'
    line_settings+='+geom_line(data='+dataset+',aes(x=mystation, y=scale(Sal),group=2,color="Sal"))'
    for env_var in env_vars:
        line_settings+='+geom_line(data='+dataset+',aes(x=mystation, y=scale('+env_var+'),group=2,color="'+env_var+'"))'

    #bar_settings+='+ggtitle("'+title+mycolor+'")' 
    line_settings+='+theme('
    line_settings+='plot.title = element_text(size='+plottitlesize+',lineheight=0.8)'
    line_settings+=',axis.text.x=element_text(size='+axistextsize+',angle=0)'
    line_settings+=',axis.text.y=element_text(size='+axistextsize+',angle=0)'
    line_settings+=',axis.title.y=element_blank()'
    line_settings+=',axis.title.x=element_blank()'
    line_settings+=',legend.key.size=unit(0.7,"lines")'
    line_settings+=',legend.title=element_blank()'
    line_settings+=',legend.position="bottom"'
    line_settings+=')'

    line_settings+='+scale_x_discrete(name="'+xname+'")'     
    #bar_settings+=''     
    #bar_settings+='+scale_y_continuous(labels=comma, name="'+yname+'")'     

#teest=robjects.r('p<-p+coord_fixed(ratio=400)')
    return line_settings

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
#teest=robjects.r('p<-p+coord_fixed(ratio=400)')
    robjects.r('oio2 <- oio2+ggtitle("'+graphtitle+'")')
    robjects.r('oio2 <- oio2+theme(plot.title = element_text(lineheight=.8, face="bold"))')
    #robjects.r.ggsave("/Users/security/science/ttest.pdf")
    robjects.r('pdf("'+outname+'", width='+str(w)+',height='+str(h)+')')
    im4=robjects.r('print(oio2)')
    robjects.r('dev.off()')

def make_model(incsv,outname,w,h,graphtitle,plotvar):
    robjects.r('library(gridExtra)')
    robjects.r('library(ggmap)')
    station_dataset=robjects.r('modeldata <- read.delim("'+incsv+'", sep="," , header=TRUE)')
    
    #robjects.r('modeldata$fraction<-modeldata$n_hit_station/modeldata$alltranscripts')
    robjects.r('modeldata$fraction<-modeldata$n_hit_station/modeldata$sig_a_alltranscripts')

    envarlist=["alltranscripts","NH4","Temp","Urea","NO2","NP","PO4","N","P","NP","DOP","pH","Sal"]
    for envar in envarlist:
        yvar_in="modeldata$"+envar
        robjects.r(envar+'stats<-cor.test(modeldata$fraction,'+yvar_in+',model="spearman")')
        robjects.r(singleplot(envar))
        robjects.r('xvars<-1:100')
        robjects.r('yvars<-1:100')
        robjects.r('textdf<-data.frame(xvars,yvars)')
        make_text_plot=envar+'Tplot<-ggplot(data=textdf, aes(x=xvars, y=yvars))'     
        make_text_plot+=make_blank_window(envar+"stats")
        robjects.r(make_text_plot)

    #robjects.r.ggsave("/Users/security/science/ttest.pdf")
    robjects.r('pdf("'+outname+'", width='+str(w)+',height='+str(h)+')')

    gridline='grid.arrange('

    for envar in envarlist:
        gridline+='arrangeGrob('+envar+'plot,'+envar+'Tplot,ncol=1,heights=c(2,1)),'
    #gridline=gridline[:-1]
    gridline+='ncol=2)'
    print gridline
    robjects.r(gridline)

    robjects.r('dev.off()')


def makepairs(incsv,outname,w,h,graphtitle,plotvar):
    gal = importr('GGally')

    robjects.r('png("'+outname+'",width=2000,heigh=2000)')
    station_dataset=robjects.r('chosendataset <- read.delim("'+incsv+'", sep="," , header=TRUE)')
    robjects.r('print(chosendataset)')
    #robjects.r('keeps <- c("n_hit_station","alltranscripts","Chla","O","Temp","Sal","pH","N","NO3","NO2","NH4","DOP","NP","P","N","PO4","Urea","Si","Sampledepth","siga_alltranscripts","sig_a_n_hit_station")')

    robjects.r('keeps <- c("n_hit_station","alltranscripts","Chla","O","Temp","Sal","pH","N","NO3","NO2","NH4","DOP","NP","P","N","PO4","Urea","Si","Sampledepth","sig_a_n_hit_station","pro_c_n_hit_station","cpc_g_n_hit_station")')

    robjects.r('factors<-chosendataset[keeps]')
    robjects.r('factors$fraction<-(factors$n_hit_station)/(factors$alltranscripts)')
    #robjects.r('factors$fraction<-factors$n_hit_station/factors$alltranscripts')
    robjects.r('factors$fractionofsiga<-factors$n_hit_station/factors$sig_a_n_hit_station')
    robjects.r('factors$fractionofproc<-factors$n_hit_station/factors$pro_c_n_hit_station')
    #robjects.r('factors$fractionofcpcg<-factors$n_hit_station/factors$cpc_g_n_hit_station')
    robjects.r('print(factors)')
    im4=robjects.r('environpairspic<-ggpairs(factors,diag=list(labelSize=3),upper = list(params = c(size = 4)),lower = list(continuous = "smooth", size=0.01,params = c(colour="red",method = "lm")))')
    envpairs=robjects.r('print(environpairspic)')
#robjects.r('ggsave("/Users/security/science/new3.pdf")')
    robjects.r('dev.off()')





def bars(incsv,outname,w,h,graphtitle,plotvar):
    #graphtitle="Total number of reads"
#robjects.r('mylocation <- c(lon=18,lat=61)')
    robjects.r('library(ggmap)')
    robjects.r('library(reshape2)')
    robjects.r('library(RColorBrewer)')
    robjects.r('library(gridExtra)')

    depth_title="surfaceORdeep"

    print incsv
    station_dataset=robjects.r('chosendataset <- read.delim("'+incsv+'", sep="," , header=TRUE)')

    robjects.r('xvars<-1:100')
    robjects.r('yvars<-1:100')
    robjects.r('plotdf<-data.frame(xvars,yvars)')
    #make_text_plot='textplot<-ggplot(data=plotdf, aes(x=xvars, y=yvars))'     
    make_text_plot='textplot<-ggplot(plotdf)'     


    # Maximum computations
    ##########
    maxhitcount=robjects.r('max_hit_count<-max(chosendataset$n_hit_station)')
    whichmaxhitcount=robjects.r('ind_max_hit_count<-which.max(chosendataset$n_hit_station)')
    maxhitstation=robjects.r('max_hit_count_name<-chosendataset[ind_max_hit_count,2]')
    #maxhitstation=robjects.r('x[which.min(x[,2]),1]')

    maxreadcount=robjects.r('max_read_count<-max(chosendataset$alltranscripts)')
    whichminhitcount=robjects.r('ind_max_read_count<-which.max(chosendataset$alltranscripts)')
    maxreadstation=robjects.r('max_read_count_name<-chosendataset[ind_max_read_count,2]')

    # fractions
    fraccount=robjects.r('chosendataset$fraction<-chosendataset$n_hit_station/chosendataset$alltranscripts')
    maxfrac=robjects.r('print(chosendataset)')
    maxfrac=robjects.r('max_frac<-max(chosendataset$fraction)')
    whichmaxfraccount=robjects.r('ind_max_frac<-which.max(chosendataset$fraction)')
    maxfracstation=robjects.r('max_frac_name<-chosendataset[ind_max_frac,2]')

    # Minimum computations
    ##########
    minhitcount=robjects.r('min_hit_count<-min(chosendataset$n_hit_station)')
    whichminhitcount=robjects.r('ind_min_hit_count<-which.min(chosendataset$n_hit_station)')
    minhitstation=robjects.r('min_hit_count_name<-chosendataset[ind_min_hit_count,2]')
    #minhitstation=robjects.r('x[which.min(x[,2]),1]')

    minreadcount=robjects.r('min_read_count<-min(chosendataset$alltranscripts)')
    whichminhitcount=robjects.r('ind_min_read_count<-which.min(chosendataset$alltranscripts)')
    minreadstation=robjects.r('min_read_count_name<-chosendataset[ind_min_read_count,2]')

    # fractions
    #fraccount=robjects.r('chosendataset$fraction<-chosendataset$n_hit_station/chosendataset$alltranscripts')
    #minfrac=robjects.r('print(chosendataset)')
    minfrac=robjects.r('min_frac<-min(chosendataset$fraction)')
    whichminfraccount=robjects.r('ind_min_frac<-which.min(chosendataset$fraction)')
    minfracstation=robjects.r('min_frac_name<-chosendataset[ind_min_frac,2]')

    text_size=4

    make_text_plot+='+geom_text(aes(x = 5, y = 95, label = "Investigation of Insertion Sequences in the transcriptome"), size=6,data = plotdf, colour = I("black"),hjust=0)'

    # Is counts
    make_text_plot+='+geom_text(aes(x = 5, y = 85, label = "Highest surface is count='+str(maxhitcount[0])+'"), size=6,data = plotdf, colour = I("black"),hjust=0)'
    make_text_plot+='+geom_text(aes(x = 4, y = 75, label = max_hit_count_name), size='+str(text_size)+',data = plotdf, colour = I("black"),hjust=0)'

    # Transcript counts
    make_text_plot+='+geom_text(aes(x = 5, y = 70, label = "Highest '+depth_title+' transcript count='+str(maxreadcount[0])+'"), size='+str(text_size)+',data = plotdf, colour = I("black"),hjust=0)'
    make_text_plot+='+geom_text(aes(x = 4, y = 65, label = max_read_count_name), size='+str(text_size)+',data = plotdf, colour = I("black"),hjust=0)'

    # Frac counts
    make_text_plot+='+geom_text(aes(x = 5, y = 60, label = "Highest '+depth_title+' fraction='+str(maxfrac[0])+'"), size='+str(text_size)+',data = plotdf, colour = I("black"),hjust=0)'
    make_text_plot+='+geom_text(aes(x = 4, y = 55, label = max_frac_name), size='+str(text_size)+',data = plotdf, colour = I("black"),hjust=0)'

    # Is counts
    make_text_plot+='+geom_text(aes(x = 5, y = 45, label = "Lowest '+depth_title+' is count='+str(minhitcount[0])+'"), size=6,data = plotdf, colour = I("black"),hjust=0)'
    make_text_plot+='+geom_text(aes(x = 4, y = 40, label = min_hit_count_name), size='+str(text_size)+',data = plotdf, colour = I("black"),hjust=0)'

    # Transcript counts
    make_text_plot+='+geom_text(aes(x = 5, y = 35, label = "Lowest '+depth_title+' transcript count='+str(minreadcount[0])+'"), size='+str(text_size)+',data = plotdf, colour = I("black"),hjust=0)'
    make_text_plot+='+geom_text(aes(x = 4, y = 30, label = min_read_count_name), size='+str(text_size)+',data = plotdf, colour = I("black"),hjust=0)'

    # Frac counts
    make_text_plot+='+geom_text(aes(x = 5, y = 25, label = "Lowest '+depth_title+' fraction='+str(minfrac[0])+'"), size='+str(text_size)+',data = plotdf, colour = I("black"),hjust=0)'
    make_text_plot+='+geom_text(aes(x = 4, y = 20, label = min_frac_name), size='+str(text_size)+',data = plotdf, colour = I("black"),hjust=0)'
    make_text_plot+='+scale_x_continuous(limits=c(1, 100))'
    make_text_plot+='+scale_y_continuous(limits=c(1, 100))'
    make_text_plot+='+theme('
    make_text_plot+='legend.position="none"'
    make_text_plot+=',axis.title.x=element_blank()'
    make_text_plot+=',axis.title.y=element_blank()'
    make_text_plot+=',axis.text.x=element_blank()'
    make_text_plot+=',axis.text.y=element_blank()'
    make_text_plot+=',axis.ticks = element_blank()'
    make_text_plot+=',panel.grid.minor=element_blank()'
    make_text_plot+=',panel.grid.major=element_blank()'
    make_text_plot+=',panel.border=element_blank()'
    make_text_plot+=',panel.background = element_rect(fill="white", colour="black")'
    
    make_text_plot+=')'
    robjects.r(make_text_plot)


     
    robjects.r('cols <- colorRampPalette(brewer.pal(9, "Set1"))')
    #station_dataset=robjects.r('molten_dataset <- melt(chosendataset)')
    station_dataset=robjects.r('molten_dataset <- melt(chosendataset,id.vars = "mystation", measure.vars = c("alltranscripts", "n_hit_station"))')
    

    make_plot='oio2<-ggplot(data=molten_dataset, aes(x=mystation, y=value,fill=variable,order = -as.numeric(variable)))'     
    make_plot+=' + geom_bar(stat="identity",color="black")'
    make_plot+=' + scale_fill_manual(values=cols(2))'
    #make_plot+=' + geom_bar(stat="identity",position=position_dodge())'
    make_plot+=' + theme(axis.text.x=element_text(angle=90))'



#    robjects.r('print(cols)')
#
#    allISfiltersgraph=robjects.r('allplot<-ggplot(data=testtt, aes(x=testtt$mystation, y=testtt$value,fill=testtt$variable))') 
#    allISfiltersgraph=robjects.r('allplot<-allplot+geom_bar(color="black",stat="identity",position="fill",size=0.2)')
#
#    allISfiltersgraph=robjects.r('allplot<-allplot+theme(axis.text.x=element_text(angle=90))')
#    #allISfiltersgraph=robjects.r('allplot<-allplot+theme(axis.text.x=element_text(angle=90),axis.text.x = element_blank())')
#





    robjects.r(make_plot)
#teest=robjects.r('p<-p+coord_fixed(ratio=400)')
    robjects.r('oio2 <- oio2+ggtitle("'+graphtitle+'")')
    robjects.r('oio2 <- oio2+theme(plot.title = element_text(lineheight=.8, face="bold"))')
    #robjects.r.ggsave("/Users/security/science/ttest.pdf")
    robjects.r('pdf("'+outname+'", width='+str(w)+',height='+str(h)+')')
    robjects.r('grid.arrange(oio2,textplot,ncol=2)')
    #im4=robjects.r('print(oio2)')
    robjects.r('dev.off()')

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
    robjects.r('oio <- ggmap(map) + geom_point(geom_point(position = "jitter"),data=chosendataset,color="red",aes(x=Lon, y = Lat,size=alltranscripts1))')
    robjects.r('oio <- oio + geom_text(geom_point(position = "jitter"),data=chosendataset,size=3,color="red",aes(x=Lon,y=Lat,label=mystation),hjust=0.5, vjust=-0.7)')
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

    graphtitle="IS fraction of transcripts"
    station_dataset=robjects.r('chosendataset_surf_shallow <- read.delim("'+incsv_surf+'", sep="," , header=TRUE)')
    allISfiltersgraph=robjects.r('bargall_shallow_frac<-ggplot(data=chosendataset_surf_shallow, aes(x=mystation, y='+"n_hit_station/alltranscripts"+')) ')
    lol2='bargall_shallow_frac<-bargall_shallow_frac+'+get_bar_settings(graphtitle,"IS frac. of total")
    robjects.r(lol2)

    #graphtitle="Number of IS transcripts"
    #station_dataset=robjects.r('chosendataset_surf_shallow <- read.delim("'+incsv_surf+'", sep="," , header=TRUE)')
    #allISfiltersgraph=robjects.r('bargall_shallow_is<-ggplot(data=chosendataset_surf_shallow, aes(x=mystation, y='+"n_hit_station"+',fill="deeppink4")) ')
    #lol3='bargall_shallow_is<-bargall_shallow_is+'+get_bar_settings(graphtitle,"Total reads")
    #robjects.r(lol3)


    #graphtitle="Total number of transcripts"
    #plotline='bargall_shallow<-ggplot(data=chosendataset_surf, aes(x=mystation, y='+"alltranscripts"+'),fill="blue",colour="blue")'
    #robjects.r(plotline)
    #lol='bargall_shallow<-bargall_shallow+'+get_bar_settings(graphtitle,"Number of reads")
    #robjects.r(lol)

    graphtitle="Environmental variables"
    plotline='environment_shallow<-ggplot(data=chosendataset_surf)'
    robjects.r(plotline)
    #lol2='environment_shallow<-environment_shallow'
    lol2='environment_shallow<-environment_shallow+'+get_line_settings("chosendataset_surf",graphtitle,"normalized")
    robjects.r(lol2)


# Three bar plots for medium 
    station_dataset=robjects.r('chosendataset_mid <- read.delim("'+incsv_mid+'", sep="," , header=TRUE)')

    graphtitle="IS fraction of transcripts"
    robjects.r('bargall_mid_frac<-ggplot(data=chosendataset_mid, aes(x=mystation, y='+"n_hit_station/alltranscripts"+',fill="'+mycolor+'"))')
    lol='bargall_mid_frac<-bargall_mid_frac+'+get_bar_settings(graphtitle,"IS frac. of total")
    robjects.r(lol)

    #graphtitle="Number of IS transcripts"
    #robjects.r('bargall_mid_hits<-ggplot(data=chosendataset_mid, aes(x=mystation, y='+"n_hit_station"+',fill="green"))')
    #lol='bargall_mid_hits<-bargall_mid_hits+'+get_bar_settings(graphtitle,"Reads with ISs")
    #robjects.r(lol)

    #graphtitle="Total numbers of transcripts"
    #robjects.r('bargall_mid_total_transcripts<-ggplot(data=chosendataset_mid, aes(x=mystation, y='+"alltranscripts"+',fill="'+mycolor+'"))')
    #lol='bargall_mid_total_transcripts<-bargall_mid_total_transcripts+'+get_bar_settings(graphtitle,"Number of reads")

    graphtitle="Environmental variables"
    plotline='environment_mid<-ggplot(data=chosendataset_mid)'
    robjects.r(plotline)
    #lol2='environment_shallow<-environment_shallow'
    lol2='environment_mid<-environment_mid+'+get_line_settings("chosendataset_mid",graphtitle,"normalized")
    robjects.r(lol2)

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

    sizekind="n_hit_station/alltranscripts"
    oncolor="#70AF7C"
    fill="green"
    shape=16
    pointline='surf_points<-mymap +geom_point(data=chosendataset_surf,shape='+str(shape)+',color="'+oncolor+'",aes(x=Lon, y = Lat'
    pointline+=',size='+sizekind
    pointline+=')'
    pointline+=''
    pointline+=')+scale_size_area(max_size=7)'
    robjects.r(pointline)

    #robjects.r('surf_points<-mymap +geom_point(data=chosendataset_surf,color="red",aes(geom_point(position = "jitter"),x=Lon, y = Lat,size=n_hit_station/alltranscripts))')
#scale_shape(solid = FALSE)

#    robjects.r('surf_points<-mymap +geom_point(data=chosendataset_surf,shape=2,color="blue",position="jitter",aes(x=Lon, y = Lat,size=alltranscripts))')
    #robjects.r('surf_points<-mymap +geom_point(geom_point(data=chosendataset_surf,color="blue",aes(position = "jitter"),x=Lon, y = Lat,size=alltranscripts))')
    robjects.r('surf_points <- surf_points+ geom_text(data=chosendataset_surf,size=4,color="red",aes(x=Lon,y=Lat,label=mystation),hjust=0.5, vjust=-0.7)')
    #robjects.r('surf_points<-surf_points +ggtitle("'+graphtitle+'")')


    robjects.r('mid_points<-mymap +geom_point(data=chosendataset_mid,color="#70AF7C",aes(x=Lon, y = Lat,size=alltranscripts))')

    robjects.r('mid_points <- mid_points+ geom_text(data=chosendataset_mid,size=4,color="red",aes(x=Lon,y=Lat,label=mystation),hjust=1, vjust=-0.7)')
    #robjects.r('mid_points <- mid_points+ggtitle("'+graphtitle+'")')


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
    grid_arrangement+='arrangeGrob(bargall_mid_frac,bargall_mid_total_transcripts,bargall_mid_hits,ncol=1,heights=c(2,1,12,1,1)),'
    #grid_arrangement+='deep_points,'
    #grid_arrangement+='arrangeGrob(bargall_deep,bargall_deep,ncol=1),' 
    grid_arrangement+='ncol=2,'
    grid_arrangement+='main=textGrob("Insertion sequences in the Baltic Sea 2009", gp = gpar(fontsize=18, fontface="bold.italic", fontsize=18))'
    grid_arrangement+=',widths=1:1)'

    hh='grid.arrange( arrangeGrob( surf_points,arrangeGrob(bargall_shallow_frac,environment_shallow,ncol=1,heights=c(1,1)) ,ncol=2,main="Shallow (0.5m)"),arrangeGrob(mid_points,arrangeGrob(bargall_mid_frac,environment_mid,ncol=1,heights=c(1,1)) ,ncol=2,main="Deep") ,ncol=1,main=textGrob("\nInsertion sequences in the Baltic Sea 2009\n", gp = gpar(fontsize=18, fontface="bold", fontsize=18)),widths=1:1)'
    #hh='grid.arrange( arrangeGrob( surf_points,arrangeGrob(bargall_shallow_frac,bargall_shallow_is,bargall_shallow,ncol=1,heights=c(2,1,1)) ,ncol=2,main="Shallow"),arrangeGrob(mid_points,arrangeGrob(bargall_mid_frac,bargall_mid_hits,environment_shallow,ncol=1,heights=c(2,1,1)) ,ncol=2,main="Deep") ,ncol=1,main=textGrob("\nInsertion sequences in the Baltic Sea 2009\n", gp = gpar(fontsize=18, fontface="bold", fontsize=18)),widths=1:1)'
    #hh='grid.arrange( arrangeGrob( surf_points,arrangeGrob(bargall_shallow_frac,bargall_shallow_is,bargall_shallow,ncol=1,heights=c(2,1,1)) ,ncol=2,main="Shallow"),arrangeGrob(mid_points,arrangeGrob(bargall_mid_frac,bargall_mid_hits,bargall_mid_total_transcripts,ncol=1,heights=c(2,1,1)) ,ncol=2,main="Deep") ,ncol=1,main=textGrob("\nInsertion sequences in the Baltic Sea 2009\n", gp = gpar(fontsize=18, fontface="bold", fontsize=18)),widths=1:1)'




    #print grid_arrangement
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
#csvh=open("/Users/security/Desktop/headers2.csv","r")
#csread=csvh.read()
#css=csread.split("\n")
#
#
#heads=[]
#for u in css:
#    if u=="":
#        continue
#    uu=u.split(",")
#    heads.append(uu[1])
#
#print "heads"
#print heads
#
#
#stations=[]
#for i in range(33,87):
#   # print heads[i]
#    stations.append(heads[i])
#
#
#
#
#dethg=""
#for d in heads:
#    if "3p0" in d:
#        dethg+="SUM("+d+")+"
#        #dethg+=d+","
#
#
#dethg=dethg[:-1]
##columns="SUM("+dethg+")"
##columns=dethg
#print "dethg"
#print dethg
#
#
#
#

robjects.r('library(reshape2)')
robjects.r('library(RColorBrewer)')
robjects.r('library(stats)')

summarycsv=""


#selstring="SELECT "+"GS670_0p8"+" FROM allist WHERE group_='"+groupthing+"'"
#selstring="SELECT "+columns+" FROM allist WHERE group_='"+groupthing+"'"

# Sample one column
############
#columns="SUM("+"GS670_0p8"+")"
#selstring="SELECT "+columns+" FROM allist"
#with con:
#        cur.execute(selstring)
#reply=cur.fetchall()
#onesampl=reply[0][0]
#summarycsv+="Sample: sum of all transcripts at GS670_0p8:,"+str(onesampl)+"\n"


# Sample one column with hits
############
#columns="SUM("+"GS670_0p8"+")"
#selstring="SELECT "+columns+" FROM allist WHERE hit_abb!='nod'"
##selstring="SELECT "+columns+" FROM allist WHERE isname!=''"
#with con:
#        cur.execute(selstring)
#reply=cur.fetchall()
#onesampl=reply[0][0]
#summarycsv+="Sample: sum of all hit transcripts at GS670_0p8:,"+str(onesampl)+"\n"


# compile a list with stations
###########
stationcsv="name,mystation,filtersize,alltranscripts,n_hit_station,"+metahead.replace("\t",",")+"\n"
#if not os.path.isfile(args.outcsv+"stats.csv"):
#    for mystation in stations:
#        print "getting hit counts for "+mystation
#        stat_loc=mystation.split("_")[0]
#        stat_filt=mystation.split("_")[1]
#
#        columns="SUM("+mystation+")"
#        selstring="SELECT "+columns+" FROM allist"
#        with con:
#                cur.execute(selstring)
#        reply=cur.fetchall()
#        onesampl=reply[0][0]
#
#        columns="SUM("+mystation+")"
#        selstring="SELECT "+columns+" FROM allist WHERE hit_abb!='nod'"
#        with con:
#                cur.execute(selstring)
#        reply=cur.fetchall()
#        hitcount=reply[0][0]
#
#        stationcsv+=mystation+","+stat_loc+","+stat_filt+","+str(onesampl)+","+str(hitcount)+","+getfromcc(mystation.split("_")[0])+"\n"
##
#    open(args.outcsv,"w").write(summarycsv)
#    fhh=open(args.outcsv+"stats.csv","w")
#    fhh.write(stationcsv)
#    fhh.close()

# Main station info csv
###################
station_dataset=robjects.r('statdataset <- read.delim("'+args.outcsv+"stats.csv"+'", sep="," , header=TRUE)')
#robjects.r('attach(statdataset)')


grdevices = importr('grDevices')


# Make stationsdb
#################
stations_db_place=basepicdir+"stationsdb.db"
con4 = sqlite3.connect(stations_db_place)
cur4 = con4.cursor()

dropstring="DROP TABLE IF EXISTS "+"t"
cur4.execute(dropstring)
cur4.execute('CREATE TABLE t (name,mystation,filtersize,n_hit_station,alltranscripts,Basin,Size,ChlMax,Lon,Lat,"Water","Water (Ft)",Thermocline,ChlorophyllMaxDepth,O,Sampledepth,Temp,Sal,pH,Chla,NO3,NO2,NH4,Urea,PO4,Si,DOP,NP,P,N,OXY,sig_a_alltranscripts,sig_a_n_hit_station,pro_c_alltranscripts,pro_c_n_hit_station,cpc_g_alltranscripts,cpc_g_n_hit_station);')

#"siga_alltranscripts","sig_a_n_hit_station","pro_c_alltranscripts","pro_c_n_hit_station"
#'siga_alltranscripts','sig_a_n_hit_station','pro_c_alltranscripts','pro_c_n_hit_station'
#SUM(siga_alltranscripts),SUM(sig_a_n_hit_station),SUM(pro_c_alltranscripts),SUM(pro_c_n_hit_station)
#i['siga_alltranscripts'],i[sig_a_n_hit_station'],i[pro_c_alltranscripts'],i[pro_c_n_hit_station']


with open(args.outcsv+"stats.csv",'r') as fin: # `with` statement available in 2.5+
    # csv.DictReader uses first line in file for column headings by default
    dr = csv.DictReader(fin) # comma is default delimiter
    to_db = [(i['name'],i['mystation'],i['filtersize'], i['n_hit_station'],i['alltranscripts'],i['Basin'],i['Size'],i['ChlMax'],i['Lon'], i['Lat'],i['Water'],i['Water (Ft)'],i['Thermocline'],i['ChlorophyllMaxDepth'],i['O'],i['Sampledepth'],i['Temp'],i['Sal'],i['pH'],i['Chla'],i['NO3'],i['NO2'],i['NH4'],i['Urea'],i['PO4'],i['Si'],i['DOP'],i['NP'],i['P'],i['N'],i['OXY'],i['sig_a_alltranscripts'],i['sig_a_n_hit_station'],i['pro_c_alltranscripts'],i['pro_c_n_hit_station'],i['cpc_g_alltranscripts'],i['cpc_g_n_hit_station']) for i in dr]



cur4.executemany('INSERT INTO t (name,mystation,filtersize, n_hit_station,alltranscripts,Basin,Size,ChlMax,Lon, Lat,"Water","Water (Ft)",Thermocline,ChlorophyllMaxDepth,O,Sampledepth, Temp,Sal,pH,Chla,NO3,NO2,NH4,Urea,PO4,Si,DOP,NP,P,N,OXY,sig_a_alltranscripts,sig_a_n_hit_station,pro_c_alltranscripts,pro_c_n_hit_station,cpc_g_alltranscripts,cpc_g_n_hit_station) VALUES (?,?,?, ?,?,?, ?,?,?, ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);', to_db)
con4.commit()


# Made stations db
############

#make csv for all depths, added filtersizes
have_lats_filtadd_csv=basepicdir+"have_lats_filtersadd.csv"
data1 = cur4.execute('SELECT name,mystation, sum(n_hit_station),sum(alltranscripts), ChlMax,Lon, Lat,"Water","Water (Ft)",Thermocline,ChlorophyllMaxDepth, O,Sampledepth ,Temp, Sal,pH,Chla,NO3,NO2,NH4,Urea,PO4,Si,DOP,NP,P,N,OXY,SUM(sig_a_alltranscripts),SUM(sig_a_n_hit_station),SUM(pro_c_alltranscripts),SUM(pro_c_n_hit_station),SUM(cpc_g_alltranscripts),SUM(cpc_g_n_hit_station) FROM t WHERE Lat!="nod2" GROUP BY mystation')
#data = cur4.execute('SELECT name,mystation, sum(n_hit_station),sum(alltranscripts), ChlMax,Lon, Lat,"Water(m)" AS Water,"Water (Ft)",Thermocline,ChlorophyllMaxDepth, O,Sampledepth Temp, Sal,pH,Chla,NO3,NO2,NH4,Urea,PO4,Si,DOP,NP,P,N,OXY,SUM(siga_alltranscripts),SUM(sig_a_n_hit_station),SUM(pro_c_alltranscripts),SUM(pro_c_n_hit_station) FROM t WHERE Lat!="nod2" GROUP BY mystation')
with open(have_lats_filtadd_csv, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['name','mystation','n_hit_station','alltranscripts','ChlMax','Lon', 'Lat','Water','Water (Ft)','Thermocline','ChlorophyllMaxDepth','O','Sampledepth','Temp','Sal','pH','Chla','NO3','NO2','NH4','Urea','PO4','Si','DOP','NP','P','N','OXY','sig_a_alltranscripts','sig_a_n_hit_station','pro_c_alltranscripts','pro_c_n_hit_station','cpc_g_alltranscripts','cpc_g_n_hit_station'])
    #writer.writerow(['name','mystation','n_hit_station','alltranscripts','ChlMax','Lon', 'Lat','Water','Water (Ft)','Thermocline','ChlorophyllMaxDepth','O','Sampledepth','Temp','Sal','pH','Chla','NO3','NO2','NH4','Urea','PO4','Si','DOP','NP','P','N','OXY','siga_alltranscripts','sig_a_n_hit_station','pro_c_alltranscripts','pro_c_n_hit_station'])
    writer.writerows(data1)



# Make bar graphs with filtersizes added
graphtitle="Total number of reads"
#bar_allsize(have_lats_filtadd_csv,basepicdir+"allbars_alltranscripts.pdf",10,10,graphtitle,"alltranscripts")
#graphtitle="Total number of rpod"
#bar_allsize(have_lats_filtadd_csv,basepicdir+"allbars_rpods.pdf",10,10,graphtitle,"n_hit_station")
#graphtitle="n_hit_station_alltranscripts"
#bar_allsize(have_lats_filtadd_csv,basepicdir+"allbars_rpods_div_alltranscripts.pdf",10,10,graphtitle,"n_hit_station/alltranscripts")
#
## Make more bar graphs with added filtersizes
#bars(have_lats_filtadd_csv,basepicdir+"bars_ouput.pdf",10,10,graphtitle,"alltranscripts")
makepairs(have_lats_filtadd_csv,basepicdir+"environmental_pairs.png",20,20,graphtitle,"alltranscripts")

make_model(have_lats_filtadd_csv,basepicdir+"model.pdf",5,20,"graphtitle","model")
#
#
#
##make csv for all depths, all filtersizes
#have_lats_csv=basepicdir+"have_lats.csv"
#data = cur4.execute("SELECT name,mystation, filtersize,n_hit_station,alltranscripts, Lon, Lat, O,Sampledepth, Temp,Sal,pH,Chla,NO3,NO2,NH4,Urea,PO4,Si,DOP,NP,P,N,OXY FROM t WHERE Lat!='nod2'")
#with open(have_lats_csv, 'w') as f:
#    writer = csv.writer(f)
#    writer.writerow(['name','mystation','n_hit_station','alltranscripts','Lon', 'Lat','O','Sampledepth','Temp','Sal','pH','Chla','NO3','NO2','NH4','Urea','PO4','Si','DOP','NP','P','N','OXY'])
#    writer.writerows(data)
#
##graphtitle="Total number of reads"
##bar_allsize(have_lats_csv,basepicdir+"allbars_alltranscripts.pdf",10,10,graphtitle,"alltranscripts")
##graphtitle="Total number of rpod"
##bar_allsize(have_lats_csv,basepicdir+"allbars_rpods.pdf",10,10,graphtitle,"n_hit_station")
##graphtitle="n_hit_station_alltranscripts"
##bar_allsize(have_lats_csv,basepicdir+"allbars_rpods_div_alltranscripts.pdf",10,10,graphtitle,"n_hit_station/alltranscripts")
##
## Collect data with filter sizes
##data = cur4.execute("SELECT name,mystation, alltranscripts, Lon, Lat, Sampledepth FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)<=0.30487806")
#
#surfLats_csv=basepicdir+"surf_lats.csv"
## Make database with latitudes with sampledepth < 0.3
########################
##data = cur4.execute("SELECT * FROM t")
##data = cur4.execute("SELECT * FROM t WHERE Lat!='nod2'")
#data = cur4.execute("SELECT name,mystation, SUM(n_hit_station),SUM(alltranscripts), Lon, Lat, O,Sampledepth ,Temp, Sal,pH,Chla,NO3,NO2,NH4,Urea,PO4,Si,DOP,NP,P,N,OXY FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)<=0.30487806 GROUP BY mystation")
#with open(surfLats_csv, 'w') as f:
#    writer = csv.writer(f)
#    writer.writerow(['name','mystation','n_hit_station','alltranscripts','Lon', 'Lat','O','Sampledepth','Temp','Sal','pH','Chla','NO3','NO2','NH4','Urea','PO4','Si','DOP','NP','P','N','OXY'])
#    writer.writerows(data)
#
#    
#bars(surfLats_csv,basepicdir+"surface_bars_ouput.pdf",10,10,graphtitle,"alltranscripts")
#
#surfLats_csv2=basepicdir+"surf_lats2.csv"
## Make database with latitudes with sampledepth < 0.3
########################
##data = cur4.execute("SELECT * FROM t")
##data = cur4.execute("SELECT * FROM t WHERE Lat!='nod2'")
#data = cur4.execute("SELECT name,mystation, SUM(n_hit_station),SUM(alltranscripts), Lon, Lat, O,Sampledepth, Temp,Sal ,pH,Chla,NO3,NO2,NH4,Urea,PO4,Si,DOP,NP,P,N,OXY FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)<=0.30487806 GROUP BY mystation")
#with open(surfLats_csv2, 'w') as f:
#    writer = csv.writer(f)
#    writer.writerow(['name','mystation','n_hit_station','alltranscripts','Lon', 'Lat','O','Sampledepth','Temp','Sal','pH','Chla','NO3','NO2','NH4','Urea','PO4','Si','DOP','NP','P','N','OXY'])
#    writer.writerows(data)
#
#
## write surface maps
##gogmap_google(surfLats_csv,basepicdir+"surface_hybrid_map.pdf",10,10)
##with_rworld(surfLats_csv,basepicdir+"surface_rworld.pdf",10,10)
##gogmap_google('/Users/security/science/output.csv',"/Users/security/science/surfacemap.pdf")
#
#
#midLats_csv=basepicdir+"mid_lats.csv"
## Make database with latitudes with sampledepth < middle
########################
#data = cur4.execute("SELECT name,mystation, sum(n_hit_station),SUM(alltranscripts), Lon, Lat, O,Sampledepth, Temp, Sal ,pH,Chla,NO3,NO2,NH4,Urea,PO4,Si,DOP,NP,P,N,OXY FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)>0.30487806 AND cast(Sampledepth as float)<20 GROUP BY mystation")
#with open(midLats_csv, 'w') as f:
#    writer = csv.writer(f)
#    writer.writerow(['name','mystation','n_hit_station','alltranscripts','Lon', 'Lat','O','Sampledepth','Temp','Sal','pH','Chla','NO3','NO2','NH4','Urea','PO4','Si','DOP','NP','P','N','OXY'])
#    writer.writerows(data)
#
##gogmap_google(midLats_csv,basepicdir+"middle_hybrid_map.pdf",10,10)
##with_rworld(midLats_csv,basepicdir+"middlemap_rworld.pdf",10,10)
#
#
#bars(midLats_csv,basepicdir+"deep_bars_ouput.pdf",10,10,graphtitle,"alltranscripts")
#
##deepLats_csv=basepicdir+"deep_lats.csv"
### Deep lvl
##################
##data = cur4.execute("SELECT name,mystation, SUM(alltranscripts), Lon, Lat, Sampledepth FROM t WHERE Lat!='nod2' AND cast(Sampledepth as float)>20 AND cast(Sampledepth as float)<200 GROUP BY mystation")
##with open(deepLats_csv, 'w') as f:
##    writer = csv.writer(f)
##    writer.writerow(['name','mystation','alltranscripts','Lon', 'Lat','Sampledepth'])
##    writer.writerows(data)
#
#
##gogmap_google(deepLats_csv,basepicdir+"deep_hybrid_map.pdf",10,10)
##with_rworld (deepLats_csv,basepicdir+"deep_rworld_map.pdf",10,10)
#
#
#gogmap3(surfLats_csv,midLats_csv,basepicdir+"deep_map_multi.pdf",10,10)
#
#
#

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




