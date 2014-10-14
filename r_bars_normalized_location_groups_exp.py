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
#parser.add_argument("-in_csv_all", help="")
parser.add_argument("-out", help="")

args=parser.parse_args()


robjects.r('library(reshape2)')
robjects.r('library(RColorBrewer)')
robjects.r('library(stats)')

dataf=robjects.r('mydataset <- read.delim("'+args.in_csv+'", sep="," , header=TRUE)')
robjects.r('attach(mydataset)')

#robjects.r('x <- matrix(runif(100), ncol = 5)')
robjects.r('print(mydataset)')
#robjects.r('group <- sample(1:8, 20, TRUE)')
robjects.r('group <- c(FALSE,TRUE,FALSE)')
robjects.r('(xsum <- rowsum(mydataset, na.rm=TRUE,group))')
robjects.r('print(xsum)')

#o3p0=robjects.r('only_depth<-mydataset[grep("'+depth+'", mydataset$mystation), ]')



#for depth in ["0p1","0p8","3p0"]:
#
##only_w_filters=robjects.r('onlyfilters<-subset(mydataset,mystation %in% c("0p1","0p8","3p0"))')
#    #o0p1=robjects.r('only0p1<-mydataset[grep("GS667_3p0", mydataset$GS667_3p0), ]')
#    #o0p8=robjects.r('only0p8<-mydataset[grep("0p8", mydataset$mystation), ]')
#    #o3p0=robjects.r('only3p0<-mydataset[grep("3p0", mydataset$mystation), ]')
#
#
#    o3p0=robjects.r('only_depth<-mydataset[grep("'+depth+'", mydataset$mystation), ]')
#
#    #olls=robjects.r('allfiltersizes<-only0p1$Cyanobacteria+only0p8$Cyanobacteria+only3p0$Cyanobacteria')
#
#    #stats=robjects.r('statss<-only0p1$mystation')
#
#    #robjects.r('newframe<-data.frame(statss,allfiltersizes)')
#
#    #robjects.r('y_name <- "transcripts"')
#    #robjects.r('x_name <- "stations"')
#
#    #robjects.r('names(newframe)<-c(x_name,y_name)')
#
#
#    robjects.r('library(reshape2)')
#    robjects.r('library(RColorBrewer)')
#
## remove landsort
#    robjects.r('new_d <- only_depth[! grepl("GS678",only_depth$mystation),]')
#
## only landsort
#    #robjects.r('onlyls <- only_depth[grepl("GS678",only_depth$mystation),]')
#    
#    #robjects.r('print(onlyls)')
#
#
#    robjects.r('testtt<-melt(new_d)')
#
#    #robjects.r('print(testtt)')
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
#    allISfiltersgraph=robjects.r('allplot<-allplot+labs(title="'+depth+'")')
#
#
##allISfiltersgraph=robjects.r('allplot<-allplot+guides(fill = guide_legend(reverse=TRUE,title = "LEFT", title.position = "left"))')
##allISfiltersgraph=robjects.r('allplot<-allplot+guide_legend(title=NULL)')
#
##robjects.r('allplot<-ggplot(aes(x=statss, y=total)) ')
#
#
#
#
#    robjects.r.ggsave(args.out.replace(".pdf","")+depth+".pdf")
#
#
#    """
#    0.028 
#    archeplastida
#    677
#
#    """


