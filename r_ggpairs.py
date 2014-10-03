import rpy2.robjects as robjects
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.packages import importr
import argparse
#import pandas as pd
#import numpy as np
#import rpy2.robjects.pandas2ri
from rpy2.robjects.vectors import DataFrame
#import math
#import datetime

parser = argparse.ArgumentParser()


parser.add_argument("-in_csv", help="")
parser.add_argument("-out", help="")

args=parser.parse_args()
robjects.r('''
    panel.smooth <- function (x, y) {
        points(x, y)
        abline(lm(y~x), col="red")
        lines(stats::lowess(y~x), col="blue")
    }
''')

#dataf = DataFrame.from_csvfile(args.in_csv, sep = "\t",header=True)

grdevices = importr('grDevices')
gal = importr('GGally')

robjects.r('png("'+args.out+'",width=2000,heigh=2000)')
#robjects.r('postscript("rplot.eps")')
#grdevices.pdf(file="/Users/security/science/pairs_test3.pdf")


#dataf=robjects.r('mydata <- read.delim("/Users/security/science/correctstations4/joined4.csv", sep="," , header=TRUE)')
dataf=robjects.r('mydata <- read.delim("'+args.in_csv+'", sep="," , header=TRUE)')

yet=robjects.r('keeps <- c("ratio","alltranscripts","Chla","Lat","Lon","Sal","Temp","pH","O","N","NO3","NH4")')
#yet=robjects.r('keeps <- c("istranscripts","alltranscripts","Chla","Lat","Lon","Sal","Temp","N")')
isnum=robjects.r('attach(mydata)')
#isnum3=robjects.r('iss<-istranscripts')
#isnum3=robjects.r('alltrans<-alltranscripts')
hh=robjects.r('ne<-mydata[keeps]')
print hh

#isnum3=robjects.r('ratio<-iss/alltranscripts')
#im4=robjects.r('oio<-ggpairs(ne,upper = list(params = c(size = 10)),panel = panel.smooth)')
#im4=robjects.r('oio<-ggpairs(ne,lower=panel.smooth)')
im4=robjects.r('oio<-ggpairs(ne,diag=list(labelSize=4),upper = list(params = c(size = 4)),lower = list(continuous = "smooth", size=0.01,params = c(colour="red",method = "glm")))')
#im4=robjects.r('oio<-ggpairs(ne,diag=list(labelSize=4),upper = list(params = c(size = 4)),lower = list(continuous = "smooth", size=0.01,params = c(colour="red",method = "loess")))')
#im4=robjects.r('oio<-ggpairs(ne,upper = list(params = c(corSize = 1),axisLabels="show"),lower = list(continuous = "smooth", size=0.01,params = c(colour="red",method = "loess")))')
im4=robjects.r('print(oio)')
#gp = ggplot2.ggplot(dataf)
#gp = ggplot2.ggpairs(dataf)
#robjects.r('ggsave("/Users/security/science/new3.pdf")')

robjects.r('dev.off()')
