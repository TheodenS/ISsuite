#parser.add_argument("-place_contig_dic", help="outcsv")
import argparse

ccresfile="/Users/security/science/ccc2.csv"
metadata="/Users/security/science/RNA/mimebs_metadata/metadata.2009.csv"
outcsv="/Users/security/science/mixedcsv.csv"

parser = argparse.ArgumentParser()
parser.add_argument("-iscounts_contigcsv", help="dirwithall")
parser.add_argument("-outcsv", help="dirwithall")
args=parser.parse_args()

metafh=open(metadata,"r")
metaread=metafh.read()
metalines=metaread.split("\n")
metalines=metalines[1:]

ccresfile=args.iscounts_contigcsv
#print ccresfile

cfh=open(ccresfile,"r")
cread=cfh.read()
#print cread
clines=cread.split("\n")
#clines=clines[1:]

def getfromcc(name):
    badst= "0,no data,no data,no_data,1,0"
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

headl="Number,name,mystation,filtersize,alltranscripts,istranscripts,ratio,Station,Basin,Size,ChlMax,Lat,Lon,Water(m),Water (Ft),Thermocline,ChlorophyllMaxDepth,O,Sampledepth,Temp,Sal,pH,Chla,NO3,NO2,NH4,Urea,PO4,Si,DOP,NP,P,N,OXY\n"
newcsv=headl


for cline in clines:
#for lin in metalines:

    if cline=="":
        continue 
    #lintabs=lin.split("\t")
    clinetab=cline.split(",")
    ccname=clinetab[2]
    #me1name=lintabs[0]
    resu=getfromcc(ccname)

    res=resu.split(",")

    #totalcontigs=int(res[4])
    #iscontigs=int(res[5])
    #rationis=float(iscontigs)/totalcontigs
    #newcsv+=str(totalcontigs)+","+str(iscontigs)+","+str(rationis)+"\n"
    newcsv+=cline+","+resu+"\n"

outfh=open(args.outcsv,"w")
outfh.write(newcsv)
