
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-incsv", help="dirwithall")
parser.add_argument("-outcsv", help="dirwithall")
args=parser.parse_args()

#infh=open(args.incsv,"r")
infh=open(args.incsv,"r")
#infh=open("/Users/security/science/correctstations6/joined.csv","r")
inread=infh.read()
infh.close()
inrows=inread.split("\n")
outrow=inrows[0]+",totalc,totalis,countd\n"
inrows=inrows[1:]

stationlist=[]


for row in inrows:
    if row=="":
        continue
    tbs=row.split(",")
    station=tbs[2]
    depth=tbs[3]
    if not station in stationlist:
        stationlist.append(station)

for station in stationlist:
    totalc=0
    totalis=0
    countd=0
    statrow=""
    for row in inrows:
        if row=="":
            continue
        tbs=row.split(",")
        mystation=tbs[2]


        if mystation==station:
            statrow=row
            countd+=1
            #print mystation
            iscount=tbs[5].split("0.")[0]
            totalis+=int(iscount)
            allcount=tbs[4].split("0.")[0]
            totalc+=int(allcount)
        
    #print totalc
    #print totalis
    #print countd
    outrow=outrow+statrow+str(totalc)+","+str(totalis)+","+str(countd)+row+"\n"
    #print outrow




infh=open(args.outcsv,"w")
#infh=open("/Users/security/science/correctstations6/joined_all.csv","w")
infh.write(outrow)
#print stationlist



