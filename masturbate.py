import cPickle as pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-place_contig_dic", help="outcsv")
parser.add_argument("-outcsv", help="outcsv")
parser.add_argument("-hitlist", help="outcsv")



args=parser.parse_args()

#work_csv="/Users/security/science/trydic2.csv"
work_csv=args.place_contig_dic
wh=open(work_csv,"r")


def getWhatIS(transcname):
    return "1"

# placelist corrected
placelist=["GS675_3p0", "GS675_0p8", "GS675_0p1", "GS676_3p0", "GS676_0p8","GS676_0p1","GS677_3p0", "GS677_0p8","GS677_0p1","GS678_3p0","GS667_3p0", "GS678_0p8","GS678_0p1","GS679_0p1","GS679_0p8","GS679_3p0","GS667_0p8","GS680_3p0","GS680_0p8","GS680_0p1","GS683_3p0","GS683_0p8","GS683_0p1","GS684_3p0","GS684_0p8","GS684_0p1","GS667_0p1","GS694_3p0","GS694_0p8","GS694_0p1","GS695_3p0","GS695_0p8","GS695_0p1","GS669_3p0","GS669_0p8","GS669_0p1","GS670_3p0","GS670_0p8","GS670_0p1","GS845_ls3","GS846_ls3","GS848_ls3","GS850_ls4","GS852_ls4","GS853_ls4","GS855_ls4","GS856_ls5","GS857_ls5","GS859_ls5","GS860_ls5","LD30M_ls2","LD3200M_ls1","LD35M_ls2","LD390M_ls1"]


di=pickle.load(wh)
#print di.keys()


pifh=open(args.hitlist,"r")
#pifh=open("/Users/security/science/iscontiglibrary.dic","r")
islist=pickle.load(pifh)

csvout=""

for numkey in range(33,87):

    # place with list of transcripts, and their number of hits 
    placedic=di[numkey]

    nameidx=numkey-33
    placename=placelist[nameidx]


    readscount=0
    iscount=0


    # for each transcript in place
    for transc in placedic.keys():
        # get number of transcripts for this place
        numtransc=int(placedic[transc])

        if transc in islist:
            iscount+=numtransc
        whatIS=getWhatIS(transc)
        readscount+=numtransc

    dsp=placename.split("_")

    isfrac=int(iscount)/float(readscount)

    csvout+=str(numkey)+","+placename+","+dsp[0]+","+dsp[1]+","+str(readscount)+","+str(iscount)+","+str(isfrac)+"\n"
        


mid=open(args.outcsv,"w")
mid.write(csvout)
