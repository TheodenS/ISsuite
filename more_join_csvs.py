import argparse
import cPickle as pickle

parser = argparse.ArgumentParser()
parser.add_argument("-basedir", help="dirwithall")
parser.add_argument("-basefastafile", help="dirwithall")
parser.add_argument("-queryfastafile", help="dirwithall")
parser.add_argument("-maxeval", help="dirwithall")
args=parser.parse_args()

detailed_csv="/Users/security/science/iscsvdetailed.csv"
long_csv="/Users/security/science/iscsvlong.csv"
is_contig_place_dic_file="/Users/security/science/transcless3ke5/place_contig_dic.pydic"


detailedfh=open(detailed_csv,"r")
detailedread=detailedfh.read()
detailedrows=detailedread.split("\n")
headline=detailedrows[0]
detailedrows=detailedrows[1:]


dicfh=open(is_contig_place_dic_file,"r")
is_contig_place_dic=pickle.load(dicfh)

print is_contig_place_dic.keys()




placelist=["GS675_3p0", "GS675_0p8", "GS675_0p1", "GS676_3p0", "GS676_0p8","GS676_0p1","GS677_3p0", "GS677_0p8","GS677_0p1","GS678_3p0","GS667_3p0", "GS678_0p8","GS678_0p1","GS679_0p1","GS679_0p8","GS679_3p0","GS667_0p8","GS680_3p0","GS680_0p8","GS680_0p1","GS683_3p0","GS683_0p8","GS683_0p1","GS684_3p0","GS684_0p8","GS684_0p1","GS667_0p1","GS694_3p0","GS694_0p8","GS694_0p1","GS695_3p0","GS695_0p8","GS695_0p1","GS669_3p0","GS669_0p8","GS669_0p1","GS670_3p0","GS670_0p8","GS670_0p1","GS845_ls3","GS846_ls3","GS848_ls3","GS850_ls4","GS852_ls4","GS853_ls4","GS855_ls4","GS856_ls5","GS857_ls5","GS859_ls5","GS860_ls5","LD30M_ls2","LD3200M_ls1","LD35M_ls2","LD390M_ls1"]

def getlocation(placenum,contig):
    #print is_contig_place_dic[placenum].keys()
    if contig in is_contig_place_dic[placenum]:
        return is_contig_place_dic[placenum][contig.decode('unicode-escape')]
    else:
        return 0



for ey in range(33,87):
    outcsv=placelist[ey]+","
outcsv=+headline+"\n"

for detailedrow in detailedrows:
    if detailedrow=="":
        continue
    numkey=33
    splitrow=detailedrow.split(",")
    propername=splitrow[2]
    numsrow=""
    for numkey in range(33,87):
        ge=getlocation(numkey,propername)
        print ge
        numsrow+=str(ge)+","
    outcsv+=numsrow+","+detailedrow+"\n"
    nameidx=numkey-33
    placename=placelist[nameidx]

outfh=open("/Users/security/science/5oct.csv","w")
outfh.write(outcsv)




