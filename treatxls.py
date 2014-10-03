import xlrd
import argparse
import sqlite3
import unicodedata
import pickle

#corrected placelist
placelist=["GS675_3p0", "GS675_0p8", "GS675_0p1", "GS676_3p0", "GS676_0p8","GS676_0p1","GS677_3p0", "GS677_0p8","GS677_0p1","GS678_3p0","GS673_3p0", "GS678_0p8","GS678_0p1","GS679_0p1","GS679_0p8","GS679_3p0","GS667_0p8","GS680_3p0","GS680_0p8","GS680_0p1","GS683_3p0","GS683_0p8","GS683_0p1","GS684_3p0","GS684_0p8","GS684_0p1","GS667_0p1","GS694_3p0","GS694_0p8","GS694_0p1","GS695_3p0","GS695_0p8","GS695_0p1","GS669_3p0","GS669_0p8","GS669_0p1","GS670_3p0","GS670_0p8","GS670_0p1","GS845_ls3","GS846_ls3","GS848_ls3","GS850_ls4","GS852_ls4","GS853_ls4","GS855_ls4","GS856_ls5","GS857_ls5","GS859_ls5","GS860_ls5","LD30M_ls2","LD3200M_ls1","LD35M_ls2","LD390M_ls1"]

#work_csv="/Users/security/science/RNA/transcriptome_Theo/annotation_all.filtered.taxgrps.stats.xlsx"

#restdic="/Users/security/science/trydic2.csv"
#restcsv="/Users/security/science/trycsv2.csv"


parser = argparse.ArgumentParser()
parser.add_argument("-place_contig_dic", help="dictionary that lists number of contigs found at each station")
parser.add_argument("-place_contig_csv", help="csv file that lists number of contigs found at each station")
parser.add_argument("-metadata_csv", help="csv file that lists number of contigs found at each station")
args=parser.parse_args()
#rh=open(restdic,"r")

book = xlrd.open_workbook(args.metadata_csv)
print "the number of worksheets is", book.nsheets
print "Worksheet name(s):", book.sheet_names()
sh = book.sheet_by_index(0)
#nn=deepcopy(sh)
#pickle.dump(nn,rh)
#sh=pickle.load(rh)
print sh.name, sh.nrows, sh.ncols
print "Cell D30 is", sh.cell_value(rowx=0, colx=0)
print "sheet has "+str(sh.nrows)+" rows"
colc=0
rowc=0

k=0
csvl=""

mydic={}

for place in range(33,87):
    place_totalcount=0
    placen_name=str(placelist[k])
    place_split=placen_name.split("_")
    # add full name
    csvl+=str(placelist[k])+","
    # add station
    csvl+=str(place_split[0])+","
    # add depth
    csvl+=str(place_split[1])+","
    k+=1

    if place in mydic.keys():
        pass
    else:
        mydic[place]={}
   
    fracrow=""
    for row in range(1,sh.nrows):
        #print row
            #if rowc>maxrows:
#			continue
        rowc+=1
        #if place<sh.ncols:
        val=sh.cell_value(rowx=row, colx=place)
        name=sh.cell_value(rowx=row, colx=0)
        #else:
        #    print "failed for "+str(place)
        #print val
        #print name
        #mydic [33][contig_20004_2_3]=3 for instance
        mydic[place][name]=val 
        # add number of reads for this contig at this place to total
        place_totalcount+=int(val)
        # add number of reads for this contig at this place to place to csv
        #csvl+=str(int(val))+","
        fracrow+=str(int(val))+","


        csvl+=str(place_totalcount)+","
        csvl+=fracrow+","
        csvl+="\n"
    #print place_totalcount

rh=open(args.place_contig_csv,"w")
dichand=open(args.place_contig_dic,"w")

rh.write(csvl)
pickle.dump(mydic,dichand)


