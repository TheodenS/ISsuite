import argparse
import cPickle as pickle
import xlrd

parser = argparse.ArgumentParser()
#parser.add_argument("-basedir", help="dirwithall")
parser.add_argument("-iscsvlongf", help="dirwithall")
parser.add_argument("-detailed_csvf", help="dirwithall")
parser.add_argument("-contigxls_isinfo_join_csv", help="dirwithall")
args=parser.parse_args()

detailed_csv=args.detailed_csvf
#detailed_csv="/Users/security/science/iscsvdetailed.csv"
#long_csv="/Users/security/science/iscsvlong.csv"
long_csv=args.iscsvlongf
#is_contig_place_dic_file="/Users/security/science/place_contig_dic.pydic"
#is_contig_place_dic_file=ydic"

metadata_csv="/Users/security/science/RNA/transcriptome_Theo/annotation_all.filtered.taxgrps.stats.xlsx"

print "reading detailed csv:"+ detailed_csv
detailedfh=open(detailed_csv,"r")
detailedread=detailedfh.read()
detailedrows=detailedread.split("\n")
headline_det=detailedrows[0]
detailedrows=detailedrows[1:]
print "done reading csv:"+ detailed_csv



def get_detailed_info(contigname):
    hcount=0
    retrow=""
    for detailedrow in detailedrows:
        if detailedrow=="":
            continue
        splitrow=detailedrow.split(",")
        propername=splitrow[2]

        if propername==contigname:
            hcount+=1
            retrow=detailedrow
    return retrow





#dicfh=open(is_contig_place_dic_file,"r")
#is_contig_place_dic=pickle.load(dicfh)

#print is_contig_place_dic.keys()
print "reading xls"
book = xlrd.open_workbook(metadata_csv)
sh = book.sheet_by_index(0)
#nn=deepcopy(sh)
#pickle.dump(nn,rh)
#sh=pickle.load(rh)
print "done xls"

# make headline
headline=""
for col in range(0,sh.ncols):
    val=str(sh.cell_value(rowx=0, colx=col)).replace(",","_comma_")
    #print val
    if val=="GS673_3p0":
        print "HEREEEEEEEEEEEEEEEEEEEEEEE"
        val="GS667_3p0"
    headline+=val+","


outrow=headline+headline_det+"\n"
count=0
for row in range(1,sh.nrows):
    if count%10000==0:
        print count
    count+=1
    for col in range(0,sh.ncols):
        #print "Cell "+str(row)+","+str(col)+" is", sh.cell_value(rowx=row, colx=col)
        outrow+=str(sh.cell_value(rowx=row, colx=col)).replace(",","_comma_")+","
    cellinf=get_detailed_info(sh.cell_value(rowx=row, colx=0))
    outrow+=cellinf
        
        #print row
            #if rowc>maxrows:
#			continue
        #rowc+=1
        #if place<sh.ncols:
        #val=sh.cell_value(rowx=row, colx=col)
        #name=sh.cell_value(rowx=row, colx=0)
    outrow+="\n"




placelist=["GS675_3p0", "GS675_0p8", "GS675_0p1", "GS676_3p0", "GS676_0p8","GS676_0p1","GS677_3p0", "GS677_0p8","GS677_0p1","GS678_3p0","GS667_3p0", "GS678_0p8","GS678_0p1","GS679_0p1","GS679_0p8","GS679_3p0","GS667_0p8","GS680_3p0","GS680_0p8","GS680_0p1","GS683_3p0","GS683_0p8","GS683_0p1","GS684_3p0","GS684_0p8","GS684_0p1","GS667_0p1","GS694_3p0","GS694_0p8","GS694_0p1","GS695_3p0","GS695_0p8","GS695_0p1","GS669_3p0","GS669_0p8","GS669_0p1","GS670_3p0","GS670_0p8","GS670_0p1","GS845_ls3","GS846_ls3","GS848_ls3","GS850_ls4","GS852_ls4","GS853_ls4","GS855_ls4","GS856_ls5","GS857_ls5","GS859_ls5","GS860_ls5","LD30M_ls2","LD3200M_ls1","LD35M_ls2","LD390M_ls1"]


outfh=open("/Users/security/science/5oct2new.csv","w")
outfh.write(outrow)

outfh.close()

outfh=open(args.contigxls_isinfo_join_csv,"w")
outfh.write(outrow)
outfh.close()

print "done"


