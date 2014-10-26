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
        propername=splitrow[3]
        if propername==contigname:
            hcount+=1
            retrow=detailedrow
            return retrow
    return "nod,nod,nod,nod,nod,nod,nod,nod,nod,nod,nod,nod,nod,nod,nod,nod"





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
        print "changed faulty station name"
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





outfh=open("/Users/security/science/5oct2new.csv","w")
outfh.write(outrow)

outfh.close()

outfh=open(args.contigxls_isinfo_join_csv,"w")
outfh.write(outrow)
outfh.close()

print "done"


