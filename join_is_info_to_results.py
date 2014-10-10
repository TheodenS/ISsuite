import argparse
import cPickle as pickle
import xlrd

parser = argparse.ArgumentParser()
parser.add_argument("-is_contigs_count_csv", help="dirwithall")
parser.add_argument("-long_csv", help="dirwithall")
parser.add_argument("-out_csv_allcontigs_is_long", help="dirwithall")
args=parser.parse_args()

#detailed_csv="/Users/security/science/iscsvdetailed.csv"
long_csv=args.long_csv
#long_csv="/Users/security/science/iscsvlong.csv"
#is_contig_place_dic_file="/Users/security/science/place_contig_dic.pydic"

#data_csv="/Users/security/science/5oct2.csv"
data_csv=args.is_contigs_count_csv

detailedfh=open(long_csv,"r")
detailedread=detailedfh.read()
detailedrows=detailedread.split("\n")
headline=detailedrows[0]
detailedrows=detailedrows[1:]



def get_detailed_info(contigname):
    #hcount=0
    retrow=""
    for detailedrow in detailedrows:
        if detailedrow=="":
            continue
        splitrow=detailedrow.split(",")
        propername=splitrow[0]

        if propername==contigname:
            #hcount+=1
            retrow=detailedrow
#    if hcount>1:
       # print retrow
    return retrow





#dicfh=open(is_contig_place_dic_file,"r")
#is_contig_place_dic=pickle.load(dicfh)

#print is_contig_place_dic.keys()



datafh=open(data_csv,"r")
dataread=datafh.read()
datarows=dataread.split("\n")
headline=datarows[0]
datarows=datarows[1:]

# make headline
headline=""
for ii in range(0,100):
    headline+=str(ii)+","
headline=headline[:-1]
outrow=headline+"\n"
count=0




for row in datarows:
    tabrow=row.split(",")
    if count%10000==0:
        print count
    #if count>10000:
    #    continue
    count+=1
    cellinf=""
    #outrow+=ma_")+","
    if len(tabrow)>86:
        #print str(tabrow[87])
        cellinf=get_detailed_info(tabrow[87])
        #print cellinf
    outrow+=row+","
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


outfh=open("/Users/security/science/5oct2-2.csv","w")
outfh.write(outrow)
outfh.close()

outfh=open(args.out_csv_allcontigs_is_long,"w")
outfh.write(outrow)
outfh.close()






