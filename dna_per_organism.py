import argparse
import cPickle as pickle
import xlrd

parser = argparse.ArgumentParser()
parser.add_argument("-is_contigs_count_csv", help="dirwithall")
parser.add_argument("-group_out_dir", help="dirwithall")
args=parser.parse_args()


#data_csv=args.is_contigs_count_csv
data_csv="/Users/security/science/5oct2-2.csv"

#detailedfh=open(long_csv,"r")
#detailedread=detailedfh.read()
#detailedrows=detailedread.split("\n")
#headline=detailedrows[0]
#detailedrows=detailedrows[1:]



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

csvs={}

grops=['Actinobacteria', 'Bacteroidetes/Chlorobi group', 'Other Bacteria', 'Alphaproteobacteria', 'Archaeplastida', 'Viruses', 'Cyanobacteria', 'Betaproteobacteria', 'Deltaproteobacteria', 'Gammaproteobacteria', 'Archaea', 'Hacrobia', 'Metazoa', 'Stramenopiles', 'Epsilonproteobacteria', 'Alveolata', 'Fungi', 'Other Eukaryota']

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

domains=[]
groups=[]


for row in datarows:
    print "datarow"
    print row
    if row=="":
        continue
    tabrow=row.split(",")
    if tabrow==['','']:
        continue
    if count%10000==0:
        print count
    #if count>10000:
    #    continue
    count+=1
    cellinf=""

    dom=tabrow[21]
    grp=tabrow[22]
    if grp=="":
        grp="unknown"
    #outrow+=ma_")+","
        #print str(tabrow[87])
    if len(tabrow)<21:
        print tabrow
        #raw_input()
    if grp in csvs.keys():
        print grp
        csvs[grp]+=row+"\n"
    else:
        csvs[grp]=""
    if not dom in domains:
        domains.append(dom)
    if not grp in groups:
        groups.append(grp)
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

for typee in csvs.keys():
    outfh=open("/Users/security/science/"+args.group_out_dir+"/"+typee.replace(" ","_").replace("/","")+".csv","w")
    outfh.write(csvs[typee])
    outfh.close()
print "done"

print groups
print domains
print len(groups)
print len(domains)








