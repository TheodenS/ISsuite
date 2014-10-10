import argparse
import cPickle as pickle
import xlrd

parser = argparse.ArgumentParser()
parser.add_argument("-is_contigs_count_csv", help="dirwithall")
parser.add_argument("-is_per_group_csv", help="dirwithall")
args=parser.parse_args()




data_csv=args.is_contigs_count_csv
#data_csv="/Users/security/science/5oct2-2.csv"

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

placelist=["GS675_3p0", "GS675_0p8", "GS675_0p1", "GS676_3p0", "GS676_0p8","GS676_0p1","GS677_3p0", "GS677_0p8","GS677_0p1","GS678_3p0","GS667_3p0", "GS678_0p8","GS678_0p1","GS679_0p1","GS679_0p8","GS679_3p0","GS667_0p8","GS680_3p0","GS680_0p8","GS680_0p1","GS683_3p0","GS683_0p8","GS683_0p1","GS684_3p0","GS684_0p8","GS684_0p1","GS667_0p1","GS694_3p0","GS694_0p8","GS694_0p1","GS695_3p0","GS695_0p8","GS695_0p1","GS669_3p0","GS669_0p8","GS669_0p1","GS670_3p0","GS670_0p8","GS670_0p1","GS845_ls3","GS846_ls3","GS848_ls3","GS850_ls4","GS852_ls4","GS853_ls4","GS855_ls4","GS856_ls5","GS857_ls5","GS859_ls5","GS860_ls5","LD30M_ls2","LD3200M_ls1","LD35M_ls2","LD390M_ls1"]

grops=['Actinobacteria', 'Bacteroidetes/Chlorobi group', 'Other Bacteria', 'Alphaproteobacteria', 'Archaeplastida', 'Viruses', 'Cyanobacteria', 'Betaproteobacteria', 'Deltaproteobacteria', 'Gammaproteobacteria', 'Archaea', 'Hacrobia', 'Metazoa', 'Stramenopiles', 'Epsilonproteobacteria', 'Alveolata', 'Fungi', 'Other Eukaryota']

datafh=open(data_csv,"r")
dataread=datafh.read()
datarows=dataread.split("\n")
headline=datarows[0]
datarows=datarows[1:]

# make headline
headline=""
for ii in range(0,53):
    headline+=str(ii)+","
headline=headline[:-1]
outrow=headline+"\n"
count=0




domains=[]
groups=[]

dcc={}
for row in datarows:
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
        grp="unidentified"
    if grp==" ":
        grp="unidentified"

    k=0
    for place in range(33,87):
        placen_name=str(placelist[k])
        place_split=placen_name.split("_")
        plane=place_split[0]
        de=place_split[1]
        k+=1
        #print placen_name
        if not placen_name in dcc.keys():
            dcc[placen_name]={}
        if not grp in dcc[placen_name].keys():
            dcc[placen_name][grp]=[float(str(tabrow[place]))]
        else:
            dcc[placen_name][grp].append(float(str(tabrow[place])))

        #print tabrow[0]
        #print placen_name
        #print grp
        #print tabrow[place]
        


hads=[]
      #outrow+="\n"
header="place,"
place_org_csv=""
for p in dcc.keys():
    print "key "+p+" in dcc.keys"
    placetot=0
    place_org_csv+=p+","
#    header+=p+","
    for uu in dcc[p].keys():
        print "key "+uu+" in dcc.keys"+p
        summ=0
        if not uu in hads:
            hads.append(uu)
        #print dcc[p][uu]
        for num in dcc[p][uu]:
# expecting a list of number at that contig, place"
            summ+=num
            placetot+=num
        place_org_csv+=str(summ)+","
        print summ
    print placetot
    #raw_input()
    place_org_csv+="\n"

for g in hads:
    print hads
    header+=g+","

place_org_csv=header+"\n"+place_org_csv


outfh=open(args.is_per_group_csv,"w")
#outfh=open("/Users/security/science/place_orgn.csv","w")
outfh.write(place_org_csv.encode('ascii','ignore'))
outfh.close()
print "done"
#
#print groups
#print domains
#print len(groups)
#print len(domains)









