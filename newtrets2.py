import argparse

from Bio import SeqIO
import sqlite3

import os
import pickle

#con = sqlite3.connect("/Users/security/science/short.db")
#cur = con.cursor()

def plusfrac(dic,found_repeat_len,org_len,mylist):

    frac=float(found_repeat_len)/float(org_len)
    for l in mylist:
        maxl=l[0]
        if frac*100>maxl:
	    #print "frac *100 ="+str(frac*100)+" is more than "+str(l[0])
            line=str(l[0])+"to"+str(l[1])
            dic[line]+=1
            break
    return dic

def getISlength(isname):
    for isn in iscsvsplit:
        if isn=="":
            continue
        parts=isn.split(",")
        name=parts[0]
        if name==isname:
            return parts[2]

    
    raw_input("no is found")
    print "for"+str(isname)

def makelengthheadline(frac):
    totline=""
    for div in fraclist:
        totline+=str(div[0])+"to"+str(div[1])+","
    return totline

def make_len_dic():
    dic={}
    counter=0
    for div in fraclist:
        counter+=1
        line=str(div[0])+"to"+str(div[1])
        dic[line]=0
    return dic

#def getKind(name):
#
#	cur.execute('''SELECT Kind FROM repeats WHERE Name=?''', (name,))
#	furrl=cur.fetchone()
#	contigseq=str(furrl[0])
#
#	return contigseq
#
def getSeq(name):
        #print name

	cur2.execute('''SELECT Seq FROM genomeinfo WHERE Myid=?''', (name,))
	#cur2.execute('''SELECT Seq FROM genomeinfo WHERE Name=?''', (name,))
	furrl=cur2.fetchone()
	contigseq=str(furrl[0])

	return contigseq

#def getSeq2(name):
#    for i in fastal:
#        if i.id==name:
#            return len(i.seq)
#    raw_input("prob in getseq2")
#    return 0

def getdir(d):
	retlist=[]
	fileslist=os.listdir(d)
	for u in fileslist:
	    if u[0]==".":
                continue
            else:
		    retlist.append(u)
	return retlist


def getsa(name):
#
        #print name
	#cur2.execute('''SELECT * FROM shortie WHERE Myid=?''', (name,))
	cur2.execute('''SELECT * FROM genomeinfo WHERE Myid=?''', (name,))
	furrl=cur2.fetchone()
	#genomelen=int(furrl[0])
	contigseq=str(furrl[0])

	return contigseq

def makefracstring(contl,fraclist):
    retstring=""
    for div in fraclist:
        line=str(div[0])+"to"+str(div[1])
        retstring+= str(contl[line])+","
    return retstring

fraclist=[[0, 5], [5, 15], [15, 25], [25, 35], [35, 45], [45, 55], [55, 65], [65, 75], [75, 85], [85, 95], [95, 105], [105, 115], [115, 125], [125, 135], [135, 145], [145, 955]]
fraclist.reverse()


parser = argparse.ArgumentParser()
parser.add_argument("-gbfilesdir", help="outcsv")
parser.add_argument("-genomedatabase", help="outcsv")
parser.add_argument("-allgenomesinfo", help="outcsv")
parser.add_argument("-maxeval", help="outcsv")
parser.add_argument("-iscontiglist", help="outcsv")
parser.add_argument("-isdicout", help="outcsv")
parser.add_argument("-iscsv", help="outcsv")
args=parser.parse_args()

totalreps=0
total_length=0


#hugedb="/Users/security/science/correctstations2/databases/"+"genomeinfo.db"
#hugedb="/Users/security/science/short.db"
#hugedb=db"
con2 = sqlite3.connect(args.genomedatabase)
cur2 = con2.cursor()  

csvstring="name,myid,start,end,hitlen,"+makelengthheadline(fraclist)+"\n"

# investigation wide info
inv_n_allf=0
inv_n_nohit=0
inv_n_repeathit=0
inv_n_footprint=0
inv_n_notfound=0
inv_n_treated=0
inv_n_is=0
inv_n_mge=0
inv_n_phage=0
inv_n_putative=0
inv_n_dna=0


tempIS_csv_handle=open(args.iscsv,"r")
iscsvread=tempIS_csv_handle.read()
iscsvsplit=iscsvread.split("\n")
iscsvsplit=iscsvsplit[1:]
#fastal=list(SeqIO.parse(tempIS_csv_handle,"fasta"))


isgenomelist=[]

isdic={}



for genome in getdir(args.gbfilesdir):
        # are these all the subject genomes?
        genome_lendic=make_len_dic()

	# Genome wide info

	#print genome
	#print getsa(genome.replace(".gb",""))
	propername=getsa(genome.replace(".gb",""))
	id_gb_location=args.gbfilesdir
        idfilepath=id_gb_location+genome
        id_file_handle=open(idfilepath,"r")
        gbparse=list(SeqIO.parse(id_file_handle,"genbank"))
	gbrecord=gbparse[0]
	for feature in gbrecord.features:
		if feature.type =="nohit":
			continue	
		if feature.type =="repeathit":
			continue	
		if feature.type =="footprint":
			continue	
		if feature.type =="not found":
			continue	




                #print str(feature)
                #print feature.qualifiers["expected"]

                evalue=float(feature.qualifiers["expected"][0])


                if evalue<float(args.maxeval):

                    # append to list containing genomes with is elements
                    if not propername in isgenomelist:
                        isgenomelist.append(propername)
                    istype=str(feature.type)
                    if istype in isdic.keys():
                        isdic[istype]+=1
                    else:
                        isdic[istype]=1
#
                hit_start=int(feature.location.start.position)
                hit_end=int(feature.location.end.position)

                #if hit_start>hit_end:
                #    raw_input("hit_start is larger than hit_end")

                hitlen=1+hit_end-hit_start

		contig_seq=getSeq(genome.replace(".gb",""))
		#contig_seq=getSeq(propername)



		#
                hitseq=contig_seq[hit_start:hit_end]

                #strand=int(feature.strand)

		total_length+=hitlen


                islent=getISlength(istype)


		fracoforig=str(round(len(hitseq)/float(islent),3))

		genome_lendic=plusfrac(genome_lendic,len(hitseq),islent,fraclist)
                me2= makefracstring(genome_lendic,fraclist)


		csvstring+=propername+","+str(genome)+","+str(hit_start)+","+str(hit_end)+","+str(len(hitseq))+","+me2+"\n"
print total_length

fh=open(args.allgenomesinfo,"w")
fh.write(csvstring)
fh.close()

print isgenomelist
print len(isgenomelist)
#pifh=open("/Users/security/science/iscontiglibrary.dic","w")
pifh=open(args.iscontiglist,"w")
pickle.dump(isgenomelist,pifh)


print isdic
print len(isgenomelist)
pifh2=open(args.isdicout,"w")
pickle.dump(isdic,pifh2)






