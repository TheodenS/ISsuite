import sqlite3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

def countweirds(myseq):
	weirds=0
	for lett in myseq:
		if not lett in "agctAGCT":
			weirds+=1
	return weirds

def countgc(myseq):
	gc=0
	for lett in myseq:
		if  lett in "gcGC":
			gc+=1
	return gc 


#connect to database
parser = argparse.ArgumentParser()
parser.add_argument("-fastafile", help="fastafile with seqs to be tested for ISISIS")
parser.add_argument("-databasefile", help="databasefile")
parser.add_argument("-databasename", help="databasename")
parser.add_argument("-renamedfastafile", help="databasename")
parser.add_argument("-sbjct_summary_csv", help="databasename")


args=parser.parse_args()

con = sqlite3.connect(args.databasefile)
cur = con.cursor()   


#"CREATE TABLE Gen2(Name TEXT, Seq TEXT,Length INT,GC FLOAT)")
dropstring="DROP TABLE IF EXISTS "+args.databasename
makestring="CREATE TABLE "+args.databasename+"("

makestring+="Name"+" "+"TEXT"+","
makestring+="Length"+" "+"INT"+","
makestring+="Seq"+" "+"TEXT"+","
makestring+="GC"+" "+"FLOAT"+","
makestring+="Weirds"+" "+"INT"+","
makestring+="Myid"+" "+"TEXT"+","

makestring=makestring[:-1]
makestring+=")"

with con:
    cur.execute(dropstring)
    cur.execute(makestring)




print "Reading fasta"
fastafile=args.fastafile
fastafh=open(fastafile,"r")
fastagroup=SeqIO.parse(fastafh,"fasta")
     		  
genomes_total_length=0
genomescount=0

csvstring="Name,Myid,Length,Seq,GC,Weirds\n"

fasta_recs_all=[]

for org in fastagroup:
	#print org
        #print "makingdb"
        genomescount+=1
        seq_len=len(org.seq)
        weirds=countweirds(str(org.seq))
        countGC=countgc(str(org.seq))
	insertstring="INSERT INTO "
	insertstring+=args.databasename
	insertstring+=" VALUES('"
	insertstring+=org.name
	insertstring+="',"
	insertstring+=str(seq_len)
	insertstring+=",'"+str(org.seq)+"',"+str(countGC/float(seq_len))+","+str(weirds)+","+"'sbjct_"+str(genomescount)+"')"


        csvstring+=org.name+","+"sbjct_"+str(genomescount)+","+str(seq_len)+","+str(org.seq)+","+str(countGC/float(seq_len))+","+str(weirds)+"\n"
	
	myid=str(genomescount)
	isname=myid
	isseq=str(org.seq)
       	s=Seq(isseq,IUPAC.IUPACUnambiguousDNA())
        newrec=SeqRecord(s)
        newrec.id="sbjct_"+str(myid)
        newrec.description=""
        fasta_recs_all.append(newrec)

        with con:
        	cur.execute(insertstring)
                #cur.execute("CREATE INDEX sbjct_db_index on "+args.databasename+" (Myid)")

fasta_fh=open(args.renamedfastafile,"w")
SeqIO.write(fasta_recs_all,fasta_fh,"fasta")
fasta_fh.close()

summary_csvh=open(args.sbjct_summary_csv,"w")
summary_csvh.write(csvstring)
print "done making genomedb"



