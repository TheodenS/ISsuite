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
parser.add_argument("-query_summary_csv", help="databasename")


args=parser.parse_args()


#"CREATE TABLE Gen2(Name TEXT, Seq TEXT,Length INT,GC FLOAT)")
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
        genomescount+=1
        seq_len=len(org.seq)
        weirds=countweirds(str(org.seq))
        countGC=countgc(str(org.seq))

        csvstring+=org.name+","+"query_"+str(genomescount)+","+str(seq_len)+","+str(org.seq)+","+str(countGC/float(seq_len))+","+str(weirds)+"\n"
	

summary_csvh=open(args.query_summary_csv,"w")
summary_csvh.write(csvstring)



