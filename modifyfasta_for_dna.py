from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-incsv", help="dirwithall")
parser.add_argument("-outcsv", help="dirwithall")
args=parser.parse_args()

infh=open(args.incsv,"r")
#infh=open("/Users/security/science/correctstations6/joined.csv","r")

#fa_file="/Users/security/science/RNA/"+"ISonly3630st.fa"
fa_file_out=args.outcsv
#fa_file_out="/Users/security/science/"+"GS659_0p1_mod.fna"

#fa_file="/Users/security/science/"+"GS659_0p1.fna"

#fa_fh=open(fa_file,"r")
seqs=SeqIO.parse(infh,"fasta")

newrecx=[]
c=0
for s in seqs:

    myseq=str(s.seq)
    myseq_clean=myseq.replace("N","")

    #if len(s.seq)>100:
    newseq=Seq(myseq_clean,IUPAC.IUPACUnambiguousDNA())
    s.seq=newseq


    newrecx.append(s)


tempfh=open(fa_file_out,"w")
SeqIO.write(newrecx,tempfh,"fasta")
tempfh.close()
print c
