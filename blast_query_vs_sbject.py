# blast_reps_vs_genomes_5aug_2012.py
# 5 aug
# based on Copied 141/home/theo/Workfiles_nov/programs/blast_reps_vs_genomes/blast_reps_vs_genomomes.py modiefied mar 2012
# copy of this program in 141/home/theo/july2012/reps_vs_genomes_collected_progs/

import os
from Bio.Blast.Applications import NcbiblastnCommandline
import argparse
import datetime
import time

sep=","
end="\n"
text_log=""

now = datetime.datetime.now()

parser = argparse.ArgumentParser()

parser.add_argument("-genomefasta", help="directory of genome fastas")
parser.add_argument("-genomeblastdb", help="location of blast database with genomes")
parser.add_argument("-blastlocation", help="location of blast prog")
parser.add_argument("-blastoutput", help="blastoutput")
parser.add_argument("-evalue", help="blastoutput")
parser.add_argument("-queryfile", help="fastaorcsv")


args=parser.parse_args()
genomedbdir=args.genomeblastdb
blast_location=args.blastlocation
genomedbfile= genomedbdir+"genomedb"

#logfh=open(args.logpath,"w")
text_log+=str(now.year)+"_"+str(now.month)+"_"+str(now.day)+end

blast_output_location=args.blastoutput


evalu=float(args.evalue)



makedbstring=blast_location+"./makeblastdb -in "+args.genomefasta+" -title "+"genomefasta"+" -out "+genomedbfile+" -hash_index -dbtype nucl"
print makedbstring

os.system(makedbstring)

			
resultname="rvg"
contigxmlresultspath=args.blastoutput


xmloutpath=contigxmlresultspath
blastn_cline = NcbiblastnCommandline(query=args.queryfile, db=genomedbfile, evalue=evalu, outfmt=5, out=contigxmlresultspath,max_target_seqs=1000,dust="no",task="blastn")
newbline=blast_location+str(blastn_cline)
print newbline
os.system(newbline)
#logfh.write(text_log)
#logfh.close()



