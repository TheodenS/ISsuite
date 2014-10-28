# blast_reps_vs_genomes_5aug_2012.py
# 5 aug
# based on Copied 141/home/theo/Workfiles_nov/programs/blast_reps_vs_genomes/blast_reps_vs_genomomes.py modiefied mar 2012
# copy of this program in 141/home/theo/july2012/reps_vs_genomes_collected_progs/

import os
from Bio.Blast.Applications import NcbiblastnCommandline
import argparse
import datetime
import time
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import sqlite3
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import datetime
import time



parser = argparse.ArgumentParser()

parser.add_argument("-genomecsv", help="directory of genome fastas")
parser.add_argument("-blastoutput", help="blastoutput")
parser.add_argument("-evalue", help="blastoutput")
parser.add_argument("-gbfilesdir", help="blastoutput")
parser.add_argument("-genomedb", help="fastaorcsv")


args=parser.parse_args()

gbfilesdir=args.gbfilesdir


blast_output_location=args.blastoutput



def treatstring(strl,rrecord):

	# REMEMBER, strl uses list numbering
	featlist=[]
	waitingfor="st"
	cpos=0

	# from position strl[0] and so on to len-1
	for onepos in strl:
		
		if waitingfor=="st":
			if onepos=="x":
				print "position start at "+str(cpos)
				# featstart is in list notation
				featst=cpos+1
				waitingfor="en"

		# Last position can not be stopped by a "-" after it
		# if we are at the end of the list
		# in a started feature
		elif waitingfor=="en" and cpos==len(strl)-1:
				
				# end is last position in list	
				featen=cpos

				# Make a feature to the last nuc in genome
				# list values have to be converted to gb values
				# This is where genbank might decide to fuck up the files
				# in this case the solution might be the line of 23nov,
				#feature = SeqFeature(FeatureLocation(featst-1,featen), strand=1,type="footprint")

				# Genbank magic adds one to start position?

				feature = SeqFeature(FeatureLocation(featst-1,featen+1), strand=1,type="footprint")
				featlist.append(feature)
		# Normal ending of x stretch
		else:
			if onepos=="-":
				# last x was position before this one
				featen=cpos-1
				#print "position end at "+str(cpos)

				# See explanation on tirade above for explanation.
				# line as of 23nov was
				#feature = SeqFeature(FeatureLocation(featst-1,featen-1), strand=1,type="footprint")

				# modification from list to gb	
				feature = SeqFeature(FeatureLocation(featst-1,featen+1), strand=1,type="footprint")
				featlist.append(feature)
				waitingfor="st"
		cpos+=1
		
	for f in featlist: 
		rrecord.features.append(f)

	return rrecord

strend="\n"
#now = datetime.datetime.now()
#logfh=open(args.logpath,"w")
#text_log=str(now.year)+"_"+str(now.month)+"_"+str(now.day)+strend

mineval=float(args.evalue)
hitcount=0
xmlfh=open(blast_output_location,"r")

genomesdic={}

#if feature.type in contig_repeat_dic.keys():
#                        contig_repeat_dic[feature.type].append([hit_number,cont,feature,hitseq,abb])
#                    else:
#                        contig_repeat_dic[feature.type]=[[hit_number,cont,feature,hitseq,abb]]
try:
	blast_records = list(NCBIXML.parse(xmlfh))	
except:
	raw_input("no file")
		#print isxmlfilename1+"/"+aa,
		#continue



con = sqlite3.connect(args.genomedb)
cur = con.cursor()

blast_recordcounter=0		
blast_alg_counter=0
hsp_counter=0
if len(blast_records)>1:
	for rec in blast_records:
                blast_recordcounter+=1
		rename=str(rec.query)
                #print rec.query
                #print "rename="+rename
		algs=rec.alignments
#if len(blast_records)<1:
#	print isxprint("no records in gb file")
#	raw_input("no records in xml file")
		algs=rec.alignments


        	for alg in algs:
                        blast_alg_counter+=1
        	        for hsp in alg.hsps:
        	                hexpect=float(hsp.expect)
                                #print alg.hit_def
                                hsp_counter+=1
        	        
        	                # Start of treatment for algs making the cut	

        	                if hexpect<mineval:
					#print alg.hit_def +" was under eval"
					hitset=[rename,hsp];
					if alg.hit_def in genomesdic.keys():
                        			genomesdic[alg.hit_def].append(hitset)
						#print "it was"
                        		else:
						genomesdic[alg.hit_def]=[hitset]
					hitcount+=1
					#print rename
					#print alg.id
#
print "records, alignemnets, hsps"
print blast_recordcounter		
print blast_alg_counter
print hsp_counter
#contigxmlresultspath=blast_output_location
#
#
#xmloutpath=contigxmlresultspath
#blastn_cline = NcbiblastnCommandline(query=is_fa, db=genomedbfile, evalue=evalu, outfmt=5, out=contigxmlresultspath,max_target_seqs=1000,dust="no",task="blastn")
#newbline=blast_location+str(blastn_cline)
#os.system(newbline)
print "len genomesdic keys"
print len(genomesdic.keys())
raw_input("len")

gen_in_genomesdic=0
for gen in genomesdic.keys():
        
        gen_in_genomesdic+=1
	allfeatures=[]

        #print gen

	#print gen 
	cur.execute('''SELECT Length,Seq,Name FROM genomeinfo WHERE Myid=?''', (gen,))
	furrl=cur.fetchone()
	genomelen=int(furrl[0])
	contigseq=str(furrl[1])
	myname=str(furrl[2])
        print "correct name"
        print myname

	totalstring=["-"]*genomelen

	bestscorehsp=genomesdic[gen][0][1]
	bestexphsp=genomesdic[gen][0][1]

	for hitl in genomesdic[gen]:
		isname=hitl[0]
		myhsp=hitl[1]
		framepart1=int(myhsp.frame[0])
		framepart2=int(myhsp.frame[1])

		if framepart1!=1:
			raw_input("frampart 1 negative")
#							
		framestrand=framepart1*framepart2
						
		#gb positions

		st=int(myhsp.sbjct_start)
		en=int(myhsp.sbjct_end)
        	                

		abs_st=min(st,en)
		abs_en=max(st,en)
		#print "x ing from "+str(abs_st-1)+" to "+str(abs_en)
						 
		totalstring[abs_st-1:abs_en]=["x"]*(1+abs_en-abs_st)


#isname=isxmlfile.split(".")[0]
		#isname=aa.split(".")[0]

		# Make feature for record file 
		# st and en are unmodified from the blast file
		# they can be backwards
		# original row as of nov23 read
		#feature = SeqFeature(FeatureLocation(st-1,en), strand=framepart2,type="repeathit")
		# start location should not be modified by -1 in this case


		# changed Featurelocation start, older line:
		#feature = SeqFeature(FeatureLocation(st,en), strand=framepart2,type="repeathit")

		# Genbank magic adds one to start position?

		feature = SeqFeature(FeatureLocation(st-1,en), strand=framepart2,type="repeathit")
		#feature.id=rename
		#feature.qualifiers["hit"]=rename
		feature.id=isname
		feature.qualifiers["hit"]=isname
		feature.qualifiers["sequence_length"]=str(len(contigseq[myhsp.sbjct_start-1:myhsp.sbjct_end]))
#
		feature.qualifiers["score"]=str(myhsp.score)
		if myhsp.score>bestscorehsp.score:
			bestscorehsp=myhsp

		if myhsp.expect<bestexphsp.expect:
			bestexphsp=myhsp

		# hit length is disance from hits start plus hit end
		feature.qualifiers["hitlength"]=str(1+max(st,en)-min(st,en))
#
		feature.qualifiers["expect"]=str(myhsp.expect)
		feature.qualifiers["query_start"]=str(myhsp.query_start)
		feature.qualifiers["query_end"]=str(myhsp.query_end)
		feature.qualifiers["sbjct_start"]=str(myhsp.sbjct_start)
		feature.qualifiers["sbjct_end"]=str(myhsp.sbjct_end)
		############
		# Kept for later
		################

		#if hsp.sbjct_start<hsp.sbjct_end:
		#	iscsv.append(len(str(contigseq[hsp.sbjct_start-1:hsp.sbjct_end])))
		#	iscsv.append(str(contigseq[hsp.sbjct_start-1:hsp.sbjct_end]))
		#allcsvs.append(str(contigseq[hsp.sbjct_start-1:hsp.sbjct_end]))
		#allcsvs.append(len(str(contigseq[hsp.sbjct_start-1:hsp.sbjct_end])))
		#else:
		#	turnseq=contigseq[hsp.sbjct_end-1	:hsp.sbjct_start]
		#	tseq=Seq(str(turnseq),IUPAC.IUPACAmbiguousDNA())
		#	neseq=tseq.reverse_complement()
		#	iscsv.append(len(str(neseq)))
		#	iscsv.append(str(neseq))
		#allcsvs.append(len(str(neseq)))
		#allcsvs.append(str(neseq))


		allfeatures.append(feature)

        	# End of treating hsp in dic
	# End of gen
	s=Seq(str(contigseq),IUPAC.IUPACUnambiguousDNA())
	newrec=SeqRecord(s)
	newrec.id=gen
	newrec.description="is search on "+args.blastoutput
	newrec.features=allfeatures

	#first_hits_gb_record_dir=gbfilesdir+gen
	if not(os.path.isdir(gbfilesdir)):
    		os.makedirs(gbfilesdir) 
	first_hits_gb_record_file=gbfilesdir+gen+".gb"
	newrec= treatstring(totalstring,newrec)
	gbrecord_fh=open(first_hits_gb_record_file,"w")
	SeqIO.write([newrec],gbrecord_fh,"genbank")
	gbrecord_fh.close()
