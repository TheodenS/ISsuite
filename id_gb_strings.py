import sys
import os
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import pickle
import cStringIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIStandalone
import operator
from Bio.Blast.Applications import NcbiblastnCommandline
import argparse
import datetime

# Make a qualifiers dictionary for None hits
def makeNoHitDic():
    dic={}
    dic["hit_def"]="No_blast_hit"
    dic["sbjct_start"]="ND"
    dic["sbjct_end"]="ND"
    dic["expect"]="ND"
    dic["frame"]=(1,1)
    dic["bits"]="ND"
    dic["query"]="ND"
    dic["query_end"]="ND"
    dic["query_start"]="ND"
    dic["sbjct"]="ND"
    dic["score"]="ND"
    dic["strand"]="ND"
    dic["align_length"]="ND"
    return dic

def getdir(d):
	retlist=[]
	fileslist=os.listdir(d)
	for u in fileslist:
	    if u[0]==".":
                continue
            else:
		    retlist.append(u)
	return retlist




def doblast(isseq,sentevalue,blastdb):
    #print("entered doblast function")
    
    #print("Blasting seq\n"+ isseq)
    returnvalue=""
    #print tempfile
    writequery=open(tempfile,"w")
    writequery.write(str(isseq))
    writequery.close()

    # blast the sequence
    blastn_cline = NcbiblastnCommandline(query=tempfile, db=blastdb, evalue=sentevalue, outfmt=5, out=blastoutfiletemp,dust="no",task="blastn")
    cmd=blast_location+str(blastn_cline)
    os.system(cmd)
    #raw_input("cmdline blast")

    try:
        blastout_handle=open(blastoutfiletemp,"r")
        #raw_input("waiting after printing blast results")
        blastout_handle.seek(0,0)
        blast_records = list(NCBIXML.parse(blastout_handle))	
    except:
        raw_input("blast parse trouble in function getbestblast")

    # Check that we only have one blast record
    if len(blast_records)>1:
        raw_input("several blast records in blast for "+condfile+" st: "+sta+"  en: "+end)
    if len(blast_records)<1:
        raw_input("no blast records in blast for "+condfile+" st: "+sta+"  en: "+end)

    #collect all hsp hits in one file
    hsps=[]
    for record in blast_records:
        for alg in record.alignments:
            for hsp in alg.hsps:
                dic={}
                dic["hit_def"]=str(alg.hit_def)
                dic["sbjct_start"]=hsp.sbjct_start
                dic["sbjct_end"]=hsp.sbjct_end
                dic["expected"]=str(hsp.expect)
                dic["expect"]=hsp.expect
                dic["frame"]=hsp.frame
                dic["bits"]=hsp.bits
                dic["query"]=str(hsp.query)
                dic["query_end"]=hsp.query_end
                dic["query_start"]=hsp.query_start
                dic["sbjct"]=str(hsp.sbjct)
                dic["score"]=hsp.score
                dic["align_length"]=hsp.align_length
                hsps.append(dic)

       
    hsps.sort(key=operator.itemgetter("score"),reverse=True)
    if len(hsps)>0:
        return hsps[0]
    else:
        #print("len of hsps is not above 0 (no results)")
        return None


#blastdb="/Users/security/science/frankiaproj/frankiajul2014/reps_and_IS_db"

#blastoutfile="/Users/security/science/frankiaproj/frankiajul2014/blast_output_temp"


#dirpath="/Users/security/science/frankiaproj/frankiajul2014/parseresults/worked_files/"



parser = argparse.ArgumentParser()
parser.add_argument("-repeatsdb", help="dirwithall")
parser.add_argument("-tempfile", help="dirwithall")
parser.add_argument("-blastlocation", help="dirwithall")
parser.add_argument("-evalue", help="dirwithall")
parser.add_argument("-repsfa", help="dirwithall")
parser.add_argument("-gbfiles", help="dirwithall")
args=parser.parse_args()

tempfile=args.tempfile
blastoutfiletemp=args.tempfile+"idrepsrestemp"
blast_location=args.blastlocation
expectvalue=args.evalue


dirpath=args.gbfiles

genomerowcounter=0

#strend="\n"
#now = datetime.datetime.now()
#logfh=open(args.logpath,"w")
#text_log=str(now.year)+"_"+str(now.month)+"_"+str(now.day)+strend

if True:
    makedbstring=blast_location+"./makeblastdb -in "+args.repsfa+" -title "+"repsfasta"+" -out "+args.repeatsdb+" -hash_index -dbtype nucl"
    os.system(makedbstring)

    contigs=getdir(args.gbfiles)


    contigcounter=0
    for cunt in contigs:
        genomerowcounter+=1
        contigcounter+=1
  
        contig_dir_path=args.gbfiles+cunt

        gbfile_handle=open(contig_dir_path,"r")
        parsed_gbfile=list(SeqIO.parse(gbfile_handle,"genbank"))

        gbfile=parsed_gbfile[0]

        gbhitslist=[]

        for feature in gbfile.features:
	    fetype=feature.type
	    if fetype=="footprint":
		pass
	    else:
		continue
	

            genseq=gbfile.seq

            feature_start=int(min(feature.location.start.position,feature.location.end.position))
            feature_end=int(max(feature.location.start.position,feature.location.end.position))
            
            real_feature_start=feature_start+1
            real_feature_end=feature_end


            
            # biopython reports coords that can be used directly as list indices
            featureseq=genseq[feature_start:feature_end]


            firstsearch=[featureseq,real_feature_start,real_feature_end]
            
            remainsearches=[firstsearch]

            while len(remainsearches)>0:
                    
                currentsearch=remainsearches.pop()
                currentseq=currentsearch[0]
                searched_seq_genomic_start=currentsearch[1]
                searched_seq_genomic_end=currentsearch[2]
              	#print "got here" 
                doblast_results=doblast(currentseq,expectvalue,args.repeatsdb)

                # If no results, produce a noresult entry
                if doblast_results==None:
                    #print("blast results was None")
                    
                    listToAppend=[searched_seq_genomic_start,searched_seq_genomic_end,"nohit",makeNoHitDic()]
                    gbhitslist.append(listToAppend)


                #elif float(doblast_results["expect"])>1.0e-5:
                elif float(doblast_results["expect"])>expectvalue:
                    listToAppend=[searched_seq_genomic_start,searched_seq_genomic_end,"nohit",makeNoHitDic()]
                    gbhitslist.append(listToAppend)
                    
                else:


                    hitstart=min(int(doblast_results["query_start"]),int(doblast_results["query_end"]))
                    hitend=max(int(doblast_results["query_start"]),int(doblast_results["query_end"]))

                    found_hit_genomic_start=searched_seq_genomic_start+hitstart-1
                    found_hit_genomic_end=searched_seq_genomic_start+hitend-1


                    listToAppend=[found_hit_genomic_start,found_hit_genomic_end,doblast_results["hit_def"],doblast_results]
                    gbhitslist.append(listToAppend)

                    if not hitstart-1<=0:
                        leftremainseq=str(featureseq[0:hitstart-1])


                        left_remains_start=searched_seq_genomic_start
                        left_remains_end=searched_seq_genomic_start+hitstart-2

                        leftremains_list=[leftremainseq,left_remains_start,left_remains_end]

                        remainsearches.append([leftremainseq,left_remains_start,left_remains_end])

                    right_start=hitend
                    right_end=len(featureseq)
                    right_end=len(currentseq)

                    if not right_end-right_start<=0:
                        rightremainseq=str(featureseq[right_start:right_end])
                        
                        # The errors are probably here, where is the -1?
                        right_remains_start=currentsearch[1]+hitend
                        right_remains_end=currentsearch[2]
                        
                        rightremains_list=[rightremainseq,right_remains_start,right_remains_end]
                    
                        remainsearches.append([rightremainseq,right_remains_start,right_remains_end])


        featureslist=[]
        for feat in gbhitslist:
	    #print feat
            mystrand=feat[3]["frame"][0]*feat[3]["frame"][1]
            
            feature = SeqFeature(FeatureLocation(feat[0]-1,feat[1]), strand=mystrand,type=feat[2])


            feature.qualifiers=feat[3]
            featureslist.append(feature)

	for f in featureslist:
		gbfile.features.append(f)

        outfilepath=contig_dir_path
        #print(outfilepath)
        tempfh=open(outfilepath,"w")
        SeqIO.write([gbfile],tempfh,"genbank")
    #print "treated "+str(contigcounter)+ " contigs"


print "treated "+str(genomerowcounter)+ "genome rows"


#logfh.write(text_log)
#logfh.close()
