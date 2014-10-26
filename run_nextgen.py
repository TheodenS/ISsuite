import os
import os.path
import argparse
import datetime
import sqlite3

# Appends text to file
def report(mytext,myfile):
     restext=open(myfile,"a")
     restext.write(str(now.year)+"_"+str(now.month)+"_"+str(now.day)+"\n")
     restext.write(str(now.hour)+":"+str(now.minute)+"\n")
     restext.write(mytext+"\n")
     restext.close()

# Parse args
parser = argparse.ArgumentParser()
parser.add_argument("-basedir", help="dirwithall")
parser.add_argument("-basefastafile", help="dirwithall")
parser.add_argument("-queryfastafile", help="dirwithall")
parser.add_argument("-maxeval", help="dirwithall")
args=parser.parse_args()

basedir=args.basedir

now = datetime.datetime.now()
text_results_file=basedir+"runprogs_results_"+str(now.year)+"_"+str(now.month)+"_"+str(now.day)+".txt"

# Report stuff about run
firstdesc= ""
firstdesc+= "running runprogs"+"\n"
firstdesc+= "basedir "+args.basedir+"\n"
firstdesc+= "basefastafile "+args.basefastafile+"\n"
firstdesc+= "queryfastafile "+args.queryfastafile+"\n"
firstdesc+= "maxeval "+str(args.maxeval)+"\n"


# Make dirs
if not os.path.isdir(args.basedir):
     os.makedirs(args.basedir)

report(firstdesc,text_results_file)

# BLAST location 
# Files location
###############
blastlocation="/Users/security/science/frankiaproj/frankiajul2014/software/ncbi-blast-2.2.29+/bin/"
metadata_csv="/Users/security/science/RNA/transcriptome_Theo/annotation_all.filtered.taxgrps.stats.xlsx"
environmentaldata="/Users/security/science/RNA/mimebs_metadata/metadata.2009.csv"

genomedbfilebase=basedir+"genomeblastdbs/"
blast_output_location=basedir+"reps_vs_genomes_blast_output/"
helperfiles_dir=basedir+"helperfiles/"
genomefasta_dir=basedir+"fastafiles/"
genomefastafile_base=args.basefastafile
genomefastafile=genomefasta_dir+"genomes_org.fa"
queryfastafile=genomefasta_dir+"queries.fa"
genomedatabase_dir=basedir+"databases/"
genomedatabase=genomedatabase_dir+"genomeinfo.db"
gbfilesdir=basedir+"gb_files/"
place_contig_dic_location="/Users/security/science/place_contig_dic.pydic"
#place_contig_dic_location=basedir+"place_contig_dic.pydic"
place_contig_csv_location=basedir+"place_contig_csv.csv"
iscounts_contigcsv=basedir+"iscounts_contigcsv.csv"
#allgenomesinfo=basedir+"allgenomesinfo.csv"
contigs_with_is_csv=basedir+"contigs_with_is.csv"
foundishistogram=basedir+"foundis_histogram.eps"
joinedcsv=basedir+"joined.csv"
iscontiglist=basedir+"iscontiglist.picl"
isdicout=basedir+"isdicout.picl"

# Make directories
####################
if not os.path.isdir(helperfiles_dir):
     os.makedirs(helperfiles_dir)

if not os.path.isdir(genomedbfilebase):
     os.makedirs(genomedbfilebase)
#
if not os.path.isdir(blast_output_location):
     os.makedirs(blast_output_location)

if not os.path.isdir(genomefasta_dir):
     os.makedirs(genomefasta_dir)

if not os.path.isdir(genomedatabase_dir):
        os.makedirs(genomedatabase_dir)

if not os.path.isdir(gbfilesdir):
     os.makedirs(gbfilesdir)

# copy the original fastafile to folder
###################
if not os.path.isfile(genomefastafile):
	cpline="cp "+ genomefastafile_base+ " " +genomefastafile
	os.system(cpline)
        report("copied "+genomefastafile_base +" to "+genomefastafile+"\n",text_results_file)

# copy is file to current dir
###################
if not os.path.isfile(queryfastafile):
	cpline="cp "+ args.queryfastafile+ " " +queryfastafile
	os.system(cpline)
	report("copied "+args.queryfastafile +" to "+   queryfastafile+"\n",text_results_file)

# Make genome database and a renamed fasta file
# database named databases/genomeinfo.db/genomeinfo
# with columns
#Name
#Length
#Seq
#GC 
#Weirds
#Myid
# and 
#subject_summary.csv
# with name, myid, length,seq,weirds,gc
##################
if not os.path.isfile(genomedatabase):
        cmstr="python /Users/security/science/software/ISsuite/make_genomedb_and_renamed_fasta.py -fastafile="+genomefastafile+" -databasefile="+genomedatabase+" -databasename="+"genomeinfo" +" -renamedfastafile="+genomefasta_dir+"renamed.fa"+" -sbjct_summary_csv="+basedir+"subject_summary.csv"

        cmstr+=" -minseqlen " + "30"
        cmstr+=" -maxseqlen " + "1000000"
        cmstr+=" -abb_prefix" + "subj_"

        print "making contig info database and renamed fasta file (propername->myid) "+genomedatabase
        os.system(cmstr)

# Make csv summarizying query (for histogram)
# Name,Myid,Length,Seq,GC,Weirds
###################
if not os.path.isfile(basedir+"query_summary.csv"):
    cmstr_querycsv="python /Users/security/science/software/ISsuite/make_query_csv.py -fastafile="+queryfastafile+" -query_summary_csv="+basedir+"query_summary.csv"

    cmstr_querycsv+=" -renamedqueryfile " + basedir+"renamed_query.fa"
    cmstr_querycsv+=" -minseqlen " + "1000"
    cmstr_querycsv+=" -maxseqlen " + "2000"


    os.system(cmstr_querycsv)
    print "Made csv file for queries"

# Make histograms
####################
print "query info"
if not os.path.isfile(basedir+"query_histogram.eps"):
    cmd_make_query_histo="python /Users/security/science/software/ISsuite/makehistogram.py -in_csv="+basedir+"query_summary.csv"+" -out="+basedir+"query_histogram.eps"
    cmd_make_query_histo+=" -title " + "Queries_lengths"
    os.system(cmd_make_query_histo)
    print "making histogram of queries"

print "subject info"
if not os.path.isfile(basedir+"subject_histogram.eps"):
    cmd_make_sbject_histo="python /Users/security/science/software/ISsuite/makehistogram.py -in_csv="+basedir+"subject_summary.csv"+" -out="+basedir+"subject_histogram.eps"
    cmd_make_sbject_histo+=" -title " + "Subjects_lengths"
    print "making historgram of subject seqs"
    os.system(cmd_make_sbject_histo)

# Set destination name for blast results
###########
resultname="rvg"
contigxmlresultspath=blast_output_location+resultname+".xml"

# Run blastn
############
prog="/Users/security/science/software/ISsuite/"+"blast_query_vs_sbject.py"
if not os.path.isfile(blast_output_location+resultname+".xml"):
    runstring="python "+prog+" -genomeblastdb "+genomedbfilebase+" -blastlocation "+blastlocation+" -blastoutput "+contigxmlresultspath+" -evalue "+str(float(args.maxeval)*10)
    #runstring+=" -queryfile " + args.queryfastafile
    runstring+=" -queryfile " + basedir+"renamed_query.fa"
    runstring+=" -genomefasta " +genomefasta_dir+"renamed.fa" 
    print "running blast outputting "+blast_output_location+resultname+".xml"
  #  os.system(runstring)
    report("ran blast_reps_vs_genomes_fasta.py\n",text_results_file)


# Parse blastn results to genbank files
# using databases/genomeinfo.db/genomeinfo
################
prog="/Users/security/science/software/ISsuite/"+"parse_qu_vs_sb_to_gbstring.py"
runstring="python "+prog+" -blastoutput "+contigxmlresultspath
runstring+=" -gbfilesdir "+gbfilesdir
runstring+=" -evalue "+str(args.maxeval)
runstring+=" -genomedb "+genomedatabase
#os.system(runstring)

# Identify the footprints by blasting against the queries
#################
prog="/Users/security/science/software/ISsuite/"+"id_gb_strings.py"
runstring="python "+prog+" -gbfiles "+gbfilesdir +" -tempfile "+helperfiles_dir+"idrepstemp"+" -blastlocation "+blastlocation+" -evalue "+str(args.maxeval)
#runstring+=" -repsfa "+args.queryfastafile
runstring+=" -repsfa "+args.basedir+"renamed_query.fa"
runstring+=" -repeatsdb "+genomedbfilebase+"repsblastbase"
print runstring
#os.system(runstring)

# Make a dictionary and csv connecting the stations to the subject seqs 
# mydic[place][sbjct]=count
# not normally run, as it uses a dictionary in
# place_contig_dic_location="/Users/security/science/place_contig_dic.pydic
#################
if not os.path.isfile(place_contig_dic_location):
    prog="/Users/security/science/software/ISsuite/"+"treatxls.py"
    runstring="python "+prog
    runstring+=" -place_contig_dic="+place_contig_dic_location
    runstring+=" -place_contig_csv="+place_contig_csv_location
    runstring+=" -metadata_csv="+metadata_csv
    print runstring
    os.system(runstring)


# outputs detainledcsv, with
# isname,abbcontig,propercontig,loc_start,loc_end,loc_strand,score,expected,sbjct_start,sbject_end,sbjctseq,queryseq
# islongcsv is a csv with extra is information, from isfinder
# iscontiglist, a list of contigs with queries in them
# isdicout, an isname-keyed dic with feature, gbrecord, name of contig
###################
prog="/Users/security/science/software/ISsuite/"+"parse_gb_files.py"
runstring="python "+prog
runstring+=" -gbfilesdir="+gbfilesdir
runstring+=" -place_contig_dic="+place_contig_dic_location
runstring+=" -contigs_with_is_csv="+contigs_with_is_csv
runstring+=" -genomedatabase="+genomedatabase
runstring+=" -maxeval="+str(args.maxeval)
runstring+=" -iscontiglist="+iscontiglist
runstring+=" -isdicout="+isdicout
runstring+=" -iscsv="+basedir+"query_summary.csv"
runstring+=" -iscsvlong="+basedir+"iscsvlong.csv"
runstring+=" -iscsvdetailed="+basedir+"iscsvdetailed.csv"

print runstring
#os.system(runstring)

# Make histogram of found pieces
####################
prog="/Users/security/science/software/ISsuite/"+"foundislengthhistogram.py"
runstring="python "+prog 
runstring+=" -in_csv "+contigs_with_is_csv
runstring+=" -out "+foundishistogram
#os.system(runstring)


# outputs list with numkey,placename,station,depth,total transcripts, hittranscripts,fraction
#as iscounts_contigcsv.csv
####################### 
if not os.path.isfile(iscounts_contigcsv):
    #prog="/Users/security/science/software/ISsuite/"+"count_no_transcripts.py"
    prog="/Users/security/science/software/ISsuite/"+"masturbate.py"
    runstring="python "+prog
    runstring+=" -place_contig_dic="+place_contig_dic_location
    runstring+=" -hitlist="+iscontiglist
    runstring+=" -outcsv="+iscounts_contigcsv
    print runstring
    #os.system(runstring)


#produces csv joined.csv,with environmental data added to stations
####################
prog="/Users/security/science/software/ISsuite/"+"join_csvs.py"
runstring="python "+prog
runstring+=" -iscounts_contigcsv="+iscounts_contigcsv
runstring+=" -outcsv="+joinedcsv
print runstring
#os.system(runstring)

# produces a list of stations, where filter sizes have been joined together
####################
prog="/Users/security/science/software/ISsuite/"+"modify_csv_joindepths.py"
runstring="python "+prog
runstring+=" -incsv="+joinedcsv
runstring+=" -outcsv="+joinedcsv+"-joinedfilters.csv"
#os.system(runstring)
print runstring

# Make graph with environment variables and hit counts
# using joined.csv as in_csv
# and joinedfilters.csv as in_csv_all
####################
prog="/Users/security/science/software/ISsuite/"+"r_environgraph.py"
runstring="python "+prog
runstring+=" -in_csv="+joinedcsv
runstring+=" -plotfraction="+"False"
runstring+=" -in_csv_joined="+joinedcsv+"-joinedfilters.csv"
runstring+=" -out="+basedir+"plt.pdf"
#os.system(runstring)
print runstring

prog="/Users/security/science/software/ISsuite/"+"r_ggpairs.py"
runstring="python "+prog
runstring+=" -in_csv="+joinedcsv
runstring+=" -out="+basedir+"environmentalpairs.png"
#os.system(runstring)
print runstring





contigs_with_is_count_csv=basedir+"contigs_with_is_count.csv"



# outputs list with contigs with ISs in them.
prog="/Users/security/science/software/ISsuite/"+"more_more_join_csvs.py"
runstring="python "+prog
runstring+=" -detailed_csvf="+basedir+"iscsvdetailed.csv"
runstring+=" -iscsvlongf="+basedir+"iscsvlong.csv"
runstring+=" -outcsv="+contigs_with_is_count_csv

#print runstring
#os.system(runstring)

prog="/Users/security/science/software/ISsuite/"+"moremoremorejoincsv.py"
runstring="python "+prog
#runstring+=" -detailed_csvf="+basedir+"iscsvdetailed.csv"
runstring+=" -detailed_csvf="+contigs_with_is_csv
runstring+=" -iscsvlongf="+basedir+"iscsvlong.csv"
runstring+=" -contigxls_isinfo_join_csv="+basedir+"fullist.csv"

#print runstring

#os.system(runstring)
# joins list of iscontigs with is hits to extra info about is
# used to output 5oct2-2.csv" 
prog="/Users/security/science/software/ISsuite/"+"join_is_info_to_results.py"
runstring="python "+prog
runstring+=" -is_contigs_count_csv="+basedir+"fullist.csv"
#runstring+=" -is_contigs_count_csv="+contigs_with_is_count_csv
runstring+=" -long_csv="+basedir+"iscsvlong.csv"
runstring+=" -out_csv_allcontigs_is_long="+basedir+"out_csv_allcontigs_is_long.csv"
#print runstring
#os.system(runstring)

# Makes directories with info per organism group
# takes5oct2-2.csv"as input 
prog="/Users/security/science/software/ISsuite/"+"dna_per_organism.py"
runstring="python "+prog
runstring+=" -is_contigs_count_csv="+basedir+"fullist.csv"
#runstring+=" -is_contigs_count_csv="+contigs_with_is_count_csv
runstring+=" -group_out_dir="+"bactgroups"
#print runstring
#os.system(runstring)

#takes science/5oct2-2.csv as input
prog="/Users/security/science/software/ISsuite/"+"get_transcripts_orgn_place.py"
runstring="python "+prog
runstring+=" -is_contigs_count_csv="+basedir+"fullist.csv"
runstring+=" -is_per_group_csv="+basedir+"is_groups.csv"
#print runstring
#os.system(runstring)

#prog="/Users/security/science/software/ISsuite/"+"r_bars_normalized_location_groups.py"
#runstring="python "+prog
#runstring+=" -in_csv="+basedir+"is_groups.csv"
#print runstring
#runstring+=" -out="+basedir+"group_frequencies.pdf"
#os.system(runstring)


if not os.path.isfile(basedir+"appended_large.db"):
    prog="/Users/security/science/software/ISsuite/"+"makedb_from_resultscsv.py"
    runstring="python "+prog
    runstring+=" -fullist="+basedir+"fullist.csv"
    runstring+=" -outdb="+basedir+"appended_large.db"
    #runstring+=" -outdb="+"/users/security/science/dest2.db"

    os.system(runstring)

#if not os.path.isfile(basedir+"appended_large.db"):
if True:
    prog="/Users/security/science/software/ISsuite/"+"wrangle_results_sigA.py"
    runstring="python "+prog

    runstring+=" -in_db="+basedir+"appended_large.db"
    runstring+=" -out="+basedir+"bajs"
    runstring+=" -outcsv="+basedir+"summary.csv"
    runstring+=" -environmentaldata="+environmentaldata

    os.system(runstring)
# wrangle_results.py uses this
#    csvh=open("/Users/security/Desktop/headers2.csv","r")

print "finished program"



