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

## BLAST location 
blastlocation="/Users/security/science/frankiaproj/frankiajul2014/software/ncbi-blast-2.2.29+/bin/"
metadata_csv="/Users/security/science/RNA/transcriptome_Theo/annotation_all.filtered.taxgrps.stats.xlsx"

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
place_contig_dic_location=basedir+"place_contig_dic.pydic"
place_contig_csv_location=basedir+"place_contig_csv.csv"
iscounts_contigcsv=basedir+"iscounts_contigcsv.csv"
allgenomesinfo=basedir+"allgenomesinfo.csv"
foundishistogram=basedir+"foundis_histogram.eps"
joinedcsv=basedir+"joined.csv"
iscontiglist=basedir+"iscontiglist.picl"
isdicout=basedir+"isdicout.picl"


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

# copy the original fastafile to basegenomes
if not os.path.isfile(genomefastafile):
	cpline="cp "+ genomefastafile_base+ " " +genomefastafile
	os.system(cpline)
        report("copied "+genomefastafile_base +" to "+genomefastafile+"\n",text_results_file)
# copy is file to current dir
if not os.path.isfile(queryfastafile):
	cpline="cp "+ args.queryfastafile+ " " +queryfastafile
	os.system(cpline)
	report("copied "+args.queryfastafile +" to "+   queryfastafile+"\n",text_results_file)

# Make genome database and a renamed fasta file
if not os.path.isfile(genomedatabase):
        cmstr="python /Users/security/science/software/ISsuite/make_genomedb_and_renamed_fasta.py -fastafile="+genomefastafile+" -databasefile="+genomedatabase+" -databasename="+"genomeinfo" +" -renamedfastafile="+genomefasta_dir+"renamed.fa"+" -sbjct_summary_csv="+basedir+"subject_summary.csv"
        os.system(cmstr)

# Make csv summarizying query
if not os.path.isfile(basedir+"query_summary.csv"):
    cmstr_querycsv="python /Users/security/science/software/ISsuite/make_query_csv.py -fastafile="+queryfastafile+" -query_summary_csv="+basedir+"query_summary.csv"
    os.system(cmstr_querycsv)

if not os.path.isfile(basedir+"query_histogram.pdf"):
    cmd_make_query_fasta="python /Users/security/science/software/ISsuite/makehistogram.py -in_csv="+basedir+"query_summary.csv"+" -out="+basedir+"query_histogram.pdf"
    os.system(cmd_make_query_fasta)

if not os.path.isfile(basedir+"subject_histogram.pdf"):
    cmd_make_sbjct_fasta="python /Users/security/science/software/ISsuite/makehistogram.py -in_csv="+basedir+"subject_summary.csv"+" -out="+basedir+"subject_histogram.pdf"
    os.system(cmd_make_sbjct_fasta)


prog="/Users/security/science/software/ISsuite/"+"blast_query_vs_sbject.py"

resultname="rvg"
contigxmlresultspath=blast_output_location+resultname+".xml"

runstring="python "+prog+" -genomeblastdb "+genomedbfilebase+" -blastlocation "+blastlocation+" -blastoutput "+contigxmlresultspath+" -evalue "+str(float(args.maxeval)*10)
runstring+=" -queryfile " + args.queryfastafile
runstring+=" -genomefasta " +genomefasta_dir+"renamed.fa" 

os.system(runstring)

report("ran blast_reps_vs_genomes_fasta.py\n",text_results_file)

prog="/Users/security/science/software/ISsuite/"+"parse_qu_vs_sb_to_gbstring.py"


runstring="python "+prog+" -blastoutput "+contigxmlresultspath
runstring+=" -gbfilesdir "+gbfilesdir
runstring+=" -evalue "+str(args.maxeval)
runstring+=" -genomedb "+genomedatabase
os.system(runstring)


prog="/Users/security/science/software/ISsuite/"+"id_gb_strings.py"

runstring="python "+prog+" -gbfiles "+gbfilesdir +" -tempfile "+helperfiles_dir+"idrepstemp"+" -blastlocation "+blastlocation+" -evalue "+str(args.maxeval)
runstring+=" -repsfa "+args.queryfastafile
runstring+=" -repeatsdb "+genomedbfilebase+"repsblastbase"
print runstring
os.system(runstring)


#prog="/Users/security/science/software/ISsuite/"+"parse_oldnewtret2.py"
prog="/Users/security/science/software/ISsuite/"+"followparweolfgnewtre2.py"
#prog="/Users/security/science/software/ISsuite/"+"newtrets2.py"
runstring="python "+prog
runstring+=" -gbfilesdir="+gbfilesdir
runstring+=" -allgenomesinfo="+allgenomesinfo
runstring+=" -genomedatabase="+genomedatabase
runstring+=" -maxeval="+str(args.maxeval)
runstring+=" -iscontiglist="+iscontiglist
runstring+=" -isdicout="+isdicout
runstring+=" -iscsv="+basedir+"query_summary.csv"

print runstring
os.system(runstring)


prog="/Users/security/science/software/ISsuite/"+"foundislengthhistogram.py"
runstring="python "+prog 
runstring+=" -in_csv "+allgenomesinfo
runstring+=" -out "+foundishistogram
os.system(runstring)


# Make dictionary with number of contigs found at each site
# use treatxls.py
#if True:
if not os.path.isfile(place_contig_dic_location):
    prog="/Users/security/science/software/ISsuite/"+"treatxls.py"
    runstring="python "+prog
    runstring+=" -place_contig_dic="+place_contig_dic_location
    runstring+=" -place_contig_csv="+place_contig_csv_location
    runstring+=" -metadata_csv="+metadata_csv
    metadata_csv
    print runstring

    os.system(runstring)


#if True:
if not os.path.isfile(iscounts_contigcsv):
#if not os.path.isfile("/Users/security/science/ccc2.csv"):
    #prog="/Users/security/science/software/ISsuite/"+"count_no_transcripts.py"
    prog="/Users/security/science/software/ISsuite/"+"masturbate.py"
    runstring="python "+prog
    runstring+=" -place_contig_dic="+place_contig_dic_location
    #runstring+=" -database="+"/Users/security/science/trydic2.csv"
    runstring+=" -outcsv="+iscounts_contigcsv
    print runstring
    os.system(runstring)

prog="/Users/security/science/software/ISsuite/"+"join_csvs.py"
runstring="python "+prog
runstring+=" -iscounts_contigcsv="+iscounts_contigcsv
runstring+=" -outcsv="+joinedcsv
print runstring

os.system(runstring)

prog="/Users/security/science/software/ISsuite/"+"modify_csv_joindepths.py"
runstring="python "+prog
runstring+=" -incsv="+joinedcsv
runstring+=" -outcsv="+joinedcsv+"m.csv"
os.system(runstring)



prog="/Users/security/science/software/ISsuite/"+"r_environgraph.py"
runstring="python "+prog
runstring+=" -in_csv="+joinedcsv
runstring+=" -in_csv_all="+joinedcsv+"m.csv"
runstring+=" -out="+basedir+"plt.pdf"
os.system(runstring)

prog="/Users/security/science/software/ISsuite/"+"r_ggpairs.py"
runstring="python "+prog
runstring+=" -in_csv="+basedir+"joinedmod.csv"
runstring+=" -out="+basedir+"environmentalpairs.png"
os.system(runstring)

quit()




print "finished program"



