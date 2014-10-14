import sqlite3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import xlrd


def makeinsertline(inlist):
    adcolst=""
    if len(inlist)>0:
        c=0
        for col2 in inlist:

            c+=1

            typp=gettype(col2)
            typp2=type(typp).__name__
            #mess=str(sh.cell_value(rowx=row, colx=col2))
            mess=str(col2)
            
            
            typp="TEXT"
            if typp=="TEXT":
                #messg=mess.encode('ascii','ignore')
                messg=mess
                messg=messg.replace("'","APOS").replace(",","_comma_")
                sqstart="'"
                sqend="'"

            elif mess=="":
                messg=mess.encode('ascii','ignore')
                messg=messg.replace("'","APOS").replace(",","_comma_")
                sqstart="'"
                sqend="'"
            else:
                messg=mess.encode('ascii','ignore')
                messg=messg.replace("'","APOS").replace(",","_comma_")
                sqstart="'"
                sqend="'"
            adcolst+=sqstart+messg+sqend+","

    while len(adcolst.split(","))<100:
        adcolst+="'novalue',"
        #print adcolst
        #print c
        #print "ll"
    return adcolst

#def longerstr(inst):
    
    


def gettype(n):
	my="TEXT"
	if n in['group_reads' , 'group_orfs']:
		my="INT"
	if n in['best_hit_percent_identity','contig_length' , 'orf_length'] :
		my="FLOAT"
	if n[0:2]=='GS' or n[0:2]=='LD':
		my="INT"
	return my

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

def getxlxinfo(sh,name):
    for row in range(1,sh.nrows):
	fullname=sh.cell_value(rowx=row, colx=0)
            

        retlist=[]
        if fullname==name:
            for c in range(0,sh.ncols):
                retlist.append(str(sh.cell_value(rowx=row, colx=c)))
            return retlist
    return retlist

def gettabinfo(name):
    for u in tabsplit: 
        sp=u.split("\t")
        if sp[0]==name:
            return sp 
    return []

#connect to database
parser = argparse.ArgumentParser()
parser.add_argument("-fullist", help="fastafile with seqs to be tested for ISISIS")
parser.add_argument("-outdb", help="databasefile")
#parser.add_argument("-databasename", help="databasename")
#parser.add_argument("-renamedfastafile", help="databasename")
args=parser.parse_args()



#metadata_xlsx="/Users/security/science/RNA/transcriptome_Theo/annotation_all.filtered.taxgrps.stats.xlsx" 
#metadata_tab="/Users/security/science/RNA/transcriptome_Theo/annotation_all.com_fil_upd_low_nob.tab" 
alldata=args.fullist
#alldata="/Users/security/science/transcripts2009is/fullist.csv" 
database="allist"
destdb=args.outdb
#destdb="/Users/security/science/dest2.db" 


alldatafh=open(alldata,"r")
allread=alldatafh.read()
allsplit=allread.split("\n")
#print len(allsplit)
headers=allsplit[0]
allsplit=allsplit[1:]
allsplit=allsplit[:-1]
#:for s in allsplit:
 #   print s[:-1]
#print allsplit

splithead=headers.split(",")
##print headers
#print splithead
#print len(splithead)
#book = xlrd.open_workbook(metadata_xlsx)
#print "the number of worksheets is", book.nsheets
#print "Worksheet name(s):", book.sheet_names()
#sh = book.sheet_by_index(0)
#print sh.name, sh.nrows, sh.ncols
#print "Cell D30 is", sh.cell_value(rowx=0, colx=0)
#print "sheet has "+str(sh.nrows)+" rows"

#booktabfh = open(metadata_tab,"rU")
#tabread=booktabfh.read()
#tabsplit=tabread.split("\n")

makestring_xls=""
nonline=""
hades=""
for v in splithead:
    #v=sh.cell_value(rowx=0, colx=col)
    if v=="group":
        print "renamed"
        v="group_"
    ty="TEXT"
    #ty=gettype(v)
    makestring_xls+=v+" "+ty+","
    hades+=v+","
    nonline+="'',"
#print makestring_xls

#print makestring_xls

makestring_tab=""
#tabtabs=tabsplit[0].split("\t")

#print "========tabs headers"

notabinfoline=""
#for co in tabtabs:
#    ty="KOOOOO"
#    if co=="group":
#        co="group_"
#    ty=gettype(co)
#    makestring_tab+="TAB_" + str(co) + " " + ty + ","
#    notabinfoline+="'',"

#print "tab-------------"
#print makestring_tab
#print "-------------"


#print "long-------------"
#longst=makestring_xls+makestring_tab
longst=makestring_xls
#print longst
#print "indlong-------------"

	#makestring=makestring[:-1]
	#makestring+=")"

#makestring=longmakestring[:-1]
longst=longst.replace("-","_dash_")

#args=parser.parse_args()

con = sqlite3.connect(destdb)
#con = sqlite3.connect(args.databasefile)
cur = con.cursor()   

#metacon = sqlite3.connect("/Users/security/science/RNA/test8copy.db")

#"CREATE TABLE Gen2(Name TEXT, Seq TEXT,Length INT,GC FLOAT)")
dropstring="DROP TABLE IF EXISTS "+database
makestring="CREATE TABLE "+database+"("
#dropstring="DROP TABLE IF EXISTS "+args.databasename
#makestring="CREATE TABLE "+args.databasename+"("

#makestring+="Name"+" "+"TEXT"+","
#makestring+="Length"+" "+"INT"+","
#makestring+="Seq"+" "+"TEXT"+","
#makestring+="GC"+" "+"FLOAT"+","
#makestring+="Weirds"+" "+"INT"+","
#makestring+="Myid"+" "+"INT"+","





makestring+=longst

makestring=makestring[:-1]
makestring+=")"
#print makestring
#print len(makestring.split(","))


with con:
     	cur.execute(dropstring)
     	cur.execute(makestring)


#print "Reading fasta"
#fastafile=args.fastafile
#fastafh=open(fastafile,"r")
#fastagroup=SeqIO.parse(fastafh,"fasta")
#     		  
#genomes_total_length=0
#genomescount=0
#
#fasta_recs_all=[]
#
#
#
iterationcount=0
for org in allsplit:
        #print org
        splitrwo=org.split(",")
        if org=="":
            continue
        iterationcount+=1
        #if iterationcount>100:
        #    continue
        ##print "STARTORG//////////////////////"
	#print org.name
	#genomestr=""
        #genomescount+=1
        #seq_len=len(org.seq)
        #weirds=countweirds(str(org.seq))
        #countGC=countgc(str(org.seq))
	insertstring="INSERT INTO "
	insertstring+=database
	insertstring+="("+headers.replace(",group,",",group_,")+") "
	#insertstring+=args.databasename
	insertstring+=" VALUES('"
	#insertstring+=org.name
	#insertstring+="',"
	#insertstring+=str(seq_len)
	#insertstring+=",'"+str(org.seq)+"',"+str(countGC/float(seq_len))+","+str(weirds)+","+str(genomescount)+","

        #print org

        restrow=org[:-1].split(",")
        #restrow=getxlxinfo(allsplit,org.name)
        if len(restrow)>0:
            xls_insertrow=makeinsertline(restrow)
        else:
            print "happneddssakj"
            xls_insertrow=nonline

        #print xls_insertrow
        #print xls_insertrow 

        #newst=longerstr()
        #print xls_insertrow



        
        xlinsert=xls_insertrow.replace("-","_dash_")
        #insertstring+=xls_insertrow
        #insertstring+=xlinsert
        insertstring+=xlinsert[1:]



        insertstring= insertstring[:-1]

	insertstring+=")"

        #print insertstring
        #print len(insertstring.split(","))
        with con:
        	#print insertstring.replace("-","_dash_")
        	cur.execute(insertstring)

#fasta_fh=open(args.renamedfastafile,"w")
#SeqIO.write(fasta_recs_all,fasta_fh,"fasta")
#fasta_fh.close()
print "done making genomedb"



