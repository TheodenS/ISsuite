import xlrd
import xlrd
import argparse
import sqlite3
import unicodedata
from Bio import SeqIO

#infile="/Users/security/science/first/bigcsv.csv"
#work_csv="/Users/security/science/RNA/transcriptome_Theo/annotation_all.filtered.taxgrps.stats.xlsx"


#db="/Users/security/science/RNA/metatrans.annotation_all.filtered.taxgrps.stats.2.db"
#
#con = sqlite3.connect(db)
#cur = con.cursor()
#
#
#selstring="SELECT sum(GS670_3p0) FROM annotation WHERE \"group\"='Cyanobacteria'"
##selstring="SELECT sum(GS670_3p0) FROM annotation WHERE \"group\"='Archaeplastida'"
##selstring="SELECT orf_id  FROM annotation WHERE \"group\"='Archaeplastida' AND GS667_3p0>0"
#
#with con:
#        cur.execute(selstring)
#
#reply=cur.fetchall()
#for u in reply:
#    print u[0]


        #	cur.execute('''SELECT Kind FROM repeats WHERE Name=?''', (name,))
#	furrl=cur.fetchone()
#	contigseq=str(furrl[0])
#
#	return contigseq



#metadata_csv="/Users/security/science/RNA/transcriptome_Theo/annotation_all.filtered.taxgrps.stats.xlsx"
#
#
#
#
##print is_contig_place_dic.keys()
#
#book = xlrd.open_workbook(metadata_csv)
#print "the number of worksheets is", book.nsheets
#print "Worksheet name(s):", book.sheet_names()
#sh = book.sheet_by_index(0)
##nn=deepcopy(sh)
##pickle.dump(nn,rh)
##sh=pickle.load(rh)
#print sh.name, sh.nrows, sh.ncols
#print "Cell D30 is", sh.cell_value(rowx=0, colx=0)
#print "sheet has "+str(sh.nrows)+" rows"
#
#print "Worksheet name(s):", book.sheet_names()
#sh = book.sheet_by_index(0)
##pickle.dump(nn,rh)
##sh=pickle.load(rh)
#print sh.name, sh.nrows, sh.ncols
#print "Cell D30 is", sh.cell_value(rowx=0, colx=0)
#print "sheet has "+str(sh.nrows)+" rows"
#
## make headline
#headline=""
#for row in range(0,sh.nrows):
#    val=str(sh.cell_value(rowx=row, colx=22)).replace(",","_comma_")
#    print val
