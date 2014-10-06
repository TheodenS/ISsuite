import os

def getdir(d):
	retlist=[]
	fileslist=os.listdir(d)
	for u in fileslist:
	    if u[0]==".":
                continue
            else:
		    retlist.append(u)
	return retlist


endi="/Users/security/science/metagen2009_2/"

if not os.path.isdir(endi):
     os.makedirs(endi)

dnadir="/Users/security/science/mimebs_metagenome_2009_unassembled/"



dnacsv="/Users/security/science/genomes2009.csv"

dnaline=""

dirs=getdir(dnadir)

for di in dirs:
    if ".gz" in di:
        dnaline+=di+"\n"

cc=open(dnacsv,"w")
cc.write(dnaline)
cc.close()

cc2=open(dnacsv,"r")
cread=cc2.read()
clines=cread.split("\n")


notlit=[]
for u in clines:
    if u=="":
        continue
    dones=getdir(endi)


    ufoldername=u.replace(".fna.gz","")
    uzipname=u.replace(".fna.gz",".zip")
    if ufoldername in dones:
        continue
        
        #print ufoldername+" was in dones"
        #print u+"already done"
    else:
        print ufoldername+" was NOT in dones"

        if uzipname in dones:
            pass
            #print u+"in ZIP already done"
        else:
            print "nondone"
            print u
            print "nondone"
            cpline="cp "+dnadir+u+" "+endi+u
            print cpline
            os.system(cpline)
            untline="gunzip "+endi+u
            os.system(untline)

            removeline="rm "+endi+u.replace(".fna.gz","")+".fna"
            print removeline
            #os.system(removeline)

            runline ="python /Users/security/science/software/ISsuite/modifyfasta_for_dna.py"
            runline +=" -incsv="+endi+u.replace(".gz","")
            runline +=" -outcsv="+endi+u.replace(".gz","")+"clean"
            os.system(runline)

            removeline2="rm "+endi+u.replace(".gz","")
            print removeline2
            os.system(removeline2)


            runline="python issuite_rerun30nov.py -basedir=/Users/security/science/metagen2009_2/"+u.replace(".fna.gz","").replace(".fna","")+"/ -basefastafile="+endi+u.replace(".gz","")+"clean"+" -queryfastafile=/Users/security/science/RNA/ISonly3630st_less3k.fa -maxeval=0.00001"
            print runline
            os.system(runline)
            
            removeline3="rm "+endi+u.replace(".gz","")+"clean"
            os.system(removeline3)

            zipline="zip "+endi+u.replace(".fna.gz","")
            os.system(zipline)

            

