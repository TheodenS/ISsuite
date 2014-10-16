import os
import pickle

def getdir(d):
	retlist=[]
	fileslist=os.listdir(d)
	for u in fileslist:
	    if u[0]==".":
                continue
            else:
		    retlist.append(u)
	return retlist


indir="/Users/security/science/metagen2009_2/"

folders=getdir(indir)

outputdir="/Users/security/science/genome2009/"

allrow=""
rowlist=[]


county=0

namedic={}
for folder in folders:
    if "done.tar.gz" in folder:

        #county+1

        #if county>2:
        #    continue
        print folder
        ungzcom="tar -xzf " +indir+folder+" -C /Users/security/science/"
        print ungzcom
        os.system(ungzcom)

        untarred_folder="/Users/security/science/Users/"

        untarred_folder_fh=open("/Users/security/science/Users/security/science/metagen2009_2/"+folder.replace("_done.tar.gz","")+"/subject_summary.csv","r")
        untarred_read=untarred_folder_fh.read()
        splitt=untarred_read.split("\n")
        splitt=splitt[1:]
        for u in splitt:
            if u=="":
                continue
            uu=u.split(",")

            loc=folder.replace("_done.tar.gz","")

            if not uu[1] in namedic.keys():
                namedic[uu[1]]=[loc]
            else:
                print loc
                namedic[uu[1]].append(loc)
        
         
            #outline+=loc+","+uu[0]+","+uu[2]+"\n"
            
        
        untarred_folder="/Users/security/science/Users/"
        rmline="rm -R "+untarred_folder
        os.system(rmline)

yy=open("/Users/security/science/pickleddic.dic","w")
pickle.dump(namedic,yy)
        #destdir="/Users/security/science/Users/security/science/metagen2009_2/"+folder.replace("_done.tar.gz","")

        #subjectcsv=destdir
        #csvdir=destdir+"/subject_summary.csv"
        #print csvdir

        #newplace=outputdir+folder.replace("_done.tar.gz","")+".csv"
        #print newplace
        #copyrow="cp "+csvdir+" "+newplace
        #print copyrow
        #os.system(copyrow)
        

#
#
#
#        sbjctsfh=open(csvdir,"r")
#        sbjctread=sbjctsfh.read()
#        sbjctsplit=sbjctread.split("\n")
#        header=sbjctsplit[0]
#        sbjctsplit=sbjctsplit[1:]
#
#        o=0
#        for row in sbjctsplit:
#            print o
#            o+=1
#            
#            allrow=folder.replace("_done.tar.gz","")+","+row+"\n"
#            if row in rowlist:
#                raw_input("already have row "+row)
#            else:
#                rowlist.append(row)
#
#        
#
#
#oo=open("/Users/security/science/longlist.csv","w")
#oo.write(allrow)
