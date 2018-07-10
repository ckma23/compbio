import os
import time

class HbPlusProcesser(object):
    def hbplusprocessedreader(self,hborvdw):
        os.system('mkdir ~/bioresearch/compbio/files_wip/hbplus_sorted_%s_files_baseset' %hborvdw)
        # os.chdir(os.path.expanduser('~/bioresearch/hbplusprocessed%sfiles' %hborvdw))
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_processed_%s_files_baseset' %hborvdw))
        listofprocessedfiles = os.listdir('.')
        print listofprocessedfiles
        # listofprocessedfiles=["pdb1mnb.hb2","pdb1b2m.hb2"]
        # listofprocessedfiles=["pdb1mnb.hb2"]
      # listofprocessedfiles=["pdb1mnb.nb2","pdb1b2m.nb2"]
        for hbplusprocessedfile in listofprocessedfiles:
            print hbplusprocessedfile
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_processed_%s_files_baseset' %hborvdw))
            info1=[]
            info2=[]
            info3=[]
            info4=[]
            info5=[]
            info6=[]
            info7=[]
            info8=[]
            individualfileobject=open(hbplusprocessedfile)
            lines=individualfileobject.readlines()
            i=8 #set the counter to skip the 8 lines of the header
            totallinesinfile=len(lines) #sum up the total lines in each file
        # totallinesinfile=20
            while i < totallinesinfile: #iterate through the file for each line in it
        #    print lines[i]
        # if info8="HH" then sum up the HH
                # print lines[i]
                store=''.join(map(str,lines[i]))
        #    print store
        #    print len(store)
                i+=1
        #    info1.append(store[0:6])
        #    info2.append(store[6:9])
        #    info3.append(store[9:12])
        #    info4.append(store[14:20])
        #    info5.append(store[20:23])
        #    info6.append(store[24:27])
        #upon investigation of this output hbplus output comes out in a single line and can only sepearte values by position in column.
                # print store
                # print store[0:6]
                # print store[6:9]
                # print store[9:13]
                # print store[14:20]
                # print store[20:23]
                # print store[24:27]
                # print store[28:32]
                # print store[33:35]
           #this is a typical expected output.
           #A0003-ARG NH2 A0014-ARG O   3.37 SM  11  5.39 161.4  2.39 147.7 153.2    46

           #A0003-
           #ARG
           # NH2
           #A0014-
           #ARG
           #O
           #3.37
           #SM
           #B0027-  C O2' B0028-  C O4' 3.32 HH  -2 -1.00 140.9  2.48 154.5  99.0    45

           #B0027-
           #  C
           # O2'
           #B0028-
           #  C
           #O4'
           #3.32
           #HH

           # we may want to be careful here because the O2' may deliberately have spaces.
           #first value is the amino acid or rna residue number
           #second value is thhe amino acid or rna atom can see that we want to strip out the empty spaces to put this into a csv
           #third value is the atom type can see that we want to strip out the empty spaces to put this into a csv
           #fourth value is the amino acid or rna residue number
           #fifth value is the amino acid or rna reside can see that we want to strip out the empty spaces to put this into a csv
           #sixth value is the amino acid or rna atom we can see that we want to strip otu the empty spaces

                #not storing these at the moment
                info7.append(store[28:32])
                info8.append(store[33:35])

           # want to switch it up here to if store is in the set and store is in the set. This will speed this up by alot
                aminoacidlist=set(["ARG","ALA","ASN","ASP","GLN","GLU","GLY","CYS","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"])
                nucleotidebase=set(["U","G","C","A"])
                for aa in aminoacidlist:
                    for nb in nucleotidebase:
                # print aa
                # print nb
                # if store[6:9].strip() == aa and store[20:23].strip() == nb:
                        if (store[6:9].strip() == aa and store[20:23].strip() == nb):
                            info4.append(store[0:6].strip())
                            info5.append(store[6:9].strip())
                            info6.append(store[9:13].strip())
                            info1.append(store[14:20].strip())
                            info2.append(store[20:23].strip())
                            info3.append(store[24:27].strip())
                #appears that we need to account if the amino acid or the nucleotide base is either side. however we need to fix the files so hbplus sorted always have nucelotide base is always on first column.
                        elif (store[6:9].strip() == nb and store[20:23].strip() == aa):
                            info1.append(store[0:6].strip())
                            info2.append(store[6:9].strip())
                            info3.append(store[9:13].strip())
                            info4.append(store[14:20].strip())
                            info5.append(store[20:23].strip())
                            info6.append(store[24:27].strip())
            self.hbplusstorefilewriter(hbplusprocessedfile,info1,info2,info3,info4,info5,info6,hborvdw)
      # for i in range(len(info1)):
      #   print "%s %s %s %s %s %s" %(info1[i],info2[i],info3[i],info4[i],info5[i],info6[i])

    def hbplusstorefilewriter(self,proteinname,col1,col2,col3,col4,col5,col6,hborvdw):
    # def hbplushbstorefilewriter(proteinname,col1,col2,col3,col4,col5,col6):
    #   os.chdir("/Users/curtisma/bioresearch")
    #   os.system('mkdir hbplus%ssortedfiles' %hborvdw)
      # os.system('mkdir hbplushbsortedfiles')
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_sorted_%s_files_baseset' %hborvdw))
    #   os.chdir("/Users/curtisma/bioresearch/hbplus%ssortedfiles" %hborvdw)
        filenamestring="%s.%ssorted" %(proteinname,hborvdw)
        os.system("rm %s" %filenamestring)
        file = open(filenamestring,"a")
        print "%s has been created" %filenamestring
        for i in range(len(col1)):
            #we need to clean the col1 and col 2 to make hbplus compatibile to dssr i.e. B0018- A..... means B strand Argenine18"
            # thus we need to clean it to be A0018 on the B strand"
            # remove the first letter from col because that is the strand information and we need the residue type information
            removedstrandstring=col1[i][1:]

            cleanednb=col2[i]+removedstrandstring.strip()
            # combine this into one and it'l be A00018-
            strandnb=col1[i]
            # still keep the strand information!
            strandnb=strandnb[0]
            #strip the space and the -
            cleanednb=self.hbplustodssrnbstringcleaner(cleanednb)
            file.write("%s,%s,%s,%s,%s,%s,%s\n" %(cleanednb,strandnb,col3[i],col4[i],col5[i],col6[i],hborvdw))
        # file.write("%s %s %s %s %s %s %s\n" %(cleanednb,strandnb,col3[i],col4[i],col5[i],col6[i],hborvdw))

    def hbplustodssrnbstringcleaner(self,nbtobecleaned):
      nbclean=nbtobecleaned
      nbclean=nbclean.strip()
      # nbclean=nbclean.replace("B","")
      nbclean=nbclean.strip('-')
      return nbclean

    def hbplushbandvdwcombiner(self):
        print "combining the hbplushb sorted and hbplusvdw sorted into one file"
        os.system("mkdir ~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_baseset")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_sorted_hb_files_baseset'))
    #   os.chdir("/Users/curtisma/bioresearch")
    #   os.system("mkdir hbplushbvdwcombined")
        # os.chdir("/Users/curtisma/bioresearch/hbplushbsortedfiles")
        listofprocessedhbplusfiles = os.listdir('.')
        listofprocessedhbplusfiles = ["1mnb.bondcategorized"]
        print listofprocessedhbplusfiles
        for i in listofprocessedhbplusfiles:
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_sorted_hb_files_baseset'))
            # os.chdir("/Users/curtisma/bioresearch/hbplushbsortedfiles")
            print os.getcwd()
            individualfileobject=open(i)
            lhs,rhs=i.split(".",1)
            print lhs
            filenamestring="%s.hbplushbvdwsorted" %(lhs)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_baseset'))
            # os.chdir("/Users/curtisma/bioresearch/hbplushbvdwcombined")
            os.system("rm %s" %filenamestring)
            file = open(filenamestring,"a")
            for line in individualfileobject:
                file.write(line)
                # os.chdir("/Users/curtisma/bioresearch/hbplusvdwsortedfiles")
                os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_sorted_vdw_files_baseset'))
                vdwfile="%s.nb2.vdwsorted" %(lhs)
                individualfileobject=open(vdwfile)
                filenamestring="%s.hbplushbvdwsorted" %(lhs)
                os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_baseset'))
                # os.chdir("/Users/curtisma/bioresearch/hbplushbvdwcombined")
                file = open(filenamestring,"a")
            for line in individualfileobject:
                file.write(line)
