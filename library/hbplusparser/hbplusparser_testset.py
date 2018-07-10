import os
import time

class HbPlusProcesserTestSet(object):
    # this is going to take in either a hb or vdw file from HBplus and process it to a format that will be used later on.
    # this time lets give it the file directory.

    def hbplusprocessed_file_prepper_reader(self,folder_path,hborvdw):
        os.system('mkdir ~/bioresearch/compbio/files_wip/hbplus_sorted_%s_files_testset' %hborvdw)
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s' %folder_path))
        list_of_processed_protein_files = os.listdir('.')
        for protein in list_of_processed_protein_files:
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_sorted_%s_files_testset' %hborvdw))
            os.system('mkdir %s'%protein)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s/%s' %(folder_path,protein)))
            list_of_pose_files = os.listdir('.')
            for pose in list_of_pose_files:
                os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s/%s' %(folder_path,protein)))
                info1=[]
                info2=[]
                info3=[]
                info4=[]
                info5=[]
                info6=[]
                info7=[]
                info8=[]
                individualfileobject=open(pose)
                lines=individualfileobject.readlines()
                i=8 #set the counter to skip the 8 lines of the header
                totallinesinfile=len(lines) #sum up the total lines in each file
                # print len(lines)
                while i < totallinesinfile: #iterate through the file for each line in it
                    store=''.join(map(str,lines[i]))
                    i +=1
                    info7.append(store[28:32])
                    info8.append(store[33:35])
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
                            # print store[9:13]
                self.hbplusstorefilewriter(pose,info1,info2,info3,info4,info5,info6,hborvdw,protein)

    def hbplustodssrnbstringcleaner(self,nbtobecleaned):
        nbclean=nbtobecleaned
        nbclean=nbclean.strip()
        # nbclean=nbclean.replace("B","")
        nbclean=nbclean.strip('-')
        return nbclean

    def hbplusstorefilewriter(self,pose_name,col1,col2,col3,col4,col5,col6,hborvdw,proteinname):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_sorted_%s_files_testset/%s' %(hborvdw,proteinname)))
        filenamestring="%s.%ssorted" %(pose_name,hborvdw)
        # need to remove the preivous file so it's not adding to a preivously made file
        os.system("rm %s" %filenamestring)
        file = open(filenamestring,"a")
        # print "%s has been created" %filenamestring
        for i in range(len(col1)):
        #we need to clean the col1 and col 2 to make hbplus compatibile to dssr i.e. B0018- A..... means B strand Argenine18"
        # thus we need to clean it to be A0018 on the B strand"
        # remove the first letter from col because that is the strand information and we need the residue type information
            removedstrandstring = col1[i][1:]
            cleanednb = col2[i]+removedstrandstring.strip()
        # combine this into one and it'l be A00018-
            strandnb = col1[i]
        # still keep the strand information!
            strandnb = strandnb[0]
        #strip the space and the -
            cleanednb = self.hbplustodssrnbstringcleaner(cleanednb)
            file.write("%s,%s,%s,%s,%s,%s,%s\n" %(cleanednb,strandnb,col3[i],col4[i],col5[i],col6[i],hborvdw))
        file.close()

    def hbplushbandvdwcombiner(self):
        print "combining the hbplushb sorted and hbplusvdw sorted into one file"
        os.system("mkdir ~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_testset")
        #lets just use the command cat hb.txt >> hb_vdw_combined.txt
        #lets just use the command cat vdw.txt >> hb_vdw_combined.txt
        # this is much easier accomplished in bash.
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_sorted_hb_files_testset'))
        proteins = os.listdir('.')
        for protein in proteins:
            os.system("mkdir ~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_testset/%s" %protein)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_sorted_hb_files_testset/%s' %protein))
            hb_files = os.listdir('.')
            #lets grab the first Complex_*string
            for hbfile_to_be_merged in hb_files:
                #the file name is this "Complex_997g.nb2.vdwsorted", strip out the Complex_997g
                complex_name,rhs=hbfile_to_be_merged.split(".",1)
                os.system("cat ~/bioresearch/compbio/files_wip/hbplus_sorted_hb_files_testset/%s/%s >  ~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_testset/%s/%s" %(protein,hbfile_to_be_merged,protein,complex_name))
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_sorted_vdw_files_testset/%s' %protein))
            vdw_files = os.listdir('.')
            for vdwfile_to_be_merged in vdw_files:
                complex_name,rhs=vdwfile_to_be_merged.split(".",1)
                os.system("cat ~/bioresearch/compbio/files_wip/hbplus_sorted_vdw_files_testset/%s/%s >> ~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_testset/%s/%s" %(protein,vdwfile_to_be_merged,protein,complex_name))


    #   os.chdir("/Users/curtisma/bioresearch/hbplushbsortedfiles")
    #   listofprocessedhbplusfiles = os.listdir('.')
    #   print listofprocessedhbplusfiles
    #   for i in listofprocessedhbplusfiles:
    #     os.chdir("/Users/curtisma/bioresearch/hbplushbsortedfiles")
    #     print os.getcwd()
    #     individualfileobject=open(i)
    #     lhs,rhs=i.split(".",1)
    #     print lhs
    #     filenamestring="%s.hbplushbvdwsorted" %(lhs)
    #     os.chdir("/Users/curtisma/bioresearch/hbplushbvdwcombined")
    #     os.system("rm %s" %filenamestring)
    #     file = open(filenamestring,"a")
    #     for line in individualfileobject:
    #       file.write(line)
    #     os.chdir("/Users/curtisma/bioresearch/hbplusvdwsortedfiles")
    #     vdwfile="%s.nb2.vdwsorted" %(lhs)
    #     individualfileobject=open(vdwfile)
    #     filenamestring="%s.hbplushbvdwsorted" %(lhs)
    #     os.chdir("/Users/curtisma/bioresearch/hbplushbvdwcombined")
    #     file = open(filenamestring,"a")
    #     for line in individualfileobject:
    #       file.write(line)
