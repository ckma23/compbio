from Bio import *                                           #use the BioPython Python library
from Bio.PDB import *                                       #more specifically import the BioPython PDB library
import os
import numpy as np
import math


class FTDockChainStripping(object):
    def pdbfile_splitter_rna_protein(self,protein_seperated_files_directory,rna_seperated_files_directory,stored_pdbfiles_database):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip'))
        os.system("mkdir %s" %protein_seperated_files_directory)
        os.system("mkdir %s" %rna_seperated_files_directory)
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s' %protein_seperated_files_directory))
        #clean up the directory so it is done again
        os.system("rm *.pdb")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s' %rna_seperated_files_directory))
        #clean up the directory so it is done again
        os.system("rm *.pdb")

        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s' %stored_pdbfiles_database))

        list_of_pdb_files = os.listdir('.')
        for pdbfile in list_of_pdb_files:
            lhs,rhs=pdbfile.split(".",1)
            lhs=lhs[3:]
            print lhs
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s'%stored_pdbfiles_database))
            with open(pdbfile) as pdblines:
                # need to double over unfortunately at the moment to get it to terminate on the first chain.
                for line in pdblines:
                    filenamestring="%s_%s.pdb" %(lhs,line[21])
                    # line[21] is the chain name
                    if line[17:20] in set(["ARG","ALA","ASN","ASP","GLN","GLU","GLY","CYS","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]):
                        if line[0:4] == "ATOM":
                            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s' %protein_seperated_files_directory))
                            file = open(filenamestring,"a")
                            file.write(line)
                        elif line[0:4] == "TER ":
                            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s' %protein_seperated_files_directory))
                            file = open(filenamestring,"a")
                            file.write(line)
                            break
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s'%stored_pdbfiles_database))
            with open(pdbfile) as pdblines:
                for line in pdblines:
                    filenamestring="%s_%s.pdb" %(lhs,line[21])
                    if line[17:20].strip() in set(["A","U","C","G"]):
                        if line[0:4] == "ATOM":
                            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s' %rna_seperated_files_directory))
                            file = open(filenamestring,"a")
                            file.write(line)
                        elif line[0:4] == "TER ":
                            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s' %rna_seperated_files_directory))
                            file = open(filenamestring,"a")
                            file.write(line)
                            break

    def preprocess_ftdock(self,rna_or_protein):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/'))
        os.system('mkdir %s_seperated_pdbfiles_preprocessperl_testset' %rna_or_protein)
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s_seperated_pdbfiles_preprocessperl_testset' %rna_or_protein))
        os.system('rm *.parsed')
        os.system('rm *.fasta')
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s_seperated_pdbfiles_testset' %rna_or_protein))
        list_of_pdb_files = os.listdir('.')
        os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/ftdock-2-dev2/scripts-2.0.3'))
        for pdbfile in list_of_pdb_files:
            os.system("/usr/bin/perl preprocess-pdb.perl -pdb /Users/curtisma/bioresearch/compbio/files_wip/%s_seperated_pdbfiles_testset/%s" %(rna_or_protein,pdbfile))
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s_seperated_pdbfiles_testset' %rna_or_protein))
        os.system("mv -i *.parsed ~/bioresearch/compbio/files_wip/%s_seperated_pdbfiles_preprocessperl_testset" %rna_or_protein)
        os.system("mv -i *.fasta ~/bioresearch/compbio/files_wip/%s_seperated_pdbfiles_preprocessperl_testset" %rna_or_protein)

    def pdb_file_combine_rms_calc(self):
        print "Combining the separated pdb file rna chain and protein chain"
        os.system("mkdir ~/bioresearch/compbio/files_wip/combined_pdbfiles_testset")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset'))
        proteins = os.listdir('.')
        for protein in proteins:
            protein =  protein[0:4]
            print protein
            os.system("cat ~/bioresearch/compbio/files_wip/rna_seperated_pdbfiles_testset/%s* >  ~/bioresearch/compbio/files_wip/combined_pdbfiles_testset/%s_combined.pdb" %(protein,protein))
            os.system("cat ~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset/%s* >> ~/bioresearch/compbio/files_wip/combined_pdbfiles_testset/%s_combined.pdb" %(protein,protein))
        #lets match based on the protein pdb name sample output: 1jbr_A.pdb


    def hbplushbandvdwcombiner(self):
        print "combining the hbplushb sorted and hbplusvdw sorted into one file"
        os.system("mkdir ~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_testset")
        #lets just use the command cat hb.txt >> hb_vdw_combined.txt
        #lets just use the command cat vdw.txt >> hb_vdw_combined.txt
        # this is much easier accomplished in bash.
        # os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_sorted_hb_files_testset'))
        # proteins = os.listdir('.')
        # for protein in proteins:
        #     os.system("mkdir ~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_testset/%s" %protein)
        #     os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_sorted_hb_files_testset/%s' %protein))
        #     hb_files = os.listdir('.')
        #     #lets grab the first Complex_*string
        #     for hbfile_to_be_merged in hb_files:
        #         #the file name is this "Complex_997g.nb2.vdwsorted", strip out the Complex_997g
        #         complex_name,rhs=hbfile_to_be_merged.split(".",1)
        #         os.system("cat ~/bioresearch/compbio/files_wip/hbplus_sorted_hb_files_testset/%s/%s >  ~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_testset/%s/%s" %(protein,hbfile_to_be_merged,protein,complex_name))
        #     os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_sorted_vdw_files_testset/%s' %protein))
        #     vdw_files = os.listdir('.')
        #     for vdwfile_to_be_merged in vdw_files:
        #         complex_name,rhs=vdwfile_to_be_merged.split(".",1)
        #         os.system("cat ~/bioresearch/compbio/files_wip/hbplus_sorted_vdw_files_testset/%s/%s >> ~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_testset/%s/%s" %(protein,vdwfile_to_be_merged,protein,complex_name))
        #



    def file_copier_to_sjsu_cluster_testset(self,baseset_or_testset):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s_complexes_pdb' %baseset_or_testset))
        os.system('scp *.ent cma@spartan02.sjsu.edu:bioresearch/compbio/files_wip/%s_complexes_pdb' %baseset_or_testset)

    def file_copier_to_sjsu_cluster(self,rna_or_protein):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s_seperated_pdbfiles_preprocessperl_testset' %rna_or_protein))
        os.system('scp *.parsed cma@spartan02.sjsu.edu:curtisma/bioresearch/%s_pre_processedperl_seperated_pdbfiles_testset' %rna_or_protein)
