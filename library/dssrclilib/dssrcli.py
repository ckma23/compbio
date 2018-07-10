import os as os #import the python os library
import csv as csv #import the python csv library
import re #provides regular expression matching
import json

class dssrclihelper(object):
    # def dssrcli(self):
    #     os.chdir('/Users/curtisma/bioresearch/compbio/pdbfiles')
    #     listoffiles = os.listdir('.')
    #     os.chdir('/Users/curtisma/bioresearch/')
    #     os.system('mkdir DSSRprocessedfiles')
    #     for f in listoffiles:
    #         stringprep=('~/bioresearch/compbio/pdbfiles/%s' %f)
    #         os.chdir('/Users/curtisma/bioresearch/')
    #         lhs,rhs=f.split(".",1)
    #         os.system('./x3dna-dssr input=%s -o=%s.json -json' %(stringprep,lhs))
    #         # os.system('./x3dna-dssr input=%s' %stringprep)
    #         pdblhs,pdbrhs=f.split(".",1)
    #         # print pdblhs
    #         os.system('mv %s.json /Users/curtisma/bioresearch/DSSRprocessedfiles' %lhs)
    #         os.system('mv dssr* /Users/curtisma/bioresearch/DSSRprocessedfiles')
    #         os.chdir('/Users/curtisma/bioresearch/DSSRprocessedfiles')
    #         dssrlistoffiles = os.listdir('.')
    #         newdssrlist = filter(re.compile("dssr-").match,dssrlistoffiles)
    #         for g in newdssrlist:
    #             dssrlhs,dssrrhs=g.split("-",1)
    #             dssrstringprep=('%s%s-%s' %(dssrlhs,pdblhs,dssrrhs))
    #             os.system('cp %s %s' %(g,dssrstringprep))

    def dssrcli_test(self):
        #make the directory for the DSSR processed files
        os.system('mkdir ~/bioresearch/compbio/files_wip/DSSRprocessed_baseset')
        #go to the baseset complexes pdb that we need to DSSR on
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/base_complexes_pdb'))
        #get the list of proteins in the baseset
        base_complexes_pdbs = os.listdir('.')
        #for each of the proteins in the baseset
        for base_complex_pdb in base_complexes_pdbs:
            #lets prep the string to feed into DSSR
            folder_stringprep=('~/bioresearch/compbio/pdbfiles/%s' %base_complex_pdb)
            # move into the bin folder where the executable dssr is
            os.chdir(os.path.expanduser('~/bioresearch/compbio/bin'))
            #prepare the json output file name
            baseset_pdb_filename,throwaway=base_complex_pdb.split(".",1)
            #execute dssr with the -json output flag
            os.system('./x3dna-dssr input=%s -o=%s.json -json' %(folder_stringprep,baseset_pdb_filename))
            #move the file to the folder we care about
            os.system('mv %s.json ~/bioresearch/compbio/files_wip/DSSRprocessed_baseset' %baseset_pdb_filename)




    #   print f


# def dssrcli():
#   os.chdir('/Users/curtisma/bioresearch/compbio/pdbfiles')
#   listoffiles = os.listdir('.')
#   os.chdir('/Users/curtisma/bioresearch/')
#   os.system('mkdir DSSRprocessedfiles')
#   for f in listoffiles:
#     stringprep=('~/bioresearch/compbio/pdbfiles/%s' %f)
#     os.chdir('/Users/curtisma/bioresearch/')
#     os.system('./x3dna-dssr input=%s' %stringprep)
#     pdblhs,pdbrhs=f.split(".",1)
#     # print pdblhs
#     os.system('mv dssr* /Users/curtisma/bioresearch/DSSRprocessedfiles')
#     os.chdir('/Users/curtisma/bioresearch/DSSRprocessedfiles')
#     dssrlistoffiles = os.listdir('.')
#     newdssrlist = filter(re.compile("dssr-").match,dssrlistoffiles)
#     for g in newdssrlist:
#       dssrlhs,dssrrhs=g.split("-",1)
#       dssrstringprep=('%s%s-%s' %(dssrlhs,pdblhs,dssrrhs))
#       os.system('cp %s %s' %(g,dssrstringprep))
