import os as os #import the python os library
import csv as csv #import the python csv library
import re #provides regular expression matching

class HbPlusCliHelper(object):
    def __init__(self):
        self.listoffiles = []

    def foldersetup(self,folder):
        # os.chdir('/Users/curtisma/bioresearch/compbio/pdbfiles')
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/base_complexes_pdb'))
        self.listoffiles = os.listdir('.')
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip'))
        # os.chdir('/Users/curtisma/bioresearch/')
        os.system('mkdir %s'%folder)
        os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/hbplus'))
        # os.chdir('/Users/curtisma/bioresearch/hbplus')

    # run hbplus in hydrogen bond mode
    def hbplushbcli(self):
        self.foldersetup("hbplus_processed_hb_files_baseset")
        print self.listoffiles
        for f in self.listoffiles:
            stringprep=('~/bioresearch/compbio/files_wip/base_complexes_pdb/%s' %f)
            os.system('./hbplus %s' %stringprep)
            lhs,rhs=f.split(".",1)
            print lhs
            #output file is an .hb2
            os.system('mv %s.hb2 ~/bioresearch/compbio/files_wip/hbplus_processed_hb_files_baseset' %lhs)
            print f

    # run hbplus in van der Waal mode
    def hbplusvdwcli(self):
        # self.foldersetup("hbplusprocessedvdwfiles_baseset")
        self.foldersetup("hbplus_processed_vdw_files_baseset")
        for f in self.listoffiles:
            # stringprep=('~/bioresearch/compbio/pdbfiles/%s' %f)
            stringprep=('~/bioresearch/compbio/files_wip/base_complexes_pdb/%s' %f)
            os.system('./hbplus %s -N' %stringprep)
            lhs,rhs=f.split(".",1)
            print lhs
            #output file is an .nb2
            os.system('mv %s.nb2 ~/bioresearch/compbio/files_wip/hbplus_processed_vdw_files_baseset' %lhs)
            print f
