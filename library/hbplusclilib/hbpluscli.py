import os as os #import the python os library
import csv as csv #import the python csv library
import re #provides regular expression matching

class hbplusclihelper(object):
    def __init__(self):
        self.listoffiles = []
    # run hbplus in hydrogen bond mode
    def foldersetup(self,folder):
        os.chdir('/Users/curtisma/bioresearch/compbio/pdbfiles')
        self.listoffiles = os.listdir('.')
        os.chdir('/Users/curtisma/bioresearch/')
        os.system('mkdir %s'%folder) 
        os.chdir('/Users/curtisma/bioresearch/hbplus')

    def hbplushbcli(self):
        self.foldersetup("hbplusprocessedhbfiles")
        # os.chdir('/Users/curtisma/bioresearch/compbio/pdbfiles')
        # listoffiles = os.listdir('.')
        # os.chdir('/Users/curtisma/bioresearch/')
        # os.system('mkdir hbplusprocessedhbfiles')
        # os.chdir('/Users/curtisma/bioresearch/hbplus')
        print self.listoffiles
        for f in self.listoffiles:
            stringprep=('~/bioresearch/compbio/pdbfiles/%s' %f)
            os.system('./hbplus %s' %stringprep)
            lhs,rhs=f.split(".",1)
            print lhs
            #output file is an .hb2
            os.system('mv %s.hb2 ~/bioresearch/hbplusprocessedhbfiles' %lhs)
            print f

    # run hbplus in van der Waal mode
    def hbplusvdwcli(self):
        self.foldersetup("hbplusprocessedvdwfiles")
        # os.chdir('/Users/curtisma/bioresearch/compbio/pdbfiles')
        # listoffiles = os.listdir('.')
        # os.chdir('/Users/curtisma/bioresearch/')
        # os.system('mkdir hbplusprocessedvdwfiles')
        # os.chdir('/Users/curtisma/bioresearch/hbplus')
        for f in self.listoffiles:
            stringprep=('~/bioresearch/compbio/pdbfiles/%s' %f)
            os.system('./hbplus %s -N' %stringprep)
            lhs,rhs=f.split(".",1)
            print lhs
            #output file is an .nb2
            os.system('mv %s.nb2 ~/bioresearch/hbplusprocessedvdwfiles' %lhs)
            print f

# hbpluscli().hbplushbcli()