import os as os #import the python os library
import csv as csv #import the python csv library
import re #provides regular expression matching

class hbplusclihelper(object):
    def __init__(self):
        self.listoffiles = []

    def foldersetup(self,folder):
        os.chdir('/Users/curtisma/bioresearch/compbio/pdbfiles')
        self.listoffiles = os.listdir('.')
        os.chdir('/Users/curtisma/bioresearch/')
        os.system('mkdir %s'%folder)
        os.chdir('/Users/curtisma/bioresearch/hbplus')

    # run hbplus in hydrogen bond mode
    def hbplushbcli(self):
        self.foldersetup("hbplusprocessedhbfiles")
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
        for f in self.listoffiles:
            stringprep=('~/bioresearch/compbio/pdbfiles/%s' %f)
            os.system('./hbplus %s -N' %stringprep)
            lhs,rhs=f.split(".",1)
            print lhs
            #output file is an .nb2
            os.system('mv %s.nb2 ~/bioresearch/hbplusprocessedvdwfiles' %lhs)
            print f
