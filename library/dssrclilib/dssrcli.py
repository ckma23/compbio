import os as os #import the python os library
import csv as csv #import the python csv library
import re #provides regular expression matching
import json

class dssrclihelper(object):
    def dssrcli(self):
        os.chdir('/Users/curtisma/bioresearch/compbio/pdbfiles')
        listoffiles = os.listdir('.')
        os.chdir('/Users/curtisma/bioresearch/')
        os.system('mkdir DSSRprocessedfiles')
        for f in listoffiles:
            stringprep=('~/bioresearch/compbio/pdbfiles/%s' %f)
            os.chdir('/Users/curtisma/bioresearch/')
            lhs,rhs=f.split(".",1)
            os.system('./x3dna-dssr input=%s -o=%s.json -json' %(stringprep,lhs))
            # os.system('./x3dna-dssr input=%s' %stringprep)
            pdblhs,pdbrhs=f.split(".",1)
            # print pdblhs
            os.system('mv %s.json /Users/curtisma/bioresearch/DSSRprocessedfiles' %lhs)
            os.system('mv dssr* /Users/curtisma/bioresearch/DSSRprocessedfiles')
            os.chdir('/Users/curtisma/bioresearch/DSSRprocessedfiles')
            dssrlistoffiles = os.listdir('.')
            newdssrlist = filter(re.compile("dssr-").match,dssrlistoffiles)
            for g in newdssrlist:
                dssrlhs,dssrrhs=g.split("-",1)
                dssrstringprep=('%s%s-%s' %(dssrlhs,pdblhs,dssrrhs))
                os.system('cp %s %s' %(g,dssrstringprep))
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
