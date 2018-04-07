from dssrparserstructure import DssrParserjson
import os as os #import the python os library
import csv as csv #import the python csv library
import re #provides regular expression matching
import json
import glob #file directory patter matcher

class DssrParser(object):

    def dssrprocessedreader(self):
        os.chdir('/Users/curtisma/bioresearch/DSSRprocessedfiles')
        listofdssrprocessedfiles = glob.glob('*.json')
        print listofdssrprocessedfiles
        # listofdssrprocessedfiles = ["pdb1mnb.json","pdb1b2m.json"]
        # listofdssrprocessedfiles = ["pdb1mnb.json","pdb1b2m.json"]
        # listofdssrprocessedfiles = ["pdb2vqe.json"]
        for i in listofdssrprocessedfiles:
            os.chdir('/Users/curtisma/bioresearch/DSSRprocessedfiles')
            print os.getcwd()
            dssrstore=[]
            lhs,rhs=i.split(".",1)
            print lhs
            data = json.load(open(i))
            # dssroutputcategories = ["stems","junction","hairpins","torsions","stacks","splays","pairs","multiplets","helices","bulges","atom2bases"]
            dssroutputcategories = ["helices","iloops","hairpins","bulges"]
            for j in dssroutputcategories:
                print j
                try:
                    # print(data[j])
                    if j == "hairpins":
                        hairpinresult=DssrParserjson().hairpindssrparser(data,j)
                        for line in hairpinresult:
                            dssrstore.append(line)
                    # deprecating stems
                    # elif j == "stems":
                    #     stemresult=DssrParserjson().dssrstemParser(data,j)
                    #     for line in stemresult:
                    #         dssrstore.append(line)
                    elif j == "helices":
                        helixresult=DssrParserjson().dssrhelixParser(data,j)
                        for line in helixresult:
                            dssrstore.append(line)
                    # elif j == "bulges":
                    #     bulgesresult=DssrParserjson().dssrhelixParser(data,j)
                    #     for line in helixresult:
                    #         dssrstore.append(line)
                    # elif j == "iloops":
                    #     helixresult=DssrParserjson().dssrhelixParser(data,j)
                    #     for line in helixresult:
                    #         dssrstore.append(line)
                  # reminder when parsing helix structures must take into the account strand 1 and strand 2 these residue structures are the helices
                except:
                    print "There was an exception likely null"
            # print dssrstore
            self.dssrstorefilewriter(lhs,dssrstore)

    def dssrstorefilewriter(self,proteinname,dssrstore):
      os.chdir("/Users/curtisma/bioresearch")
      os.system('mkdir DSSRprocessedfiles')
      os.chdir("/Users/curtisma/bioresearch/DSSRparsedfiles")
      # os.system('touch dssrparsed.dsr')
      filenamestring="%s.dsr" %proteinname
      os.system("rm %s" %filenamestring )
      file = open(filenamestring,"a")
      print "%s has been created" %filenamestring
      for line in dssrstore:
        file.write("%s\n" %line)
      os.chdir('/Users/curtisma/bioresearch')
