from dssrparserstructure import DssrParserjson
import os as os #import the python os library
import csv as csv #import the python csv library
import re #provides regular expression matching
import json
import glob #file directory patter matcher

class DssrParserTestSet(object):

    def dssrprocessedreader(self):
        os.system('mkdir ~/bioresearch/compbio/files_wip/DSSRparsed_testset')
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/DSSRprocessed_testset'))
        list_of_testset_proteins = os.listdir('.')
        for testset_protein in list_of_testset_proteins:
            os.system('mkdir ~/bioresearch/compbio/files_wip/DSSRparsed_testset/%s' %testset_protein)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/DSSRprocessed_testset/%s'%testset_protein))
            list_of_testset_protein_poses = os.listdir('.')
            for i in list_of_testset_protein_poses:
                os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/DSSRprocessed_testset/%s'%testset_protein))
                print os.getcwd()
                dssrstore=[]
                lhs,rhs=i.split(".",1)
                print lhs
                data = json.load(open(i))
                # dssroutputcategories = ["stems","junction","hairpins","torsions","stacks","splays","pairs","multiplets","helices","bulges","atom2bases"]
                dssroutputcategories = ["helices","iloops","hairpins","bulges","junctions"]
                for j in dssroutputcategories:
                    print j
                    try:
                        # print(data[j])
                        if j == "hairpins":
                            hairpinresult=DssrParserjson().dssrhairpinParser(data,j)
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
                        elif j == "bulges":
                            bulgesresult=DssrParserjson().dssrbulgeParser(data,j)
                            for line in bulgesresult:
                                dssrstore.append(line)
                        elif j == "iloops":
                            iloopsresult=DssrParserjson().dssriloopsParser(data,j)
                            for line in iloopsresult:
                                dssrstore.append(line)
                        elif j == "junctions":
                            junctionsresult=DssrParserjson().dssrjunctionsParser(data,j)
                            for line in junctionsresult:
                                dssrstore.append(line)
                      # reminder when parsing helix structures must take into the account strand 1 and strand 2 these residue structures are the helices
                    except Exception as e:
                        print "There was an exception likely null %s" %e
                # print dssrstore
                self.dssrstorefilewriter(testset_protein,lhs,dssrstore)

    def dssrstorefilewriter(self,protein_name,pose_name,dssrstore):
      os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/DSSRparsed_testset/%s'%protein_name))
      filenamestring="%s.dsr" %pose_name
      os.system("rm %s" %filenamestring)
      file = open(filenamestring,"a")
      print "%s has been created" %filenamestring
      print os.getcwd()
      for line in dssrstore:
        file.write("%s\n" %line)
