from dssrparserstructure import DssrParserjson
import os as os #import the python os library
import csv as csv #import the python csv library
import re #provides regular expression matching
import json
import glob #file directory patter matcher

class DssrParser(object):

    def dssrprocessedreader(self):
        # make the directory for the parsed files
        os.system('mkdir ~/bioresearch/compbio/files_wip/DSSRparsed_baseset')
        # move into the directory where the procssed DSSR files are to begin parsing!
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/DSSRprocessed_baseset'))
        dssr_processed_files = glob.glob('*.json')
        # for each file in the list
        for dssr_process_file in dssr_processed_files:
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/DSSRprocessed_baseset'))
            dssrstore=[]
            baseset_pdb_filename,throwaway=dssr_process_file.split(".",1)
            print baseset_pdb_filename
            data = json.load(open(dssr_process_file))
            # dssroutputcategories = ["stems","junction","hairpins","torsions","stacks","splays","pairs","multiplets","helices","bulges","atom2bases"]
            dssroutputcategories = ["helices","iloops","hairpins","bulges","junctions"]
            for rna_secondary_structure_type in dssroutputcategories:
                print rna_secondary_structure_type
                try:
                    # print(data[j])
                    if rna_secondary_structure_type == "hairpins":
                        hairpinresult=DssrParserjson().dssrhairpinParser(data,rna_secondary_structure_type)
                        for line in hairpinresult:
                            dssrstore.append(line)
                    # deprecating stems
                    # elif j == "stems":
                    #     stemresult=DssrParserjson().dssrstemParser(data,j)
                    #     for line in stemresult:
                    #         dssrstore.append(line)
                    elif rna_secondary_structure_type == "helices":
                        helixresult=DssrParserjson().dssrhelixParser(data,rna_secondary_structure_type)
                        for line in helixresult:
                            dssrstore.append(line)
                    elif rna_secondary_structure_type == "bulges":
                        bulgesresult=DssrParserjson().dssrbulgeParser(data,rna_secondary_structure_type)
                        for line in bulgesresult:
                            dssrstore.append(line)
                    elif rna_secondary_structure_type == "iloops":
                        iloopsresult=DssrParserjson().dssriloopsParser(data,rna_secondary_structure_type)
                        for line in iloopsresult:
                            dssrstore.append(line)
                    elif rna_secondary_structure_type == "junctions":
                        junctions_result=DssrParserjson().dssrjunctionsParser(data,rna_secondary_structure_type)
                        for line in junctions_result:
                            dssrstore.append(line)
                  # parsing helix structures must take into the account strand 1 and strand 2 these residue structures are the helices
                except Exception as e:
                    print "There was an exception likely null %s" %e
            self.dssrstorefilewriter(baseset_pdb_filename,dssrstore)

    def dssrstorefilewriter(self,proteinname,dssrstore):
        # move into the directory to write the file.
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/DSSRparsed_baseset'))
        filenamestring="%s.dsr" %proteinname
        os.system("rm %s" %filenamestring )
        file = open(filenamestring,"a")
        print "%s has been created" %filenamestring
        for line in dssrstore:
            file.write("%s\n" %line)
            os.chdir(os.path.expanduser('~/bioresearch'))
