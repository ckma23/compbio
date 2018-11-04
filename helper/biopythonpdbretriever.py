import os as os
import csv as csv
import Bio
from Bio.PDB import PDBList

class FileRetriever (object):
    def fileretrieving(self,string,foldername):
        print "Downloading the %s pdb files now!" %string
        filedir= "%s" %(os.path.expanduser('~/bioresearch/compbio/%s')) %foldername
        pdbl = PDBList()
        # note the file_format can be taken in as http://biopython.org/DIST/docs/api/Bio.PDB.PDBList%27.PDBList-class.html cif,pdb,etc
        pdbl.retrieve_pdb_file(string,pdir=filedir,file_format="pdb")

    def structure_of_interest_return(self):
        #this method opens the structure_of_interest.csv file and returns an object array containing the objects of interest to the caller
        with open('structures_of_interest.csv','rb') as csvfile:
            csvstore = csv.reader(csvfile,delimiter = ',')
            structure_of_interest_array = []
            for structure_id in csvstore:
                structure_of_interest_array.append(''.join(structure_id))
            return structure_of_interest_array
