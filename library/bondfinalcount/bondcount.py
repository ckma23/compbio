from Bio import *                                           #use the BioPython Python library
from Bio.PDB import *                                       #more specifically import the BioPython PDB library
import os
import numpy as np
import math

class Bondcounter(object):
    def bondcounter(self):
        # returns a prepared hash of based on hydrogen bonding or Vander Waals, which category, which nucleotide_base, which amino_acid
        nucleotide_base_amino_acid = {}
        category=["CAT_1","CAT_2","CAT_3","CAT_4","CAT_5","CAT_6","CAT_7","CAT_8","CAT_9"]
        nucleotide_base =["A","C","U","G"]
        amino_acid=["ARG","ALA","ASN","ASP","GLN","GLU","GLY","CYS","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
        vdw_or_hb=["hb","vdw"]
        for vdorhb in vdw_or_hb:
            nucleotide_base_amino_acid[vdorhb]={}
            for cat in category:
                nucleotide_base_amino_acid[vdorhb][cat]={}
                for nb in nucleotide_base:
                    nucleotide_base_amino_acid[vdorhb][cat][nb] = {}
                    for aa in amino_acid:
                        nucleotide_base_amino_acid[vdorhb][cat][nb][aa]={"count":0,"statistical_potential":0}
        # print nucleotide_base_amino_acid.keys()
        # print nucleotide_base_amino_acid["hb"].keys()
        # print nucleotide_base_amino_acid["hb"]["CAT_1"].keys()
        # print nucleotide_base_amino_acid["hb"]["CAT_1"]["A"].keys()
        return nucleotide_base_amino_acid

    def file_opener(self):
        final_bond_counted_store = Bondcounter().bondcounter()
        os.chdir(os.path.expanduser('~/bioresearch/bondcategorized'))
        listofcategorizedfiles = os.listdir('.')
        for categorized_file in listofcategorizedfiles:
            categorized_file_open = open(categorized_file,"r")
            for line in categorized_file_open:
                category = line[0:5]
                hborvdw = line[34:37].strip(' ')
                nucleotide_base = line[6]
                amino_acid = line[25:28]
                # we should switch this to spaces..
                # final_bond_counted_store[hborvdw][category][nucleotide_base][amino_acid]["count"]+=1
                try:
                    final_bond_counted_store[hborvdw][category][nucleotide_base][amino_acid]["count"]+=1
                except:
                    print "Clean up string: %s" %line
        print final_bond_counted_store
        return final_bond_counted_store
