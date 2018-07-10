from Bio import *                                           #use the BioPython Python library
from Bio.PDB import *                                       #more specifically import the BioPython PDB library
import os
import numpy as np
import math

class StructureCounter(object):
    def structure_hasher(self):
        # returns a prepared hash of based on which nucleotide_base, which amino_acid
        nucleotide_base_amino_acid = {}
        nucleotide_base =["A","C","U","G"]
        amino_acid=["ARG","ALA","ASN","ASP","GLN","GLU","GLY","CYS","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
        for nb in nucleotide_base:
            nucleotide_base_amino_acid[nb]=0
        for aa in amino_acid:
            nucleotide_base_amino_acid[aa]=0
        return nucleotide_base_amino_acid

    def structure_counter(self):
        structure_store = StructureCounter().structure_hasher()
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/base_complexes_pdb'))
        # os.chdir(os.path.expanduser('~/bioresearch/compbio/pdbfiles'))
        listofcategorizedfiles = os.listdir('.')
        for pdbfile in listofcategorizedfiles:
            lhs,rhs = pdbfile.split('.')
            pdb_name = lhs.strip('pdb')
            parser = PDBParser()
            structure = parser.get_structure(pdb_name,pdbfile)
            for model in structure:
                # print model
                for chain in model:
                    # print chain
                    for residue in chain:
                        try:
                            structure_store[residue.get_resname().strip()]+=1
                        except:
                            print "clean up string: %s" %residue.get_resname().strip()
        print structure_store
        return structure_store

# StructureCounter().structure_counter()
