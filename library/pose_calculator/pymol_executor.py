from Bio import *                                           #use the BioPython Python library
from Bio.PDB import *                                       #more specifically import the BioPython PDB library
import os
import numpy as np
import math

def pymol_cli_executor():
    os.system("mkdir pymol_rmsd_calculated_files")
    os.chdir(os.path.expanduser('~/bioresearch/combio/bin'))
    os.system("./pymol <input> <output>")

    os.chdir(os.path.expanduser('~/bioresearch/compbio'))
    os.chdir(os.path.expanduser('~/bioresearch/protein_seperated_pdbfiles'))
