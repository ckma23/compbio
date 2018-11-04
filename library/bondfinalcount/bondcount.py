from Bio import *                                           #use the BioPython Python library
from Bio.PDB import *                                       #more specifically import the BioPython PDB library
import os
import numpy as np
import math
import csv

class Bondcounter(object):
    def bondcounter(self):
        # returns a prepared hash of based on hydrogen bonding or Vander Waals, which category, which nucleotide_base, which amino_acid
        nucleotide_base_amino_acid = {}
        category=["CAT_1","CAT_2","CAT_3","CAT_4","CAT_5","CAT_6","CAT_7","CAT_8","CAT_9"]
        # category=["CAT_1","CAT_2","CAT_3","CAT_4","CAT_5","CAT_6","CAT_7","CAT_8","CAT_9","CAT_10","CAT_11"]
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
        # os.chdir(os.path.expanduser('~/bioresearch/bondcategorized'))
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/bondcategorized_baseset'))
        listofcategorizedfiles = os.listdir('.')
        for categorized_file in listofcategorizedfiles:
            # categorized_file_open = open(categorized_file,"r")

            with open(categorized_file,'rb') as categorized_file_meta_open:
                categorized_file_open = csv.reader(categorized_file_meta_open)

                for line in categorized_file_open:
                    #the output from a comma spearated bondcategorized_file is expected to be
                    #CAT_4,C0023,B,N4,A0010-,ARG,CD,vdw,helices,B,C0023,1,6,x,nt2,cW-W,WC
                    #first value is the category line[0]
                    #eight value is the hborvdw line[7]
                    #second value is the nucleotidebase but want the first letter line[1][0]
                    #sixth value is the amino acide line[5]
                    category = line[0]
                    hborvdw = line[7]
                    nucleotide_base = line[1][0]
                    amino_acid = line[5]

                    # category = line[0:5]
                    # hborvdw = line[34:37].strip(' ')
                    # nucleotide_base = line[6]
                    # amino_acid = line[25:28]
                    try:
                        final_bond_counted_store[hborvdw][category][nucleotide_base][amino_acid]["count"]+=1
                    except Exception as e:
                        print "Clean up string: %s" %line
                        print e

        number_by_bondtype = {"vdw":0,"hb":0}
        for vdworhb in final_bond_counted_store.keys():

            number_by_CAT = {"CAT_1":0,"CAT_2":0,"CAT_3":0,"CAT_4":0,"CAT_5":0,"CAT_6":0,"CAT_7":0,"CAT_8":0,"CAT_9":0}
            # number_by_CAT = {"CAT_1":0,"CAT_2":0,"CAT_3":0,"CAT_4":0,"CAT_5":0,"CAT_6":0,"CAT_7":0,"CAT_8":0,"CAT_9":0,"CAT_10":0,"CAT_11":0}
            number_by_NB_in_all_CAT = {"A":0,"U":0,"C":0,"G":0}
            number_by_NB_to_AA_in_all_CAT = {"A":{"ARG":0,"ALA":0,"ASN":0,"ASP":0,"GLN":0,"GLU":0,"GLY":0,"CYS":0,"HIS":0,"ILE":0,"LEU":0,"LYS":0,"MET":0,"PHE":0,"PRO":0,"SER":0,"THR":0,"TRP":0,"TYR":0,"VAL":0},"U":{"ARG":0,"ALA":0,"ASN":0,"ASP":0,"GLN":0,"GLU":0,"GLY":0,"CYS":0,"HIS":0,"ILE":0,"LEU":0,"LYS":0,"MET":0,"PHE":0,"PRO":0,"SER":0,"THR":0,"TRP":0,"TYR":0,"VAL":0},"C":{"ARG":0,"ALA":0,"ASN":0,"ASP":0,"GLN":0,"GLU":0,"GLY":0,"CYS":0,"HIS":0,"ILE":0,"LEU":0,"LYS":0,"MET":0,"PHE":0,"PRO":0,"SER":0,"THR":0,"TRP":0,"TYR":0,"VAL":0},"G":{"ARG":0,"ALA":0,"ASN":0,"ASP":0,"GLN":0,"GLU":0,"GLY":0,"CYS":0,"HIS":0,"ILE":0,"LEU":0,"LYS":0,"MET":0,"PHE":0,"PRO":0,"SER":0,"THR":0,"TRP":0,"TYR":0,"VAL":0}}
            for cat in final_bond_counted_store[vdworhb].keys():
                for nb in final_bond_counted_store[vdworhb][cat].keys():
                    for aa in final_bond_counted_store[vdworhb][cat][nb].keys():
                        number_by_bondtype[vdworhb]+=final_bond_counted_store[vdworhb][cat][nb][aa]["count"]
                        number_by_CAT[cat]+= final_bond_counted_store[vdworhb][cat][nb][aa]["count"]
                        number_by_NB_in_all_CAT[nb] += final_bond_counted_store[vdworhb][cat][nb][aa]["count"]
                        number_by_NB_to_AA_in_all_CAT[nb][aa]+=final_bond_counted_store[vdworhb][cat][nb][aa]["count"]
            print "\nBY NB %s" %vdworhb
            print number_by_NB_in_all_CAT
            print "\nBY NB TO AA %s" %vdworhb
            print number_by_NB_to_AA_in_all_CAT
            print "\nBY NB TO AA print friendly %s" %vdworhb
            for nb in sorted(number_by_NB_to_AA_in_all_CAT.keys()):
                for aa in sorted(number_by_NB_to_AA_in_all_CAT[nb].keys()):
                    print "%s,%s,%s" %(nb,aa,number_by_NB_to_AA_in_all_CAT[nb][aa])
            print "\nBY CATEGORY %s" %vdworhb
            print number_by_CAT
            for cat in sorted(number_by_CAT.keys()):
                print "%s,%s" %(cat,number_by_CAT[cat])

            print "\nBY NB TO AA BY CAT print friendly %s" %vdworhb
            # category=["CAT_1","CAT_2","CAT_3","CAT_4","CAT_5","CAT_6","CAT_7","CAT_8","CAT_9","CAT_10","CAT_11"]
            category=["CAT_1","CAT_2","CAT_3","CAT_4","CAT_5","CAT_6","CAT_7","CAT_8","CAT_9"]
            nucleotide_base =["A","C","U","G"]
            amino_acid=["ARG","ALA","ASN","ASP","GLN","GLU","GLY","CYS","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
            vdw_or_hb=["hb","vdw"]
            for cat in category:
                for nb in nucleotide_base:
                    for aa in amino_acid:
                        print "%s,%s,%s,%s" %(cat,nb,aa,final_bond_counted_store[vdworhb][cat][nb][aa]["count"])
        print "\nBY BONDTYPE %s" %vdworhb
        print number_by_bondtype

        # os.system("rm results")
        # file = open("results","a")

        # print final_bond_counted_store
        return final_bond_counted_store
