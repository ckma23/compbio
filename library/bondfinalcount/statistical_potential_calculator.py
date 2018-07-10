from Bio import *                                           #use the BioPython Python library
from Bio.PDB import *                                       #more specifically import the BioPython PDB library
import os
import numpy as np
import math
from bondcount import Bondcounter
from structurecount import StructureCounter

class StatisticalPotential(object):
    def statistical_potential(self):
        #lets get all the bonds counted!
        # final_bond_counted_store = Bondcounter().bondcounter()
        bond_counted_hash = Bondcounter().file_opener()
        print bond_counted_hash
        #lets get all the nucleotide_base and amino_acid counted!
        # structure_hash = StructureCounter().structure_hasher

        structure_counted_hash = StructureCounter().structure_counter()
        print structure_counted_hash
        for vdworhb in bond_counted_hash.keys():
            for cat in bond_counted_hash[vdworhb].keys():
                for nb in bond_counted_hash[vdworhb][cat].keys():
                    for aa in bond_counted_hash[vdworhb][cat][nb].keys():
                        count_all_pairs_cat = 0
                        count_all_nb_type_cat_type = 0
                        count_all_aa_type_cat_type = 0
                        for nucleotide_base in bond_counted_hash[vdworhb][cat].keys():
                            for amino_acid in bond_counted_hash[vdworhb][cat][nb].keys():
                                count_all_pairs_cat += bond_counted_hash[vdworhb][cat][nucleotide_base][amino_acid]["count"]

                        # for whatever nucleotide_base we are on, count all the nucleotide base in CAT_1
                        for amino_acid_cat in bond_counted_hash[vdworhb][cat][nb].keys():
                            count_all_nb_type_cat_type += bond_counted_hash[vdworhb][cat][nb][amino_acid_cat]["count"]

                        # for whichever amino_acid we are one, count all the amino_cide base in CAT_1
                        for nucleotide_base_cat in bond_counted_hash[vdworhb][cat].keys():
                            count_all_aa_type_cat_type += bond_counted_hash[vdworhb][cat][nucleotide_base_cat][aa]["count"]
                        # one of the numbers have to be a float!
                        # print bond_counted_hash[vdworhb][cat][nb][aa]["count"]
                        # print type(bond_counted_hash[vdworhb][cat][nb][aa]["count"])
                        # print count_all_nb_type_cat_type
                        # print count_all_aa_type_cat_type

                        # LETS MAKE SURE  we are dividing each of the hb by hb and the vdw by vdw
                        # and not hb/all and vdw/all
                        if bond_counted_hash[vdworhb][cat][nb][aa]["count"] == 0:
                            propensity = 0.001 #DEFAULT POTENTIAL VALUE
                        else:
                            numerator=round(float(bond_counted_hash[vdworhb][cat][nb][aa]["count"])/count_all_pairs_cat,3)
                            first_denom=round(float(count_all_nb_type_cat_type)/structure_counted_hash[nb],3)
                            second_denom=round(float(count_all_aa_type_cat_type)/structure_counted_hash[aa],3)
                            propensity = numerator/(first_denom+second_denom)
                            # print "start"
                            # print numerator
                            # print first_denom
                            # print second_denom
                        RT=.59
                        # print bond_counted_hash[vdworhb][cat][nb][aa]["count"]
                        # print propensity
                        try:
                            statistical_potential=-1*RT*math.log(propensity)
                        except Exception as e:
                            print "error in calculating statistical_potential defaulting to 0.001"
                            propensity = 0.001
                            statistical_potential=-1*RT*math.log(propensity)

                        bond_counted_hash[vdworhb][cat][nb][aa]["statistical_potential"] = statistical_potential
                        print "%s %s %s %s %s" %(vdworhb,cat,nb,aa,bond_counted_hash[vdworhb][cat][nb][aa]["statistical_potential"])
        print bond_counted_hash
        return bond_counted_hash

    def statistical_poential(self,propensity):
        RT=.59
        statistical_potential=-1*RT*math.log(propensity)
        return statistical_poential

# def all_pair_counter():
    # first denom = bond_counted_hash["hb"]["CAT_1"]["A"]["GUA"]/
    #numerator = count of arg-guanine in Cat 1/count of all pairs in Cat 1
    #first demoninator = count of arg in Cat 1/count of all Arg in learning set
    #second denominator = count of guanine in Cat 1/ count of all guanine in learning set


# statistical_potential()
