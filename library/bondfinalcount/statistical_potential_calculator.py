from Bio import *                                           #use the BioPython Python library
from Bio.PDB import *                                       #more specifically import the BioPython PDB library
import os
import numpy as np
import math
import json
from bondcount import Bondcounter
from structurecount import StructureCounter

from library.energy_formation.energy_balancer import Energyknockdown as Energyknockdown

class StatisticalPotential(object):
    def statistical_potential(self):
        #lets get all the bonds counted!
        # final_bond_counted_store = Bondcounter().bondcounter()
        bond_counted_hash = Bondcounter().file_opener()
        # print bond_counted_hash
        #lets get all the nucleotide_base and amino_acid counted!
        # structure_hash = StructureCounter().structure_hasher
        structure_counted_hash = StructureCounter().structure_counter()

        nb_learning_set_sum = 0
        aa_learning_set_sum = 0
        nucleotide_base = ["A","C","U","G"]
        amino_acid = ["ARG","ALA","ASN","ASP","GLN","GLU","GLY","CYS","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
        for nb_learningset in nucleotide_base:
            nb_learning_set_sum += structure_counted_hash[nb_learningset]
        for aa_learningset in amino_acid:
            aa_learning_set_sum += structure_counted_hash[aa_learningset]
        # print structure_counted_hash
        for vdworhb in bond_counted_hash.keys():
            for cat in sorted(bond_counted_hash[vdworhb].keys()):
                for nb in bond_counted_hash[vdworhb][cat].keys():
                    for aa in bond_counted_hash[vdworhb][cat][nb].keys():
                        count_all_pairs_cat = 0
                        count_all_nb_type_cat_type = 0
                        count_all_aa_type_cat_type = 0
                        for nucleotide_base in bond_counted_hash[vdworhb][cat].keys():
                            for amino_acid in bond_counted_hash[vdworhb][cat][nb].keys():
                                count_all_pairs_cat += bond_counted_hash[vdworhb][cat][nucleotide_base][amino_acid]["count"]

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
                            # as 10/8/18 rebuilding the potentials
                            # first_denom=round(float(count_all_nb_type_cat_type)/structure_counted_hash[nb],3)
                            # second_denom=round(float(count_all_aa_type_cat_type)/structure_counted_hash[aa],3)
                            first_denom=round(float(structure_counted_hash[nb])/nb_learning_set_sum,3)
                            second_denom=round(float(structure_counted_hash[aa])/aa_learning_set_sum,3)
                            propensity = numerator/(first_denom*second_denom)
                            # print "start"
                            # print numerator
                            # print first_denom
                            # print second_denom
                        if vdworhb == "hb":
                            RT=.59
                        # DEVELOP DIFFERENT POTENTIALS FOR HB OR VDW
                        elif vdworhb == "vdw":
                            RT=.59
                        # print bond_counted_hash[vdworhb][cat][nb][aa]["count"]
                        # print propensity
                        # log in python is actually the natural log
                        try:
                            statistical_potential=-1*RT*math.log(propensity)
                        except Exception as e:
                            # print "error in calculating statistical_potential defaulting to 0.001"
                            propensity = 0.001
                            statistical_potential=-1*RT*math.log(propensity)

                        bond_counted_hash[vdworhb][cat][nb][aa]["statistical_potential"] = statistical_potential
                        line = "%s,%s,%s,%s,%s" %(vdworhb,cat,nb,aa,bond_counted_hash[vdworhb][cat][nb][aa]["statistical_potential"])
                        print line
        stat_potential_file_path = os.path.join(os.path.expanduser('~/bioresearch/compbio/metadata_folder/'),"statistical_potential.json")
        with open (stat_potential_file_path ,'w') as stat_pot_file:
            json.dump(bond_counted_hash,stat_pot_file)
        # print bond_counted_hash
        # balance_stat_potential = self.energy_balancer(bond_counted_hash)
        print "All the AA in the learning set"
        print aa_learning_set_sum
        print "All the NB in the learning set"
        print nb_learning_set_sum
        print "Specific potentials vdw, cat_9,A "
        print bond_counted_hash["vdw"]["CAT_9"]["A"]
        return bond_counted_hash

    # this is meant to rebalance the statistical_potentials
    def energy_balancer(self,bond_counted_hash):
        balance_hash = energyknockdown().energyknockdownretriever()
        for vdworhb in bond_counted_hash.keys():
            for cat in sorted(bond_counted_hash[vdworhb].keys()):
                for nb in bond_counted_hash[vdworhb][cat].keys():
                    for aa in bond_counted_hash[vdworhb][cat][nb].keys():
                        bond_counted_hash[vdworhb][cat][nb][aa]["statistical_potential"] = bond_counted_hash[vdworhb][cat][nb][aa]["statistical_potential"] * balance_hash[vdworhb][cat]

        print balance_hash
        return bond_counted_hash

    def statistical_poential(self,propensity):
        RT=.59
        statistical_potential=-1*RT*math.log(propensity)
        return statistical_potential

# def all_pair_counter():
    # first denom = bond_counted_hash["hb"]["CAT_1"]["A"]["GUA"]/
    #numerator = count of arg-guanine in Cat 1/count of all pairs in Cat 1
    #first demoninator = count of arg in Cat 1/count of all Arg in learning set
    #second denominator = count of guanine in Cat 1/ count of all guanine in learning set
