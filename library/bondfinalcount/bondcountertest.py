from Bio import *                                           #use the BioPython Python library
from Bio.PDB import *                                       #more specifically import the BioPython PDB library
import os
import numpy as np

class Bondcounter(object):
    #I should probably make this into a hashtable instead....
    # pro bably want to assign a class variable and let all methods access this.
    def bondcounter(self):
        os.chdir(os.path.expanduser('~/bioresearch/bondcategorized'))
        listofcategorizedfiles = os.listdir('.')
        category_storehb=[0,0,0,0,0,0,0,0,0]
        category_storevdw=[0,0,0,0,0,0,0,0,0]
    ###NUCLEOTIDE BASE # in order of A C U G
        #categorized = [[[0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],[[0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],[[0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],[[0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]]
        #categorized[0][0]=A categorized[1][0]=C categorized[1][0] = G  categorized[1][0] = U
        category_1_nucleotidebase=np.array([0,0,0,0])
        category_2_nucleotidebase=np.array([0,0,0,0])
        category_3_nucleotidebase=np.array([0,0,0,0])
        category_4_nucleotidebase=np.array([0,0,0,0])
        category_5_nucleotidebase=np.array([0,0,0,0])
        category_6_nucleotidebase=np.array([0,0,0,0])
        category_7_nucleotidebase=np.array([0,0,0,0])
        category_8_nucleotidebase=np.array([0,0,0,0])
        category_9_nucleotidebase=np.array([0,0,0,0])
    ###AMINO ACID # in order of ["ARG","ALA","ASN","ASP","GLN","GLU","GLY","CYS","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
        category_1_aminoacid=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        category_2_aminoacid=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        category_3_aminoacid=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        category_4_aminoacid=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        category_5_aminoacid=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        category_6_aminoacid=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        category_7_aminoacid=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        category_8_aminoacid=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        category_9_aminoacid=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        # categorized_1_vdw=([[0,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],[0,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],[0,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],[0,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]])
        # HB
        categorized_1_hb_nb=np.array([0,0,0,0])
        categorized_2_hb_nb=np.array([0,0,0,0])
        categorized_3_hb_nb=np.array([0,0,0,0])
        categorized_4_hb_nb=np.array([0,0,0,0])
        categorized_5_hb_nb=np.array([0,0,0,0])
        categorized_6_hb_nb=np.array([0,0,0,0])
        categorized_7_hb_nb=np.array([0,0,0,0])
        categorized_8_hb_nb=np.array([0,0,0,0])
        categorized_9_hb_nb=np.array([0,0,0,0])

        categorized_1_hb_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_1_hb_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_1_hb_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_1_hb_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_2_hb_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_2_hb_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_2_hb_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_2_hb_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_3_hb_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_3_hb_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_3_hb_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_3_hb_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_4_hb_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_4_hb_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_4_hb_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_4_hb_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_5_hb_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_5_hb_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_5_hb_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_5_hb_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_6_hb_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_6_hb_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_6_hb_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_6_hb_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_7_hb_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_7_hb_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_7_hb_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_7_hb_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_8_hb_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_8_hb_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_8_hb_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_8_hb_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_9_hb_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_9_hb_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_9_hb_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_9_hb_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])


        # categorized_1_vdw=([[0,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],[0,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],[0,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],[0,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]])
        # VDW
        categorized_1_vdw_nb=np.array([0,0,0,0])
        categorized_2_vdw_nb=np.array([0,0,0,0])
        categorized_3_vdw_nb=np.array([0,0,0,0])
        categorized_4_vdw_nb=np.array([0,0,0,0])
        categorized_5_vdw_nb=np.array([0,0,0,0])
        categorized_6_vdw_nb=np.array([0,0,0,0])
        categorized_7_vdw_nb=np.array([0,0,0,0])
        categorized_8_vdw_nb=np.array([0,0,0,0])
        categorized_9_vdw_nb=np.array([0,0,0,0])

        categorized_1_vdw_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_1_vdw_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_1_vdw_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_1_vdw_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_2_vdw_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_2_vdw_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_2_vdw_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_2_vdw_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_3_vdw_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_3_vdw_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_3_vdw_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_3_vdw_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_4_vdw_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_4_vdw_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_4_vdw_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_4_vdw_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_5_vdw_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_5_vdw_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_5_vdw_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_5_vdw_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_6_vdw_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_6_vdw_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_6_vdw_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_6_vdw_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_7_vdw_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_7_vdw_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_7_vdw_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_7_vdw_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_8_vdw_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_8_vdw_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_8_vdw_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_8_vdw_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        categorized_9_vdw_nbA_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_9_vdw_nbC_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_9_vdw_nbG_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        categorized_9_vdw_nbU_aa=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

        for categorized_file in listofcategorizedfiles:
            # with open(categorized_file) as csvfile:
                # csvstore = csv.reader(csvfile, delimiter = '')
                # for line in csvstore:
            categorized_file_open = open(categorized_file,"r")
            for line in categorized_file_open:
                i=1
                #check each line what CAT the CAT is line[0]
                #check each line is hb or vdw is line[8]
                cat = line[0:5]
                # this is returning "hb " the space is not matching up
                hborvdw = line[34:37].strip(' ')
                nucleotide_base = line[6]
                amino_acid = line[25:28]

                if hborvdw == "hb":
                    if cat == "CAT_1":
                        # category_storehb[0]+=1
                        #
                        # categorized_1_hb_nb[0]+=1 if nucleotide_base == "A" else categorized_1_hb_nb[0]
                        # categorized_1_hb_nb[1]+=1 if nucleotide_base == "C" else categorized_1_hb_nb[1]
                        # categorized_1_hb_nb[2]+=1 if nucleotide_base == "G" else categorized_1_hb_nb[2]
                        # categorized_1_hb_nb[3]+=1 if nucleotide_base == "U" else categorized_1_hb_nb[3]

                        # categorized_1_hb_nb = categorized_1_hb_nb + self.nucleotide_basecounter(nucleotide_base)

                        categorized_1_hb_nbA_aa = categorized_1_hb_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_1_hb_nbA_aa
                        categorized_1_hb_nbC_aa = categorized_1_hb_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_1_hb_nbC_aa
                        categorized_1_hb_nbG_aa = categorized_1_hb_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_1_hb_nbG_aa
                        categorized_1_hb_nbU_aa = categorized_1_hb_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_1_hb_nbU_aa
                    elif cat == "CAT_2":
                        # category_storehb[1]+=1
                        #
                        # categorized_2_hb_nb[0]+=1 if nucleotide_base == "A" else categorized_2_hb_nb[0]
                        # categorized_2_hb_nb[1]+=1 if nucleotide_base == "C" else categorized_2_hb_nb[1]
                        # categorized_2_hb_nb[2]+=1 if nucleotide_base == "G" else categorized_2_hb_nb[2]
                        # categorized_2_hb_nb[3]+=1 if nucleotide_base == "U" else categorized_2_hb_nb[3]

                        # categorized_2_hb_nb = categorized_2_hb_nb + self.nucleotide_basecounter(nucleotide_base)

                        categorized_2_hb_nbA_aa = categorized_2_hb_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_2_hb_nbA_aa
                        categorized_2_hb_nbC_aa = categorized_2_hb_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_2_hb_nbC_aa
                        categorized_2_hb_nbG_aa = categorized_2_hb_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_2_hb_nbG_aa
                        categorized_2_hb_nbU_aa = categorized_2_hb_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_2_hb_nbU_aa
                    elif cat == "CAT_3":
                        # category_storehb[2]+=1
                        #
                        # categorized_3_hb_nb[0]+=1 if nucleotide_base == "A" else categorized_3_hb_nb[0]
                        # categorized_3_hb_nb[1]+=1 if nucleotide_base == "C" else categorized_3_hb_nb[1]
                        # categorized_3_hb_nb[2]+=1 if nucleotide_base == "G" else categorized_3_hb_nb[2]
                        # categorized_3_hb_nb[3]+=1 if nucleotide_base == "U" else categorized_3_hb_nb[3]

                        # categorized_3_hb_nb = categorized_3_hb_nb + self.nucleotide_basecounter(nucleotide_base)

                        categorized_3_hb_nbA_aa = categorized_3_hb_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_3_hb_nbA_aa
                        categorized_3_hb_nbC_aa = categorized_3_hb_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_3_hb_nbC_aa
                        categorized_3_hb_nbG_aa = categorized_3_hb_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_3_hb_nbG_aa
                        categorized_3_hb_nbU_aa = categorized_3_hb_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_3_hb_nbU_aa
                    elif cat == "CAT_4":
                        # category_storehb[3]+=1
                        #
                        # categorized_4_hb_nb[0]+=1 if nucleotide_base == "A" else categorized_4_hb_nb[0]
                        # categorized_4_hb_nb[1]+=1 if nucleotide_base == "C" else categorized_4_hb_nb[1]
                        # categorized_4_hb_nb[2]+=1 if nucleotide_base == "G" else categorized_4_hb_nb[2]
                        # categorized_4_hb_nb[3]+=1 if nucleotide_base == "U" else categorized_4_hb_nb[3]

                        categorized_4_hb_nbA_aa = categorized_4_hb_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_4_hb_nbA_aa
                        categorized_4_hb_nbC_aa = categorized_4_hb_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_4_hb_nbC_aa
                        categorized_4_hb_nbG_aa = categorized_4_hb_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_4_hb_nbG_aa
                        categorized_4_hb_nbU_aa = categorized_4_hb_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_4_hb_nbU_aa
                    elif cat == "CAT_5":
                        # category_storehb[4]+=1
                        #
                        # categorized_5_hb_nb[0]+=1 if nucleotide_base == "A" else categorized_5_hb_nb[0]
                        # categorized_5_hb_nb[1]+=1 if nucleotide_base == "C" else categorized_5_hb_nb[1]
                        # categorized_5_hb_nb[2]+=1 if nucleotide_base == "G" else categorized_5_hb_nb[2]
                        # categorized_5_hb_nb[3]+=1 if nucleotide_base == "U" else categorized_5_hb_nb[3]

                        # categorized_5_hb_nb = categorized_5_hb_nb + self.nucleotide_basecounter(nucleotide_base)

                        categorized_5_hb_nbA_aa = categorized_5_hb_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_5_hb_nbA_aa
                        categorized_5_hb_nbC_aa = categorized_5_hb_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_5_hb_nbC_aa
                        categorized_5_hb_nbG_aa = categorized_5_hb_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_5_hb_nbG_aa
                        categorized_5_hb_nbU_aa = categorized_5_hb_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_5_hb_nbU_aa
                    elif cat == "CAT_6":
                        # category_storehb[5]+=1
                        # categorized_6_hb_nb[0]+=1 if nucleotide_base == "A" else categorized_6_hb_nb[0]
                        # categorized_6_hb_nb[1]+=1 if nucleotide_base == "C" else categorized_6_hb_nb[1]
                        # categorized_6_hb_nb[2]+=1 if nucleotide_base == "G" else categorized_6_hb_nb[2]
                        # categorized_6_hb_nb[3]+=1 if nucleotide_base == "U" else categorized_6_hb_nb[3]

                        # categorized_6_hb_nb = categorized_6_hb_nb + self.nucleotide_basecounter(nucleotide_base)

                        categorized_6_hb_nbA_aa = categorized_6_hb_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_6_hb_nbA_aa
                        categorized_6_hb_nbC_aa = categorized_6_hb_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_6_hb_nbC_aa
                        categorized_6_hb_nbG_aa = categorized_6_hb_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_6_hb_nbG_aa
                        categorized_6_hb_nbU_aa = categorized_6_hb_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_6_hb_nbU_aa
                    elif cat == "CAT_7":
                        # category_storehb[6]+=1
                        # categorized_7_hb_nb[0]+=1 if nucleotide_base == "A" else categorized_7_hb_nb[0]
                        # categorized_7_hb_nb[1]+=1 if nucleotide_base == "C" else categorized_7_hb_nb[1]
                        # categorized_7_hb_nb[2]+=1 if nucleotide_base == "G" else categorized_7_hb_nb[2]
                        # categorized_7_hb_nb[3]+=1 if nucleotide_base == "U" else categorized_7_hb_nb[3]

                        # categorized_7_hb_nb = categorized_7_hb_nb + self.nucleotide_basecounter(nucleotide_base)

                        categorized_7_hb_nbA_aa = categorized_7_hb_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_7_hb_nbA_aa
                        categorized_7_hb_nbC_aa = categorized_7_hb_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_7_hb_nbC_aa
                        categorized_7_hb_nbG_aa = categorized_7_hb_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_7_hb_nbG_aa
                        categorized_7_hb_nbU_aa = categorized_7_hb_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_7_hb_nbU_aa
                    elif cat == "CAT_8":
                        # category_storehb[7]+=1
                        #
                        # categorized_8_hb_nb[0]+=1 if nucleotide_base == "A" else categorized_8_hb_nb[0]
                        # categorized_8_hb_nb[1]+=1 if nucleotide_base == "C" else categorized_8_hb_nb[1]
                        # categorized_8_hb_nb[2]+=1 if nucleotide_base == "G" else categorized_8_hb_nb[2]
                        # categorized_8_hb_nb[3]+=1 if nucleotide_base == "U" else categorized_8_hb_nb[3]

                        # categorized_8_hb_nb = categorized_8_hb_nb + self.nucleotide_basecounter(nucleotide_base)

                        categorized_8_hb_nbA_aa = categorized_8_hb_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_8_hb_nbA_aa
                        categorized_8_hb_nbC_aa = categorized_8_hb_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_8_hb_nbC_aa
                        categorized_8_hb_nbG_aa = categorized_8_hb_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_8_hb_nbG_aa
                        categorized_8_hb_nbU_aa = categorized_8_hb_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_8_hb_nbU_aa
                    elif cat == "CAT_9":
                        # category_storehb[8]+=1
                        # categorized_9_hb_nb[0]+=1 if nucleotide_base == "A" else categorized_9_hb_nb[0]
                        # categorized_9_hb_nb[1]+=1 if nucleotide_base == "C" else categorized_9_hb_nb[1]
                        # categorized_9_hb_nb[2]+=1 if nucleotide_base == "G" else categorized_9_hb_nb[2]
                        # categorized_9_hb_nb[3]+=1 if nucleotide_base == "U" else categorized_9_hb_nb[3]

                        # categorized_9_hb_nb = categorized_9_hb_nb + self.nucleotide_basecounter(nucleotide_base)

                        categorized_9_hb_nbA_aa = categorized_9_hb_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_9_hb_nbA_aa
                        categorized_9_hb_nbC_aa = categorized_9_hb_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_9_hb_nbC_aa
                        categorized_9_hb_nbG_aa = categorized_9_hb_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_9_hb_nbG_aa
                        categorized_9_hb_nbU_aa = categorized_9_hb_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_9_hb_nbU_aa
                elif hborvdw  == "vdw":
                    if cat == "CAT_1":
                        # category_storevdw[0]+=1

                        # categorized_1_vdw_nb = categorized_1_vdw_nb + self.nucleotide_basecounter(nucleotide_base)

                        # categorized_1_vdw_nb[0]+=1 if nucleotide_base == "A" else categorized_1_vdw_nb[0]
                        # categorized_1_vdw_nb[1]+=1 if nucleotide_base == "C" else categorized_1_vdw_nb[1]
                        # categorized_1_vdw_nb[2]+=1 if nucleotide_base == "G" else categorized_1_vdw_nb[2]
                        # categorized_1_vdw_nb[3]+=1 if nucleotide_base == "U" else categorized_1_vdw_nb[3]

                        categorized_1_vdw_nbA_aa = categorized_1_vdw_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_1_vdw_nbA_aa
                        categorized_1_vdw_nbC_aa = categorized_1_vdw_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_1_vdw_nbC_aa
                        categorized_1_vdw_nbG_aa = categorized_1_vdw_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_1_vdw_nbG_aa
                        categorized_1_vdw_nbU_aa = categorized_1_vdw_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_1_vdw_nbU_aa
                    elif cat == "CAT_2":
                        # category_storevdw[1]+=1

                        # categorized_2_vdw_nb = categorized_2_vdw_nb + self.nucleotide_basecounter(nucleotide_base)

                        # categorized_2_vdw_nb[0]+=1 if nucleotide_base == "A" else categorized_2_vdw_nb[0]
                        # categorized_2_vdw_nb[1]+=1 if nucleotide_base == "C" else categorized_2_vdw_nb[1]
                        # categorized_2_vdw_nb[2]+=1 if nucleotide_base == "G" else categorized_2_vdw_nb[2]
                        # categorized_2_vdw_nb[3]+=1 if nucleotide_base == "U" else categorized_2_vdw_nb[3]

                        categorized_2_vdw_nbA_aa = categorized_2_vdw_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_2_vdw_nbA_aa
                        categorized_2_vdw_nbC_aa = categorized_2_vdw_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_2_vdw_nbC_aa
                        categorized_2_vdw_nbG_aa = categorized_2_vdw_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_2_vdw_nbG_aa
                        categorized_2_vdw_nbU_aa = categorized_2_vdw_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_2_vdw_nbU_aa
                    elif cat == "CAT_3":
                        # category_storevdw[2]+=1

                        # categorized_3_vdw_nb = categorized_3_vdw_nb + self.nucleotide_basecounter(nucleotide_base)

                        # categorized_3_vdw_nb[0]+=1 if nucleotide_base == "A" else categorized_3_vdw_nb[0]
                        # categorized_3_vdw_nb[1]+=1 if nucleotide_base == "C" else categorized_3_vdw_nb[1]
                        # categorized_3_vdw_nb[2]+=1 if nucleotide_base == "G" else categorized_3_vdw_nb[2]
                        # categorized_3_vdw_nb[3]+=1 if nucleotide_base == "U" else categorized_3_vdw_nb[3]

                        categorized_3_vdw_nbA_aa = categorized_3_vdw_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_3_vdw_nbA_aa
                        categorized_3_vdw_nbC_aa = categorized_3_vdw_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_3_vdw_nbC_aa
                        categorized_3_vdw_nbG_aa = categorized_3_vdw_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_3_vdw_nbG_aa
                        categorized_3_vdw_nbU_aa = categorized_3_vdw_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_3_vdw_nbU_aa
                    elif cat == "CAT_4":
                        # category_storevdw[3]+=1

                        # categorized_4_vdw_nb = categorized_4_vdw_nb + self.nucleotide_basecounter(nucleotide_base)

                        # categorized_4_vdw_nb[0]+=1 if nucleotide_base == "A" else categorized_4_vdw_nb[0]
                        # categorized_4_vdw_nb[1]+=1 if nucleotide_base == "C" else categorized_4_vdw_nb[1]
                        # categorized_4_vdw_nb[2]+=1 if nucleotide_base == "G" else categorized_4_vdw_nb[2]
                        # categorized_4_vdw_nb[3]+=1 if nucleotide_base == "U" else categorized_4_vdw_nb[3]

                        categorized_4_vdw_nbA_aa = categorized_4_vdw_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_4_vdw_nbA_aa
                        categorized_4_vdw_nbC_aa = categorized_4_vdw_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_4_vdw_nbC_aa
                        categorized_4_vdw_nbG_aa = categorized_4_vdw_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_4_vdw_nbG_aa
                        categorized_4_vdw_nbU_aa = categorized_4_vdw_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_4_vdw_nbU_aa
                    elif cat == "CAT_5":
                        # category_storevdw[4]+=1

                        # categorized_5_vdw_nb = categorized_5_vdw_nb + self.nucleotide_basecounter(nucleotide_base)

                        # categorized_5_vdw_nb[0]+=1 if nucleotide_base == "A" else categorized_5_vdw_nb[0]
                        # categorized_5_vdw_nb[1]+=1 if nucleotide_base == "C" else categorized_5_vdw_nb[1]
                        # categorized_5_vdw_nb[2]+=1 if nucleotide_base == "G" else categorized_5_vdw_nb[2]
                        # categorized_5_vdw_nb[3]+=1 if nucleotide_base == "U" else categorized_5_vdw_nb[3]

                        categorized_5_vdw_nbA_aa = categorized_5_vdw_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_5_vdw_nbA_aa
                        categorized_5_vdw_nbC_aa = categorized_5_vdw_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_5_vdw_nbC_aa
                        categorized_5_vdw_nbG_aa = categorized_5_vdw_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_5_vdw_nbG_aa
                        categorized_5_vdw_nbU_aa = categorized_5_vdw_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_5_vdw_nbU_aa
                    elif cat == "CAT_6":
                        # category_storevdw[5]+=1

                        # categorized_6_vdw_nb = categorized_6_vdw_nb + self.nucleotide_basecounter(nucleotide_base)

                        # categorized_6_vdw_nb[0]+=1 if nucleotide_base == "A" else categorized_6_vdw_nb[0]
                        # categorized_6_vdw_nb[1]+=1 if nucleotide_base == "C" else categorized_6_vdw_nb[1]
                        # categorized_6_vdw_nb[2]+=1 if nucleotide_base == "G" else categorized_6_vdw_nb[2]
                        # categorized_6_vdw_nb[3]+=1 if nucleotide_base == "U" else categorized_6_vdw_nb[3]

                        categorized_6_vdw_nbA_aa = categorized_6_vdw_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_6_vdw_nbA_aa
                        categorized_6_vdw_nbC_aa = categorized_6_vdw_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_6_vdw_nbC_aa
                        categorized_6_vdw_nbG_aa = categorized_6_vdw_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_6_vdw_nbG_aa
                        categorized_6_vdw_nbU_aa = categorized_6_vdw_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_6_vdw_nbU_aa
                    elif cat == "CAT_7":
                        # category_storevdw[6]+=1

                        # categorized_7_vdw_nb = categorized_7_vdw_nb + self.nucleotide_basecounter(nucleotide_base)

                        # categorized_7_vdw_nb[0]+=1 if nucleotide_base == "A" else categorized_7_vdw_nb[0]
                        # categorized_7_vdw_nb[1]+=1 if nucleotide_base == "C" else categorized_7_vdw_nb[1]
                        # categorized_7_vdw_nb[2]+=1 if nucleotide_base == "G" else categorized_7_vdw_nb[2]
                        # categorized_7_vdw_nb[3]+=1 if nucleotide_base == "U" else categorized_7_vdw_nb[3]

                        categorized_7_vdw_nbA_aa = categorized_7_vdw_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_7_vdw_nbA_aa
                        categorized_7_vdw_nbC_aa = categorized_7_vdw_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_7_vdw_nbC_aa
                        categorized_7_vdw_nbG_aa = categorized_7_vdw_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_7_vdw_nbG_aa
                        categorized_7_vdw_nbU_aa = categorized_7_vdw_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_7_vdw_nbU_aa
                    elif cat == "CAT_8":
                        # category_storevdw[7]+=1
                        # categorized_8_vdw_nb = categorized_8_vdw_nb + self.nucleotide_basecounter(nucleotide_base)
                        # categorized_8_vdw_nb[0]+=1 if nucleotide_base == "A" else categorized_8_vdw_nb[0]
                        # categorized_8_vdw_nb[1]+=1 if nucleotide_base == "C" else categorized_8_vdw_nb[1]
                        # categorized_8_vdw_nb[2]+=1 if nucleotide_base == "G" else categorized_8_vdw_nb[2]
                        # categorized_8_vdw_nb[3]+=1 if nucleotide_base == "U" else categorized_8_vdw_nb[3]

                        categorized_8_vdw_nbA_aa = categorized_8_vdw_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_8_vdw_nbA_aa
                        categorized_8_vdw_nbC_aa = categorized_8_vdw_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_8_vdw_nbC_aa
                        categorized_8_vdw_nbG_aa = categorized_8_vdw_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_8_vdw_nbG_aa
                        categorized_8_vdw_nbU_aa = categorized_8_vdw_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_8_vdw_nbU_aa
                    elif cat == "CAT_9":
                        # category_storevdw[8]+=1
                        # categorized_9_vdw_nb = categorized_9_vdw_nb + self.nucleotide_basecounter(nucleotide_base)

                        # categorized_9_vdw_nb[0]+=1 if nucleotide_base == "A" else categorized_9_vdw_nb[0]
                        # categorized_9_vdw_nb[1]+=1 if nucleotide_base == "C" else categorized_9_vdw_nb[1]
                        # categorized_9_vdw_nb[2]+=1 if nucleotide_base == "G" else categorized_9_vdw_nb[2]
                        # categorized_9_vdw_nb[3]+=1 if nucleotide_base == "U" else categorized_9_vdw_nb[3]

                        categorized_9_vdw_nbA_aa = categorized_9_vdw_nbA_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "A" else categorized_9_vdw_nbA_aa
                        categorized_9_vdw_nbC_aa = categorized_9_vdw_nbC_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "C" else categorized_9_vdw_nbC_aa
                        categorized_9_vdw_nbG_aa = categorized_9_vdw_nbG_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "G" else categorized_9_vdw_nbG_aa
                        categorized_9_vdw_nbU_aa = categorized_9_vdw_nbU_aa + self.amino_acid_counter(amino_acid) if nucleotide_base == "U" else categorized_9_vdw_nbU_aa
            print categorized_file

        category_1_store_hb = sum( categorized_1_hb_nbA_aa + categorized_1_hb_nbC_aa + categorized_1_hb_nbG_aa + categorized_1_hb_nbU_aa)
        category_2_store_hb = sum( categorized_2_hb_nbA_aa + categorized_2_hb_nbC_aa + categorized_2_hb_nbG_aa + categorized_2_hb_nbU_aa)
        category_3_store_hb = sum( categorized_3_hb_nbA_aa + categorized_3_hb_nbC_aa + categorized_3_hb_nbG_aa + categorized_3_hb_nbU_aa)
        category_4_store_hb = sum( categorized_4_hb_nbA_aa + categorized_4_hb_nbC_aa + categorized_4_hb_nbG_aa + categorized_4_hb_nbU_aa)
        category_5_store_hb = sum( categorized_5_hb_nbA_aa + categorized_5_hb_nbC_aa + categorized_5_hb_nbG_aa + categorized_5_hb_nbU_aa)
        category_6_store_hb = sum( categorized_6_hb_nbA_aa + categorized_6_hb_nbC_aa + categorized_6_hb_nbG_aa + categorized_6_hb_nbU_aa)
        category_7_store_hb = sum( categorized_7_hb_nbA_aa + categorized_7_hb_nbC_aa + categorized_7_hb_nbG_aa + categorized_7_hb_nbU_aa)
        category_8_store_hb = sum( categorized_8_hb_nbA_aa + categorized_8_hb_nbC_aa + categorized_8_hb_nbG_aa + categorized_8_hb_nbU_aa)
        category_9_store_hb = sum( categorized_9_hb_nbA_aa + categorized_9_hb_nbC_aa + categorized_9_hb_nbG_aa + categorized_9_hb_nbU_aa)

        category_1_store_vdw = sum( categorized_1_vdw_nbA_aa + categorized_1_vdw_nbC_aa + categorized_1_vdw_nbG_aa + categorized_1_vdw_nbU_aa)
        category_2_store_vdw = sum( categorized_2_vdw_nbA_aa + categorized_2_vdw_nbC_aa + categorized_2_vdw_nbG_aa + categorized_2_vdw_nbU_aa)
        category_3_store_vdw = sum( categorized_3_vdw_nbA_aa + categorized_3_vdw_nbC_aa + categorized_3_vdw_nbG_aa + categorized_3_vdw_nbU_aa)
        category_4_store_vdw = sum( categorized_4_vdw_nbA_aa + categorized_4_vdw_nbC_aa + categorized_4_vdw_nbG_aa + categorized_4_vdw_nbU_aa)
        category_5_store_vdw = sum( categorized_5_vdw_nbA_aa + categorized_5_vdw_nbC_aa + categorized_5_vdw_nbG_aa + categorized_5_vdw_nbU_aa)
        category_6_store_vdw = sum( categorized_6_vdw_nbA_aa + categorized_6_vdw_nbC_aa + categorized_6_vdw_nbG_aa + categorized_6_vdw_nbU_aa)
        category_7_store_vdw = sum( categorized_7_vdw_nbA_aa + categorized_7_vdw_nbC_aa + categorized_7_vdw_nbG_aa + categorized_7_vdw_nbU_aa)
        category_8_store_vdw = sum( categorized_8_vdw_nbA_aa + categorized_8_vdw_nbC_aa + categorized_8_vdw_nbG_aa + categorized_8_vdw_nbU_aa)
        category_9_store_vdw = sum( categorized_9_vdw_nbA_aa + categorized_9_vdw_nbC_aa + categorized_9_vdw_nbG_aa + categorized_9_vdw_nbU_aa)

        total_nb_aa_pairs_hb  = category_1_store_hb + category_2_store_hb + category_3_store_hb + category_4_store_hb + category_5_store_hb + category_6_store_hb + category_7_store_hb + category_8_store_hb +category_9_store_hb
        total_nb_aa_pairs_vdw = category_1_store_vdw + category_2_store_vdw + category_3_store_vdw + category_4_store_vdw + category_5_store_vdw + category_6_store_vdw + category_7_store_vdw + category_8_store_vdw + category_9_store_vdw

        total_nucleotide_base_count_in_all_proteins, total_amino_acid_count_in_all_proteins = self.total_atom_counter()

        os.chdir(os.path.expanduser('~/bioresearch/statisticalpotentialfiles'))
        os.system("rm statisticalpotentialstore.csv")

        statistical_potential_cat_1_hb_nbA = self.develop_statistical_potential(categorized_1_hb_nbA_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_1_hb_nbC = self.develop_statistical_potential(categorized_1_hb_nbC_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_1_hb_nbG = self.develop_statistical_potential(categorized_1_hb_nbG_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_1_hb_nbU = self.develop_statistical_potential(categorized_1_hb_nbU_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_2_hb_nbA = self.develop_statistical_potential(categorized_2_hb_nbA_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_2_hb_nbC = self.develop_statistical_potential(categorized_2_hb_nbC_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_2_hb_nbG = self.develop_statistical_potential(categorized_2_hb_nbG_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_2_hb_nbU = self.develop_statistical_potential(categorized_2_hb_nbU_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_3_hb_nbA = self.develop_statistical_potential(categorized_3_hb_nbA_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_3_hb_nbC = self.develop_statistical_potential(categorized_3_hb_nbC_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_3_hb_nbG = self.develop_statistical_potential(categorized_3_hb_nbG_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_3_hb_nbU = self.develop_statistical_potential(categorized_3_hb_nbU_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_4_hb_nbA = self.develop_statistical_potential(categorized_4_hb_nbA_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_4_hb_nbC = self.develop_statistical_potential(categorized_4_hb_nbC_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_4_hb_nbG = self.develop_statistical_potential(categorized_4_hb_nbG_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_4_hb_nbU = self.develop_statistical_potential(categorized_4_hb_nbU_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_5_hb_nbA = self.develop_statistical_potential(categorized_5_hb_nbA_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_5_hb_nbC = self.develop_statistical_potential(categorized_5_hb_nbC_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_5_hb_nbG = self.develop_statistical_potential(categorized_5_hb_nbG_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_5_hb_nbU = self.develop_statistical_potential(categorized_5_hb_nbU_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_6_hb_nbA = self.develop_statistical_potential(categorized_6_hb_nbA_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_6_hb_nbC = self.develop_statistical_potential(categorized_6_hb_nbC_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_6_hb_nbG = self.develop_statistical_potential(categorized_6_hb_nbG_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_6_hb_nbU = self.develop_statistical_potential(categorized_6_hb_nbU_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_7_hb_nbA = self.develop_statistical_potential(categorized_7_hb_nbA_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_7_hb_nbC = self.develop_statistical_potential(categorized_7_hb_nbC_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_7_hb_nbG = self.develop_statistical_potential(categorized_7_hb_nbG_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_7_hb_nbU = self.develop_statistical_potential(categorized_7_hb_nbU_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_8_hb_nbA = self.develop_statistical_potential(categorized_8_hb_nbA_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_8_hb_nbC = self.develop_statistical_potential(categorized_8_hb_nbC_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_8_hb_nbG = self.develop_statistical_potential(categorized_8_hb_nbG_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_8_hb_nbU = self.develop_statistical_potential(categorized_8_hb_nbU_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_9_hb_nbA = self.develop_statistical_potential(categorized_9_hb_nbA_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_9_hb_nbC = self.develop_statistical_potential(categorized_9_hb_nbC_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_9_hb_nbG = self.develop_statistical_potential(categorized_9_hb_nbG_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_9_hb_nbU = self.develop_statistical_potential(categorized_9_hb_nbU_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_1_hb_nbA = self.develop_statistical_potential(categorized_1_vdw_nbA_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_1_vdw_nbC = self.develop_statistical_potential(categorized_1_vdw_nbC_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_1_vdw_nbG = self.develop_statistical_potential(categorized_1_vdw_nbG_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_1_vdw_nbU = self.develop_statistical_potential(categorized_1_vdw_nbU_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_2_vdw_nbA = self.develop_statistical_potential(categorized_2_vdw_nbA_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_2_vdw_nbC = self.develop_statistical_potential(categorized_2_vdw_nbC_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_2_vdw_nbG = self.develop_statistical_potential(categorized_2_vdw_nbG_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_2_vdw_nbU = self.develop_statistical_potential(categorized_2_vdw_nbU_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_3_vdw_nbA = self.develop_statistical_potential(categorized_3_vdw_nbA_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_3_vdw_nbC = self.develop_statistical_potential(categorized_3_vdw_nbC_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_3_vdw_nbG = self.develop_statistical_potential(categorized_3_vdw_nbG_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_3_vdw_nbU = self.develop_statistical_potential(categorized_3_vdw_nbU_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_4_vdw_nbA = self.develop_statistical_potential(categorized_4_vdw_nbA_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_4_vdw_nbC = self.develop_statistical_potential(categorized_4_vdw_nbC_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_4_vdw_nbG = self.develop_statistical_potential(categorized_4_vdw_nbG_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_4_vdw_nbU = self.develop_statistical_potential(categorized_4_vdw_nbU_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_5_vdw_nbA = self.develop_statistical_potential(categorized_5_vdw_nbA_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_5_vdw_nbC = self.develop_statistical_potential(categorized_5_vdw_nbC_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_5_vdw_nbG = self.develop_statistical_potential(categorized_5_vdw_nbG_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_5_vdw_nbU = self.develop_statistical_potential(categorized_5_vdw_nbU_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_6_vdw_nbA = self.develop_statistical_potential(categorized_6_vdw_nbA_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_6_vdw_nbC = self.develop_statistical_potential(categorized_6_vdw_nbC_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_6_vdw_nbG = self.develop_statistical_potential(categorized_6_vdw_nbG_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_6_vdw_nbU = self.develop_statistical_potential(categorized_6_vdw_nbU_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_7_vdw_nbA = self.develop_statistical_potential(categorized_7_vdw_nbA_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_7_vdw_nbC = self.develop_statistical_potential(categorized_7_vdw_nbC_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_7_vdw_nbG = self.develop_statistical_potential(categorized_7_vdw_nbG_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_7_vdw_nbU = self.develop_statistical_potential(categorized_7_vdw_nbU_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_8_vdw_nbA = self.develop_statistical_potential(categorized_8_vdw_nbA_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_8_vdw_nbC = self.develop_statistical_potential(categorized_8_vdw_nbC_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_8_vdw_nbG = self.develop_statistical_potential(categorized_8_vdw_nbG_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_8_vdw_nbU = self.develop_statistical_potential(categorized_8_vdw_nbU_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        statistical_potential_cat_9_vdw_nbA = self.develop_statistical_potential(categorized_9_vdw_nbA_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_9_vdw_nbC = self.develop_statistical_potential(categorized_9_vdw_nbC_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_9_vdw_nbG = self.develop_statistical_potential(categorized_9_vdw_nbG_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)
        statistical_potential_cat_9_vdw_nbU = self.develop_statistical_potential(categorized_9_vdw_nbU_aa,category_1_store_vdw,total_nb_aa_pairs_vdw,total_nucleotide_base_count_in_all_proteins,total_amino_acid_count_in_all_proteins)

        # print "HB  CAT_1:%s,CAT_2:%s,CAT_3:%s,CAT_4:%s,CAT_5:%s,CAT_6:%s,CAT_7:%s,CAT_8:%s,CAT_9:%s" %(category_storehb[0],category_storehb[1],category_storehb[2],category_storehb[3],category_storehb[4],category_storehb[5],category_storehb[6],category_storehb[7],category_storehb[8])
        # print "VDW CAT_1:%s,CAT_2:%s,CAT_3:%s,CAT_4:%s,CAT_5:%s,CAT_6:%s,CAT_7:%s,CAT_8:%s,CAT_9:%s" %(sum(categorized_1_vdw_nb),sum(categorized_2_vdw_nb),sum(categorized_3_vdw_nb),sum(categorized_4_vdw_nb),sum(categorized_5_vdw_nb),sum(categorized_6_vdw_nb),sum(categorized_7_vdw_nb),sum(categorized_8_vdw_nb),sum(categorized_9_vdw_nb))
        print "HB  CAT_1:%s,CAT_2:%s,CAT_3:%s,CAT_4:%s,CAT_5:%s,CAT_6:%s,CAT_7:%s,CAT_8:%s,CAT_9:%s" %(category_1_store_hb,category_2_store_hb,category_3_store_hb,category_4_store_hb,category_5_store_hb,category_6_store_hb,category_7_store_hb,category_8_store_hb,category_9_store_hb)
        print "VDW CAT_1:%s,CAT_2:%s,CAT_3:%s,CAT_4:%s,CAT_5:%s,CAT_6:%s,CAT_7:%s,CAT_8:%s,CAT_9:%s" %(category_1_store_vdw,category_2_store_vdw,category_3_store_vdw,category_4_store_vdw,category_5_store_vdw,category_6_store_vdw,category_7_store_vdw,category_8_store_vdw,category_9_store_vdw)

    def structureParser(self,structure_id,filename):
        # this function returns to the caller a structure object from the pdb file name by parsing the pdb file
        parser = PDBParser()
        structure = parser.get_structure(structure_id,filename)
        nucleotide_base_count_in_protein = np.array([0,0,0,0])
        amino_acid_count_in_protein = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        # The Structure object follows the so-called SMCRA (Structure/Model/Chain/Residue/Atom)
        for model in structure:
            # print model
            for chain in model:
                # print chain
                for residue in chain:
                    # print residue.get_resname()
                    if residue.get_resname().strip() == "G":
                      nucleotide_base_count_in_protein[0]+=1
                    elif residue.get_resname().strip() == "C":
                      nucleotide_base_count_in_protein[1]+=1
                    elif residue.get_resname().strip() == "A":
                      nucleotide_base_count_in_protein[2]+=1
                    elif residue.get_resname().strip() == "U":
                      nucleotide_base_count_in_protein[3]+=1

                    elif residue.get_resname().strip() == "ARG":
                        amino_acid_count_in_protein[0]+=1
                    elif residue.get_resname().strip() == "ALA":
                        amino_acid_count_in_protein[1]+=1
                    elif residue.get_resname().strip() == "ASN":
                        amino_acid_count_in_protein[2]+=1
                    elif residue.get_resname().strip() == "ASP":
                        amino_acid_count_in_protein[3]+=1
                    elif residue.get_resname().strip() == "GLN":
                        amino_acid_count_in_protein[4]+=1
                    elif residue.get_resname().strip() == "GLU":
                        amino_acid_count_in_protein[5]+=1
                    elif residue.get_resname().strip() == "GLY":
                        amino_acid_count_in_protein[6]+=1
                    elif residue.get_resname().strip() == "CYS":
                        amino_acid_count_in_protein[7]+=1
                    elif residue.get_resname().strip() == "HIS":
                        amino_acid_count_in_protein[8]+=1
                    elif residue.get_resname().strip() == "ILE":
                        amino_acid_count_in_protein[9]+=1
                    elif residue.get_resname().strip() == "LEU":
                        amino_acid_count_in_protein[10]+=1
                    elif residue.get_resname().strip() == "LYS":
                        amino_acid_count_in_protein[11]+=1
                    elif residue.get_resname().strip() == "MET":
                        amino_acid_count_in_protein[12]+=1
                    elif residue.get_resname().strip() == "PHE":
                        amino_acid_count_in_protein[13]+=1
                    elif residue.get_resname().strip() == "PRO":
                        amino_acid_count_in_protein[14]+=1
                    elif residue.get_resname().strip() == "SER":
                        amino_acid_count_in_protein[15]+=1
                    elif residue.get_resname().strip() == "THR":
                        amino_acid_count_in_protein[16]+=1
                    elif residue.get_resname().strip() == "TRP":
                        amino_acid_count_in_protein[17]+=1
                    elif residue.get_resname().strip() == "TYR":
                        amino_acid_count_in_protein[18]+=1
                    elif residue.get_resname().strip() == "VAL":
                        amino_acid_count_in_protein[19]+=1

        return nucleotide_base_count_in_protein, amino_acid_count_in_protein

    def total_atom_counter(self):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/pdbfiles'))
        listofcategorizedfiles = os.listdir('.')
        total_nucleotide_base_count_in_protein = [0,0,0,0]
        total_amino_acid_count_in_protein = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        for pdbfile in listofcategorizedfiles:
            lhs,rhs = pdbfile.split('.')
            pdb_name = lhs.strip('pdb')
            nucleotide_base_count_in_protein, amino_acid_count_in_protein = self.structureParser(pdb_name,pdbfile)
            total_nucleotide_base_count_in_protein = total_nucleotide_base_count_in_protein + nucleotide_base_count_in_protein
            total_amino_acid_count_in_protein = total_amino_acid_count_in_protein + amino_acid_count_in_protein
            # print pdb_name
            # print nucleotide_base_count_in_protein
            # print amino_acid_count_in_protein

        return total_nucleotide_base_count_in_protein, total_amino_acid_count_in_protein
        # print total_nucleotide_base_count_in_protein
        # print total_amino_acid_count_in_protein
        # print sum(total_nucleotide_base_count_in_protein)
        # print sum(total_amino_acid_count_in_protein)

    def develop_statistical_potential(self,categorized_1_hb_nbA_aa,category_1_store_hb,total_nb_aa_pairs_hb,total_nucleotide_base_count_in_protein,total_amino_acid_count_in_protein):
        # formula for statistical potential
        RT=8.31*273
        statistical_potential_cat = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        # print statistical_potential_cat
        # print categorized_1_hb_nbA_aa
        # print category_1_store_hb
        # print total_nb_aa_pairs_hb
        # print sum(total_nucleotide_base_count_in_protein)
        # print sum(total_amino_acid_count_in_protein)
        for i in range(0,len(categorized_1_hb_nbA_aa)):
            # print categorized_1_hb_nbA_aa[i]
            if categorized_1_hb_nbA_aa[i] == 0 :
                statistical_potential_cat[i] = round(0.000,3)
            else:
                first_denom=(float(sum(categorized_1_hb_nbA_aa))/sum(total_nucleotide_base_count_in_protein))
                second_denom=(float(categorized_1_hb_nbA_aa[i])/sum(total_amino_acid_count_in_protein))
                statistical_potential_cat[i] = round((float(categorized_1_hb_nbA_aa[i])/total_nb_aa_pairs_hb)/(first_denom+second_denom),3)
            # statistical_potential_cat [i] = (float(categorized_1_hb_nbA_aa[i])/total_nb_aa_pairs_hb)/(sum(categorized_1_hb_nbA_aa)/sum(total_nucleotide_base_count_in_protein)) + (category_1_store_hb[i]/sum(total_amino_acid_count_count_in_protein))
        # for nb_aa in categorized_1_hb_nbA_aa:
        #     (float(nb_aa)/total_nb_aa_pairs_hb)/( /total_nucleotide_base_count_in_protein) + ( / total_amino_acid_count_count_in_protein)
        # ((# of pair of the specific nucelotide_base to the specifc amino_acid in category 1) / (all nucleotide_base to amino acid pairs in learning set))/((# of specific nucleotide in category 1)/(# total of nucleotide bases in learning set ) * (# of amino_acids in category 1)/(# total of amino acids in learning set)
        # _____________________________________________________
        # a=sum(total_nucleotide_base_count_in_protein)
        # b=sum(total_amino_acid_count_in_protein)
        # statistical_potential_cat_x_hb = float()/
        # statistical_potential_cat_x_vdw = float()/

        # lets store the statistical potentials here and write it to a file here.
        # os.chdir(os.path.expanduser('~/bioresearch/statisticalpotentialfiles'))
        file = open("statisticalpotentialstore.csv","a")
        statistical_potential_line="%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n" %(statistical_potential_cat[0],statistical_potential_cat[1],statistical_potential_cat[2],statistical_potential_cat[3],statistical_potential_cat[4],statistical_potential_cat[5],statistical_potential_cat[6],statistical_potential_cat[7],statistical_potential_cat[8],statistical_potential_cat[9],statistical_potential_cat[10],statistical_potential_cat[11],statistical_potential_cat[12],statistical_potential_cat[13],statistical_potential_cat[14],statistical_potential_cat[15],statistical_potential_cat[16],statistical_potential_cat[17],statistical_potential_cat[18],statistical_potential_cat[19])
        file.write(statistical_potential_line)
        print statistical_potential_cat
        return statistical_potential_cat
        # statistical_potential_cat_1_hb_nbA
        # statistical_potential_cat_1_hb_nbC
        # statistical_potential_cat_1_hb_nbG
        # statistical_potential_cat_1_hb_nbU
        #
        # statistical_potential_cat_2_hb_nbA
        # statistical_potential_cat_2_hb_nbC
        # statistical_potential_cat_2_hb_nbG
        # statistical_potential_cat_2_hb_nbU
        #
        # statistical_potential_cat_3_hb_nbA
        # statistical_potential_cat_3_hb_nbC
        # statistical_potential_cat_3_hb_nbG
        # statistical_potential_cat_3_hb_nbU
        #
        # statistical_potential_cat_4_hb_nbA
        # statistical_potential_cat_4_hb_nbC
        # statistical_potential_cat_4_hb_nbG
        # statistical_potential_cat_4_hb_nbU
        #
        # statistical_potential_cat_5_hb_nbA
        # statistical_potential_cat_5_hb_nbC
        # statistical_potential_cat_5_hb_nbG
        # statistical_potential_cat_5_hb_nbU
        #
        # statistical_potential_cat_6_hb_nbA
        # statistical_potential_cat_6_hb_nbC
        # statistical_potential_cat_6_hb_nbG
        # statistical_potential_cat_6_hb_nbU
        #
        # statistical_potential_cat_7_hb_nbA
        # statistical_potential_cat_7_hb_nbC
        # statistical_potential_cat_7_hb_nbG
        # statistical_potential_cat_7_hb_nbU
        #
        # statistical_potential_cat_8_hb_nbA
        # statistical_potential_cat_8_hb_nbC
        # statistical_potential_cat_8_hb_nbG
        # statistical_potential_cat_8_hb_nbU
        #
        # statistical_potential_cat_9_hb_nbA
        # statistical_potential_cat_9_hb_nbC
        # statistical_potential_cat_9_hb_nbG
        # statistical_potential_cat_9_hb_nbU


    # def amino_acid_nucleotide_base_bond_type_calculator(self,nucleotide_base,amino_acid):
    #     # categorized=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    #     categorized_nb_to_aa=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    #     if nucleotide_base == "A":
    #         if amino_acid == "ARG":
    #             categorized_nb_to_aa[0]+=1
    #         elif amino_acid == "ALA":
    #             categorized_nb_to_aa[1]+=1
    #         elif amino_acid == "ASN":
    #             categorized_nb_to_aa[2]+=1
    #         elif amino_acid == "ASP":
    #             categorized_nb_to_aa[3]+=1
    #         elif amino_acid == "GLN":
    #             categorized_nb_to_aa[4]+=1
    #         elif amino_acid == "GLU":
    #             categorized_nb_to_aa[5]+=1
    #         elif amino_acid == "GLY":
    #             categorized_nb_to_aa[6]+=1
    #         elif amino_acid == "CYS":
    #             categorized_nb_to_aa[7]+=1
    #         elif amino_acid == "HIS":
    #             categorized_nb_to_aa[8]+=1
    #         elif amino_acid == "ILE":
    #             categorized_nb_to_aa[9]+=1
    #         elif amino_acid == "LEU":
    #             categorized_nb_to_aa[10]+=1
    #         elif amino_acid == "LYS":
    #             categorized_nb_to_aa[11]+=1
    #         elif amino_acid == "MET":
    #             categorized_nb_to_aa[12]+=1
    #         elif amino_acid == "PHE":
    #             categorized_nb_to_aa[13]+=1
    #         elif amino_acid == "PRO":
    #             categorized_nb_to_aa[14]+=1
    #         elif amino_acid == "SER":
    #             categorized_nb_to_aa[15]+=1
    #         elif amino_acid == "THR":
    #             categorized_nb_to_aa[16]+=1
    #         elif amino_acid == "TRP":
    #             categorized_nb_to_aa[17]+=1
    #         elif amino_acid == "TYR":
    #             categorized_nb_to_aa[18]+=1
    #         elif amino_acid == "VAL":
    #             categorized_nb_to_aa[19]+=1
    #         return nucletide_base, categorized_nb_to_aa
    #     elif nucleotide_base == "C":
    #         categorized[1][0][0]+=1
    #     elif nucleotide_base == "G":
    #         categorized[2][0][0]+=1
    #     elif nucleotide_base == "U":
    #         categorized[3][0][0]+=1
    #     return categorized

    #dont really need thi part either

    def nucleotide_basecounter(self,nucleotide_base):
        category_nucleotidebase=[0,0,0,0]
        if nucleotide_base == "A":
            category_nucleotidebase[0]+=1
        elif nucleotide_base == "C":
            category_nucleotidebase[1]+=1
        elif nucleotide_base == "G":
            category_nucleotidebase[2]+=1
        elif nucleotide_base == "U":
            category_nucleotidebase[3]+=1
        return category_nucleotidebase

    def amino_acid_counter(self,amino_acid):
        category_aminoacid=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        if amino_acid == "ARG":
            category_aminoacid[0]+=1
        elif amino_acid == "ALA":
            category_aminoacid[1]+=1
        elif amino_acid == "ASN":
            category_aminoacid[2]+=1
        elif amino_acid == "ASP":
            category_aminoacid[3]+=1
        elif amino_acid == "GLN":
            category_aminoacid[4]+=1
        elif amino_acid == "GLU":
            category_aminoacid[5]+=1
        elif amino_acid == "GLY":
            category_aminoacid[6]+=1
        elif amino_acid == "CYS":
            category_aminoacid[7]+=1
        elif amino_acid == "HIS":
            category_aminoacid[8]+=1
        elif amino_acid == "ILE":
            category_aminoacid[9]+=1
        elif amino_acid == "LEU":
            category_aminoacid[10]+=1
        elif amino_acid == "LYS":
            category_aminoacid[11]+=1
        elif amino_acid == "MET":
            category_aminoacid[12]+=1
        elif amino_acid == "PHE":
            category_aminoacid[13]+=1
        elif amino_acid == "PRO":
            category_aminoacid[14]+=1
        elif amino_acid == "SER":
            category_aminoacid[15]+=1
        elif amino_acid == "THR":
            category_aminoacid[16]+=1
        elif amino_acid == "TRP":
            category_aminoacid[17]+=1
        elif amino_acid == "TYR":
            category_aminoacid[18]+=1
        elif amino_acid == "VAL":
            category_aminoacid[19]+=1
        return category_aminoacid
