# this stores the complex (Bound) vs Protein(unbound) and Rna(unbound)
import os

class TestsetHashstore:
    def __init__(self):
        self.hashy = {
        "1WSU_A_E":{"protein":"1LVA_A","rna":"1MFK_A"},
        "2PJP_A_B":{"protein":"2PJP_A","rna":"1MFK_A"},
        "1LNG_A_B":{"protein":"1LNG_A","rna":"1Z43_A"},
        "1E7K_A_C":{"protein":"2JNB_A","rna":"1E7K_C"},
        "1WPU_A_C":{"protein":"1WPV_A","rna":"1WPU_C"},
        "2QUX_A_C":{"protein":"2QUD_A","rna":"2QUX_C"},
        "2JEA_A_C":{"protein":"2JE6_A","rna":"2JEA_C"},
        "2FMT_A_C":{"protein":"1FMT_A","rna":"3CW5_A"},
        "1MFQ_C_A":{"protein":"1QB2_B","rna":"1L9A_B"},
        "1U0B_B_A":{"protein":"1L17_A","rna":"1B23_R"},
        "1EC6_A_D":{"protein":"1DTJ_A","rna":"1EC6_D"},
        "1HC8_A_C":{"protein":"1FOY_A","rna":"1HC8_C"},
        "1JBR_B_D":{"protein":"1AQZ_A","rna":"1JBR_D"},
        "1KOG_A_I":{"protein":"1EVL_A","rna":"1KOG_I"},
        "1M8W_A_C":{"protein":"1M8Z_A","rna":"1M8W_C"},
        "1F7U_A_B":{"protein":"1BS2_A","rna":"1F7U_B"},
        "1K8W_A_B":{"protein":"1R3F_A","rna":"1K8W_B"},
        "1N78_A_C":{"protein":"1J09_A","rna":"1N78_C"},
        "1U63_A_B":{"protein":"1I2A_A","rna":"1U63_B"},
        "2BTE_A_B":{"protein":"1H3N_A","rna":"2BTE_B"},
        "2HW8_A_B":{"protein":"1AD2_A","rna":"2HW8_B"},
        }

    def testset_list_prepare(self,protein_or_rna,structure_file_name):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/metadata_folder'))
        testset_retrieve_list = TestsetHashstore().hashy
        protein_list_temp = []
        for element in testset_retrieve_list:
            value = testset_retrieve_list[element][protein_or_rna]
            value = value[0:4]
            protein_list_temp.append(value)
        os.chdir(os.path.expanduser('~/bioresearch/compbio'))

        for protein_identifier in protein_list_temp:
            final_file = open(structure_file_name,"a")
            final_file.write("\"%s\"\n" %(protein_identifier))

    def testset_list_chain_selector (self,protein_or_rna,structure_file_name):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/metadata_folder'))
        testset_retrieve_list = TestsetHashstore().hashy
        protein_list_temp = []
        for element in testset_retrieve_list:
            value =  testset_retrieve_list[element][protein_or_rna]
            value = value[0:4]
            protein_list_temp.append(value)
        os.chdir(os.path.expanduser('~/bioresearch/compbio'))

        for protein_identifier in protein_list_temp:
            final_file = open(structure_file_name,"a")
            final_file.write("\"%s\"\n" %(protein_identifier))

    def file_cleaner_based_on_chain(self,testset_directory,protein_or_rna):
        #go into the  directory currently where the testset pdb is stored and lets go strip the chains out
        os.chdir(os.path.expanduser('~/bioresearch/compbio/%s' %testset_directory))
        print testset_directory
        file_list = os.listdir('.')
        testset_retrieve_list = TestsetHashstore().hashy
        for file_to_be_stripped in file_list:
            print "processing new chain"
            os.chdir(os.path.expanduser('~/bioresearch/compbio/%s' %testset_directory))
            #n^2 its okay because it's not alot of things going on here.
            #technically dont need to use two for loops, can refactor
            for element in testset_retrieve_list.keys():
            #first match the name and find from the hash which chain to strip.
                if file_to_be_stripped[3:7].upper() == testset_retrieve_list[element][protein_or_rna][0:4]:
                    chain_to_keep = testset_retrieve_list[element][protein_or_rna][5]
                    # bound_complex_name = testset_retrieve_list[element][0:4].lower()
                    bound_complex_name = element[0:4]
                    print file_to_be_stripped
                    print bound_complex_name
                    print testset_retrieve_list[element][protein_or_rna][0:4]
                    print testset_retrieve_list[element][protein_or_rna][5]
                    print chain_to_keep
                    self.pdbfile_splitter_rna_protein(testset_directory,file_to_be_stripped,protein_or_rna,chain_to_keep,bound_complex_name.lower())

    def pdbfile_splitter_rna_protein(self,testset_directory,file_to_be_stripped,protein_or_rna,chain_to_keep,bound_complex_name):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip'))
        os.system("mkdir %s_seperated_pdbfiles_testset" %protein_or_rna)
        os.chdir(os.path.expanduser('~/bioresearch/compbio/%s' %testset_directory))
        # os.system("rm *.pdb")
        if protein_or_rna == "protein":
            atom_check = set(["ARG","ALA","ASN","ASP","GLN","GLU","GLY","CYS","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"])
        elif protein_or_rna == "rna":
            atom_check = set(["A","U","C","G"])
        with open(file_to_be_stripped) as pdblines:
            #1WSU_A_E first chain is the protein, second chain is the rna
            filenamestring="%s_%s.pdb" %(bound_complex_name,chain_to_keep)
            for line in pdblines:
                # noticing in the pdb file tha tthe nucloetide base is "   A" and needs to be stripped to match "A"
                if line[17:20].strip() in atom_check:
                    if line[0:4] == "ATOM" and line[21] == chain_to_keep:
                        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s_seperated_pdbfiles_testset' %protein_or_rna))
                        file = open(filenamestring,"a")
                        file.write(line)
                    elif line[0:4] == "TER " and line[21] == chain_to_keep:
                        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s_seperated_pdbfiles_testset' %protein_or_rna))
                        file = open(filenamestring,"a")
                        file.write(line)
                        break

    def pdbfile_native_fileprepper_rna_protein(self,testset_directory,protein_or_rna):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip'))
        os.system("mkdir cleaned_pdbfiles_testset_bound")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/%s' %testset_directory))
        filelist = os.listdir('.')
        testset_retrieve_list = TestsetHashstore().hashy
        if protein_or_rna == "protein":
            atom_check = set(["ARG","ALA","ASN","ASP","GLN","GLU","GLY","CYS","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"])
        elif protein_or_rna == "rna":
            atom_check = set(["A","U","C","G"])

        for file_to_be_stripped in filelist:
            for element in testset_retrieve_list.keys():
                if file_to_be_stripped[3:7].upper() == element[0:4]:
                    if protein_or_rna == "protein":
                        chain_to_keep = element[5]
                    elif protein_or_rna == "rna":
                        chain_to_keep = element[7]
                    bound_complex_name = element[0:4].lower()
                    os.chdir(os.path.expanduser('~/bioresearch/compbio/%s' %testset_directory))
                    with open(file_to_be_stripped) as pdblines:
                        filenamestring="%s_combined.pdb" %(bound_complex_name)
                        for line in pdblines:
                            if line[17:20].strip() in atom_check:
                                if line[0:4] == "ATOM" and line[21] == chain_to_keep:
                                    os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/cleaned_pdbfiles_testset_bound'))
                                    file = open(filenamestring,"a")
                                    file.write(line)
                                elif line[0:4] == "TER " and line[21] == chain_to_keep:
                                    os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/cleaned_pdbfiles_testset_bound'))
                                    file = open(filenamestring,"a")
                                    file.write(line)
                                    break
