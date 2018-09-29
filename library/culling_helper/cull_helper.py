import os

class CullHelper:
    def cullhelper(self,structure_file_name):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/metadata_folder'))
        #open the culled_file
        culled_file = open("culled_proteins.txt")
        #skip the header
        lines = culled_file.readlines()[1:]
        #initalize an array which is okay and then sort
        protein_list = []
        for line in lines:
            protein_list.append(line[0:4])
        protein_list.sort()
        os.chdir(os.path.expanduser('~/bioresearch/compbio'))
        for protein_identifier in protein_list:
            final_file = open(structure_file_name,"a")
            final_file.write("\"%s\"\n" %(protein_identifier))
