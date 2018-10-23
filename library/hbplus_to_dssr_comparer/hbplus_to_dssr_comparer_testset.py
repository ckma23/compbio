import os as os                                             #import the python os library
import csv as csv                                           #import the python csv library
import sys                                                  #import the sys to retrieve command line arguments
import re                                                   #provides regular expression matching
import json                                                 #import the biopython json library
import ConfigParser                                         #import the configParser library
import glob                                                 #import the glob library used for regression grabbing
import fnmatch
import time


#might be ALOT faster if we combine the hb and dssr files one more time.
class HbPlusToDssrComparerTestset(object):

    def hbplushbvdwtodssrcomparer(self):
        # make the directory that the generated files will go into.
        os.system("mkdir ~/bioresearch/compbio/files_wip/bondcategorized_testset")
        i=0
        # need to grab all the protein names from the test set
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_testset'))
        list_of_protein_testset_complexes = os.listdir('.')
        list_of_protein_testset_complexes = ["1e7k"]
        for protein_testset_complex in list_of_protein_testset_complexes:
            #make the directory that the generated files will go into individual protein folders.
            os.system("mkdir ~/bioresearch/compbio/files_wip/bondcategorized_testset/%s" %protein_testset_complex)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/bondcategorized_testset/%s' %protein_testset_complex))
            #clean up the folder everytime we run this
            os.system('rm *.bondcategorized')
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_testset/%s' %protein_testset_complex))
            listofprocessedhbplusfiles = os.listdir('.')
            for hbplusfile in listofprocessedhbplusfiles:

                #MAY NEED TO KEEP THIS??
                os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/hbplus_hb_vdw_combined_testset/%s' %protein_testset_complex))
                # hbplusfilestore = open(hbplusfile)
                # hbplusfilestore.close
                # for each of the hbplus files please go ahead and make a bondcategorized file even if nothing matches
                filematch_lhs = hbplusfile
                filenamestring = "%s.bondcategorized" %(filematch_lhs)
                file_name_altogether = os.path.expanduser('~/bioresearch/compbio/files_wip/bondcategorized_testset/%s/%s' %(protein_testset_complex,filenamestring))
                create_file = open(file_name_altogether,"a")
                create_file.close

                with open(hbplusfile,'rb') as hbplus_meta_filestore:
                    hbplusfilestore = csv.reader(hbplus_meta_filestore)
                    for hbline in hbplusfilestore:
                        #MAY NEED TO KEEP THIS??
                        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/DSSRparsed_testset/%s'%protein_testset_complex))
                        # alhs,arhs=hbplusfile.split(".",1)
                        # for filematch in os.listdir('.'):
                        #     filename_lhs = hbplusfile + "*"
                        #     if fnmatch.fnmatch(filematch,filename_lhs):
                        #         dssrfile = filematch
                        filematch_lhs = hbplusfile
                        # filematch_lhs,throwaway = hbplusfile.split('.',1)
                        dssrfile = filematch_lhs + ".dsr"
                            # os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/DSSRparsed_testset/%s'%protein_testset_complex))
                            # os.chdir('/Users/curtisma/bioresearch/DSSRparsedfiles')
                        try:
                            with open(dssrfile,'rb') as dssrfilestore:
                                dssr_filestore = csv.reader(dssrfilestore)
                                # filename_lhs = filename_lhs.strip("*")
                                # storing=dssrfile.strip("pdb")
                                # lhs,rhs=storing.split(".",1)
                                filenamestring="%s.bondcategorized" %(filematch_lhs)
                                # os.chdir('/Users/curtisma/bioresearch/bondcategorized')
                    # WE HAVE TO KEEP THIS A FOR LOOP in case an HB line matches more than once in dsr i.e. can be a hairpin or a helices
                                for dssrline in dssr_filestore:
                                    self.dssrcomparer(hbline,dssrline,filenamestring,protein_testset_complex)
                        except Exception as e:
                            print "there was an error: %s" %e


    def dssrcomparer(self,hbline,dssrline,filenamestring,protein_testset_complex):
    #   hblinecompare=hbline.split(' ')
      hblinecompare=hbline
    #   dssrcompare=dssrline.strip().split(' ')
      dssrcompare = dssrline
      # this dssrline is returning a really long black space after it
      # print "%s hi" %dssrline
      # note match the 2d structure, check that it is the same chain, then check if it's the same residue name from hbPlus and DSSR
      if (dssrcompare[0] in ["hairpins","bulges","iloops","junctions"] and str(hblinecompare[1].strip()) == str(dssrcompare[1].strip()) and str(hblinecompare[0].strip()) == str(dssrcompare[2].strip())):
        result = self.non_helix_form_comparer(hblinecompare[2].strip())
        self.bondcategorizedwriter(filenamestring,result,hbline,dssrcompare,protein_testset_complex)
      # deprecating stems
      # elif (dssrcompare[0] == "stems" and str(hblinecompare[1].strip()) == str(dssrcompare[1].strip()) and str(hblinecompare[0].strip()) == str(dssrcompare[2].strip())):
      #   print dssrcompare[0]
      elif (dssrcompare[0] == "helices" and hblinecompare[1].strip() == dssrcompare[1].strip() and hblinecompare[0].strip() == dssrcompare[2].strip()):
        try:
            # send the cW-W which is dssrcompare[7].strip()
            # send the residue type which is hblinecompare[0][0]
            residue_type = hblinecompare[0]
            residue_type = residue_type[0]
            # print residue_type
            result = self.helix_comparer(dssrcompare[5].strip(),hblinecompare[2].strip(),dssrcompare[7].strip(),residue_type)
        except Exception as e:
            result = "error: %s" %e
            print result
        self.bondcategorizedwriter(filenamestring,result,hbline,dssrcompare,protein_testset_complex)
        
    #the last thing we might want to match is if we can't find anything at all!
    #   else:
    #       result = self.non_helix_form_comparer(hblinecompare[2].strip())
    #       result = self.nocategory_comparer(hblinecompare[2].strip())
    #       self.bondcategorizedwriter(filenamestring,result,hbline,dssrcompare)
    # def nocategory_comparer(self,backboneatom):
    #     nhfplaceholder = self.backbonechecker(backboneatom)
    #     if nhfplaceholder == "backbone":
    #         return "CAT_10"
    #     elif nhfplaceholder == "base":
    #         return "CAT_11"

    def helix_comparer(self,aformmarker,nbatom,dssr_canonical_pair_determinant,nucleotide_residue_type):
        # will need to build a separate one for hydrogen bond and one for vanderwaals
        # category_store=[0,0,0,0,0,0,0,0,0]
        if aformmarker == "A":
            hbplacer = self.backbonechecker(nbatom)

            if hbplacer == "backbone":
                # category_store[2]+=1
                return "CAT_3"
            elif hbplacer == "base":
                # should implement a check for if it was a canonical_base_pair CAT 1 and CAT 2 can not be as simple as if it's cW-W... then CAT 1 or CAT 2
                # a_form_mg_placeholder = a_form_major_or_minor_groove_checker(nbatom,nucleotide_residue_type)
                # if a_form_mg_placeholder == "major_groove":
                #     # category_store[0]+=1
                #     return "CAT_1"
                # elif a_form_mg_placeholder == "not_major_groove":
                #     # category_store[1]+=1
                #     return "CAT_2"
                if dssr_canonical_pair_determinant == "cW-W":
                    a_form_mg_placeholder = self.a_form_major_or_minor_groove_checker(nbatom,nucleotide_residue_type)
                    if a_form_mg_placeholder == "major_groove":
                        # category_store[0]+=1
                        return "CAT_1"
                    elif a_form_mg_placeholder == "not_major_groove":
                        return "CAT_2"
                else:
                    # category_store[1]+=1
                    return "CAT_2"
                # return "N/A"
        elif aformmarker in ["B","Z",".","x"]:
            not_a_form_placeholder = self.backbonechecker(nbatom)
            # return "CAT_4"
            if not_a_form_placeholder == "backbone":
                # category_store[6]+=1
                return "CAT_7"
            elif not_a_form_placeholder == "base":
                not_a_form_base_placeholder = self.not_a_form_checker(dssr_canonical_pair_determinant)

                if not_a_form_base_placeholder == "canonical_base_pair":
                    not_a_form_base_mg_placeholder = self.not_a_form_major_or_minor_groove_checker(nbatom,nucleotide_residue_type)

                    if not_a_form_base_mg_placeholder == "major_groove":
                        return "CAT_4"
                    elif not_a_form_base_mg_placeholder == "not_major_groove":
                        return "CAT_6"
                elif not_a_form_base_placeholder == "not_canonical_base_pair":
                    not_a_form_base_mg_placeholder = self.not_a_form_major_or_minor_groove_checker(nbatom,nucleotide_residue_type)

                    if not_a_form_base_mg_placeholder == "major_groove":
                        return "CAT_5"
                    elif not_a_form_base_mg_placeholder == "not_major_groove":
                        return "CAT_6"
                # not_a_form_base_placeholder = not_a_form_checker(dssr_canonical_pair_determinant)
                #
                # if not_a_form_base_placeholder == "canonical_base_pair":
                #     # category_store[1]+=1
                #     return "CAT_4"
                # elif not_a_form_base_placeholder == "not_canonical_base_pair":
                #     not_a_form_base_mg_placeholder = not_a_form_major_or_minor_groove_checker(nbatom,nucleotide_residue_type)
                #
                #     if not_a_form_base_mg_placeholder == "major_groove":
                #         return "CAT_5"
                #     elif not_a_form_base_mg_placeholder == "not_major_groove":
                #         return "CAT_6"
                # not_a_form_base_placeholder = not_a_form_checker(dssr_canonical_pair_determinant)
                #
                # if not_a_form_base_placeholder == "canonical_base_pair":
                #     # category_store[1]+=1
                #     return "CAT_4"
                # elif not_a_form_base_placeholder == "not_canonical_base_pair":
                #     not_a_form_base_mg_placeholder = not_a_form_major_or_minor_groove_checker(nbatom,nucleotide_residue_type)
                #
                #     if not_a_form_base_mg_placeholder == "major_groove":
                #         return "CAT_5"
                #     elif not_a_form_base_mg_placeholder == "not_major_groove":
                #         return "CAT_6"
        elif aformmarker in ["end"]:
                return "SHEAR"

    def not_a_form_checker (self,dssr_canonical_pair_determinant):
        if dssr_canonical_pair_determinant == "cW-W":
            base_pair_placeholder = "canonical_base_pair"
        else:
            base_pair_placeholder = "not_canonical_base_pair"
        return base_pair_placeholder

    def a_form_major_or_minor_groove_checker (self,nbatom,nucleotide_residue_type):
        # this is assuming it was a canonical_base_pair
        if nucleotide_residue_type == "A":
            # major groove atoms in A are:
            # H61,H62,N6,C6,C5,N7,C8,H8
            # minor groove atoms in A are:
            # N2,C2,C5,C4
            if nbatom in set(["H61","H62","N6","C6","C5","N7","C8","H8"]):
                nb_mg_placeholder = "major_groove"
            else:
                nb_mg_placeholder = "not_major_groove"
        elif nucleotide_residue_type == "U":
            # major groove atoms in U are:
            # O4,C4,C5,H5,C6,H6
            # minor groove atoms in U are:
            # H3,N3,C2,O2,N1
            if nbatom in set(["O4","C4","C5","H5","C6","H6"]):
                nb_mg_placeholder = "major_groove"
            else:
                nb_mg_placeholder = "not_major_groove"
        elif nucleotide_residue_type == "C":
            # major groove atoms in C are:
            # H6,C6,C5,H5,C4,N4,H42,H41,N3
            # minor groove atoms in C are:
            # N1,C2,O2
            if nbatom in set(["H6","C6","C5","H5","C4","N4","H42","H41","N3"]):
                nb_mg_placeholder = "major_groove"
            else:
                nb_mg_placeholder = "not_major_groove"
        elif nucleotide_residue_type == "G":
            # major groove atoms in G are:
            # H1,N2,C6,O6,C5,N7,C8,H8
            # minor groove atoms in G are:
            # N2,C2,C5,C4
            if nbatom in set(["H1", "N2", "C6", "O6", "C5", "N7", "C8", "H8"]):
                nb_mg_placeholder = "major_groove"
            else:
                nb_mg_placeholder = "not_major_groove"
        # need to return this inside each one instead.
        # else:
        #     nb_mg_placeholder = "not_major_groove"
        return nb_mg_placeholder
        print "checking if the atom is in the major groove or minor groove"

    def not_a_form_major_or_minor_groove_checker (self,nbatom,nucleotide_residue_type):
        # we need to check the residue type.... If it's in A or C then it'll be these atoms
        if nucleotide_residue_type == "A":
            # major groove atoms in A are:
            # H61,H62,N6,C6,C5,N7,C8,H8
            # minor groove atoms in A are:
            # N2,C2,C5,C4
            if nbatom in set(["H61","H62","N6","C6","C5","N7","C8","H8"]):
                nb_mg_placeholder = "major_groove"
            else:
                nb_mg_placeholder = "not_major_groove"
        elif nucleotide_residue_type == "U":
            # major groove atoms in U are:
            # O4,C4,C5,H5,C6,H6
            # minor groove atoms in U are:
            # H3,N3,C2,O2,N1
            if nbatom in set(["O4","C4","C5","H5","C6","H6"]):
                nb_mg_placeholder = "major_groove"
            else:
                nb_mg_placeholder = "not_major_groove"
        elif nucleotide_residue_type == "C":
            # major groove atoms in C are:
            # H6,C6,C5,H5,C4,N4,H42,H41,N3
            # minor groove atoms in C are:
            # N1,C2,O2
            if nbatom in set(["H6","C6","C5","H5","C4","N4","H42","H41","N3"]):
                nb_mg_placeholder = "major_groove"
            else:
                nb_mg_placeholder = "not_major_groove"
        elif nucleotide_residue_type == "G":
            # major groove atoms in G are:
            # H1,N2,C6,O6,C5,N7,C8,H8
            # minor groove atoms in G are:
            # N2,C2,C5,C4
            if nbatom in set(["H1","N2","C6","O6","C5","N7","C8","H8"]):
                nb_mg_placeholder = "major_groove"
            else:
                nb_mg_placeholder = "not_major_groove"
        # need to return this inside each one instead.
        # else:
        #     nb_mg_placeholder = "not_major_groove"
        return nb_mg_placeholder
        print "checking if the atom is in the major groove or minor groove"

    def non_helix_form_comparer(self,backboneatom):
        nhfplaceholder = self.backbonechecker(backboneatom)
        if nhfplaceholder == "backbone":
            return "CAT_9"
        elif nhfplaceholder == "base":
            return "CAT_8"

    def backbonechecker(self,backboneatom):
        if backboneatom in set(["C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3"]):
            nb_placeholder = "backbone"
        else:
            nb_placeholder = "base"
        return nb_placeholder

    def bondcategorizedwriter(self,filenamestring,category,hbline,dssrcompare,protein_testset_complex):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/bondcategorized_testset/%s'%protein_testset_complex))
        final_file = open(filenamestring,"a")
        # hbline=hbline.strip('\n')
        hbline=','.join(hbline)
        # dssrcompare=' '.join(dssrcompare).strip('\n').strip()
        dssrcompare=','.join(dssrcompare)
        final_file.write("%s,%s,%s\n" %(category,hbline,dssrcompare))
        #need to confirm why there is a timer here
        # time.sleep(.005)
        # time.sleep(.05)
        final_file.close
