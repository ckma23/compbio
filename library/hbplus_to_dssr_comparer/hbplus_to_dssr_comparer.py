import os as os                                             #import the python os library
import csv as csv                                           #import the python csv library
import sys                                                  #import the sys to retrieve command line arguments
import re                                                   #provides regular expression matching
import json                                                 #import the biopython json library
import ConfigParser                                         #import the configParser library
import glob                                                 #import the glob library used for regression grabbing
import fnmatch

class HbPlusToDssrComparer(object):

    def hbplushbvdwtodssrcomparer(self):
      i=0
      os.chdir('/Users/curtisma/bioresearch/hbplushbvdwcombined')
      listofprocessedhbplusfiles = os.listdir('.')
      # listofprocessedhbplusfiles = ["pdb1mnb.hbplushbvdwsorted"]
      for hbplusfile in listofprocessedhbplusfiles:
        os.chdir('/Users/curtisma/bioresearch/hbplushbvdwcombined')
        print hbplusfile
        hbplusfilestore = open(hbplusfile)
        hbplusfilestore.close
        os.chdir('/Users/curtisma/bioresearch')
        os.system('mkdir bondcategorized') #look into os.makedirs
        os.chdir(os.path.expanduser('~/bioresearch/bondcategorized'))
        hbplusfile = hbplusfile.replace("pdb","")
        hbplusfile = hbplusfile.replace("hbplushbvdwsorted","")
        hbplusfile+="bondcategorized"
        os.system('rm %s' %hbplusfile)
        # for hbline in hbplusfilestore:
        #     os.chdir('/Users/curtisma/bioresearch/DSSRparsedfiles')
        #     listofprocesseddssrfiles = os.listdir('.')
        #     print listofprocesseddssrfiles
        #     listofprocesseddssrfiles = ["pdb1mnb.dsr"]
        #     # listofprocesseddssrfiles = ""
        #     alhs,arhs=hbplusfile.split(".",1)
        #     print alhs
        #     for filematch in os.listdir('.'):
        #         print filematch
        #         alhs = "*" + alhs + "*"
        #         if fnmatch.fnmatch(filematch,alhs):
        #             "HIT"
        #             print listofprocesseddssrfiles
        #             listofprocesseddssrfiles = filematch
        #     #         print listofprocesseddssrfiles
        #     for dssrfile in listofprocesseddssrfiles:
        #       os.chdir('/Users/curtisma/bioresearch/DSSRparsedfiles')
        #       dssrfilestore = open(dssrfile)
        #       dssrfile=dssrfile.strip("pdb")
        #       lhs,rhs=dssrfile.split(".",1)
        #       filenamestring="%s.bondcategorized" %(lhs)
        #       os.chdir('/Users/curtisma/bioresearch/bondcategorized')
        #       for dssrline in dssrfilestore:
        #         dssrcomparer(hbline,dssrline,filenamestring)
        for hbline in hbplusfilestore:
            os.chdir('/Users/curtisma/bioresearch/DSSRparsedfiles')
            # listofprocesseddssrfiles = os.listdir('.')
            # listofprocesseddssrfiles = ""
            alhs,arhs=hbplusfile.split(".",1)
            for filematch in os.listdir('.'):
                alhs = "*" + alhs + "*"
                if fnmatch.fnmatch(filematch,alhs):
                    dssrfile = filematch
            # for dssrfile in listofprocesseddssrfiles:
            os.chdir('/Users/curtisma/bioresearch/DSSRparsedfiles')
            dssrfilestore = open(dssrfile)
            dssrfilestore.close
            storing=dssrfile.strip("pdb")
            lhs,rhs=storing.split(".",1)
            filenamestring="%s.bondcategorized" %(lhs)
            os.chdir('/Users/curtisma/bioresearch/bondcategorized')
            # print dssrfile
            # print dssrfilestore
            # i+=1
            # print i
            # print hbline
            # print filenamestring
            # try:
            #     lines=dssrfilestore.readlines()
            #     for line in lines:
            #     # it's not going through here
            #         print "HIT"
            #         print line
            # except error as e:
            #     print "%s" %e
            # # print "is it skipping this"
            # WE HAVE TO KEEP THIS A FOR LOOP in case an HB line matches more than once in dsr i.e. can be a hairpin or a helices
            for dssrline in dssrfilestore:
              self.dssrcomparer(hbline,dssrline,filenamestring)

    # def dssr_parsed_file_opener(pdbfile):
    #     os.chdir(os.path.expanduser('~/bioresearch/DSSRparsedfiles'))
    #     pdbfile.replace("pdb","")
    #     pdbfile.replace("hbplushbvdwsorted","")
    #     pdbfile+=".dsr"
    #     dssrfilestore=open(pdbfile)


    def dssrcomparer(self,hbline,dssrline,filenamestring):
      hblinecompare=hbline.split(' ')
      dssrcompare=dssrline.strip().split(' ')
      # this dssrline is returning a really long black space after it
      # print "%s hi" %dssrline
      # note match the 2d structure, check that it is the same chain, then check if it's the same residue name from hbPlus and DSSR
      if (dssrcompare[0] in ["hairpins","bulges","loops"] and str(hblinecompare[1].strip()) == str(dssrcompare[1].strip()) and str(hblinecompare[0].strip()) == str(dssrcompare[2].strip())):
        result = self.non_helix_form_comparer(hblinecompare[2].strip())
        self.bondcategorizedwriter(filenamestring,result,hbline,dssrcompare)
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
        self.bondcategorizedwriter(filenamestring,result,hbline,dssrcompare)

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

    def not_a_form_checker(self,dssr_canonical_pair_determinant):
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

    def bondcategorizedwriter(self,filenamestring,category,hbline,dssrcompare):
        final_file = open(filenamestring,"a")
        hbline=hbline.strip('\n')
        dssrcompare=' '.join(dssrcompare).strip('\n').strip()
        final_file.write("%s %s %s\n" %(category,hbline,dssrcompare))
        final_file.close
