import os as os                                             #import the python os library
import csv as csv                                           #import the python csv library
import sys                                                  #import the sys to retrieve command line arguments
import re                                                   #provides regular expression matching
import json                                                 #import the biopython json library
import ConfigParser                                         #import the configParser library
import glob                                                 #import the glob library used for regression grabbing
import fnmatch                                              #import the fnmatch library
from pprint import pprint                                   #pprint for json
from Bio import *                                           #use the BioPython Python library
from Bio.PDB import *                                       #more specifically import the BioPython PDB library
from helper.biopythonpdbretriever import FileRetriever      #from the helper folder biopythonpdbretriever.py import FileRetriever class
from helper.PDBRestService import PDBRestServicehelper
from library.hbplusclilib.hbpluscli import hbplusclihelper
from library.dssrclilib.dssrcli import dssrclihelper
from library.dssrparser.dssrparser import DssrParser
from optparse import OptionParser # will need a option parsing library

#os.chdir(os.path.expanduser("~/bioresearch"))

def pdbfileretriever():
    with open('structures_of_interest.csv','rb') as csvfile:
        csvstore = csv.reader(csvfile,delimiter = ',')
        for i in csvstore:
            print (i)
            # note: had to strip the string as it was returning ['<string>'] for []
            FileRetriever().fileretrieving(''.join(i))

def hbplusprocessedreader(hborvdw):
  os.chdir(os.path.expanduser('~/bioresearch/hbplusprocessed%sfiles' %hborvdw))
  listofprocessedfiles = os.listdir('.')
  print listofprocessedfiles
  # listofprocessedfiles=["pdb1mnb.hb2","pdb1b2m.hb2"]
  # listofprocessedfiles=["pdb1mnb.nb2","pdb1b2m.nb2"]
  for hbplusprocessedfile in listofprocessedfiles:
    print hbplusprocessedfile
    os.chdir(os.path.expanduser('~/bioresearch/hbplusprocessed%sfiles' %hborvdw))
    info1=[]
    info2=[]
    info3=[]
    info4=[]
    info5=[]
    info6=[]
    info7=[]
    info8=[]
    individualfileobject=open(hbplusprocessedfile)
    lines=individualfileobject.readlines()
    i=8 #set the counter to skip the 8 lines of the header
    totallinesinfile=len(lines) #sum up the total lines in each file
    print len(lines)
    # totallinesinfile=20
    while i < totallinesinfile: #iterate through the file for each line in it
    #    print lines[i]
    # if info8="HH" then sum up the HH
       store=''.join(map(str,lines[i]))
    #    print store
    #    print len(store)
       i +=1
    #    info1.append(store[0:6])
    #    info2.append(store[6:9])
    #    info3.append(store[9:12])
    #    info4.append(store[14:20])
    #    info5.append(store[20:23])
    #    info6.append(store[24:27])
       info7.append(store[28:32])
       info8.append(store[33:35])
       aminoacidlist=["ARG","ALA","ASN","ASP","GLN","GLU","GLY","CYS","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
       nucleotidebase=["U","G","C","A"]
       for aa in aminoacidlist:
          for nb in nucleotidebase:
            # print aa
            # print nb
            # if store[6:9].strip() == aa and store[20:23].strip() == nb:
            if (store[6:9].strip() == aa and store[20:23].strip() == nb):
              info4.append(store[0:6])
              info5.append(store[6:9])
              info6.append(store[9:13])
              info1.append(store[14:20])
              info2.append(store[20:23])
              info3.append(store[24:27])
            #appears that we need to account if the amino acid or the nucleotide base is either side. however we need to fix the files so hbplus sorted always have nucelotide base is always on first column.
            elif (store[6:9].strip() == nb and store[20:23].strip() == aa):
              info1.append(store[0:6])
              info2.append(store[6:9])
              info3.append(store[9:13])
              info4.append(store[14:20])
              info5.append(store[20:23])
              info6.append(store[24:27])
    # print hbplusprocessedfile
    # this is messing up and when it goes fomr hb to vdw the line shrinks hence ti gets cut of
    # it's screwing up when writing the vdw sorted
    # Curtiss-MacBook-Pro:hbplushbvdwcombined curtisma$ cat pdb1h4s.hbplushbvdwsorted | grep U0029
    #U0029 T OP1 A0127- TRP  NE hb
    #U0029 T OP2 A0128- ARG  NH hb
    #U0029 T  P  A0127- TRP NE1 vdw
    #U0029 T  OP A0127- TRP NE1 vdw
    #U0029 T  OP A0127- TRP CE2 vdw
    # the spacing is messed up...
    print
    hbplusstorefilewriter(hbplusprocessedfile,info1,info2,info3,info4,info5,info6,hborvdw)
  # for i in range(len(info1)):
  #   print "%s %s %s %s %s %s" %(info1[i],info2[i],info3[i],info4[i],info5[i],info6[i])

def hbplusstorefilewriter(proteinname,col1,col2,col3,col4,col5,col6,hborvdw):
# def hbplushbstorefilewriter(proteinname,col1,col2,col3,col4,col5,col6):
  os.chdir("/Users/curtisma/bioresearch")
  os.system('mkdir hbplus%ssortedfiles' %hborvdw)
  # os.system('mkdir hbplushbsortedfiles')
  os.chdir("/Users/curtisma/bioresearch/hbplus%ssortedfiles" %hborvdw)
  filenamestring="%s.%ssorted" %(proteinname,hborvdw)
  os.system("rm %s" %filenamestring)
  file = open(filenamestring,"a")
  print "%s has been created" %filenamestring
  for i in range(len(col1)):
    #we need to clean the col1 and col 2 to make hbplus compatibile to dssr i.e. B0018- A..... means B strand Argenine18"
    # thus we need to clean it to be A0018 on the B strand"
    # remove the first letter from col because that is the strand information and we need the residue type information
    removedstrandstring=col1[i][1:]
    cleanednb=col2[i]+removedstrandstring.strip()
    # combine this into one and it'l be A00018-
    strandnb=col1[i]
    # still keep the strand information!
    strandnb=strandnb[0]
    #strip the space and the -
    cleanednb=hbplustodssrnbstringcleaner(cleanednb)
    file.write("%s %s %s %s %s %s %s\n" %(cleanednb,strandnb,col3[i],col4[i],col5[i],col6[i],hborvdw))

def hbplushbandvdwcombiner():
  print "combining the hbplushb sorted and hbplusvdw sorted into one file"
  os.chdir("/Users/curtisma/bioresearch")
  os.system("mkdir hbplushbvdwcombined")
  os.chdir("/Users/curtisma/bioresearch/hbplushbsortedfiles")
  listofprocessedhbplusfiles = os.listdir('.')
  print listofprocessedhbplusfiles
  for i in listofprocessedhbplusfiles:
    os.chdir("/Users/curtisma/bioresearch/hbplushbsortedfiles")
    print os.getcwd()
    individualfileobject=open(i)
    lhs,rhs=i.split(".",1)
    print lhs
    filenamestring="%s.hbplushbvdwsorted" %(lhs)
    os.chdir("/Users/curtisma/bioresearch/hbplushbvdwcombined")
    os.system("rm %s" %filenamestring)
    file = open(filenamestring,"a")
    for line in individualfileobject:
      file.write(line)
    os.chdir("/Users/curtisma/bioresearch/hbplusvdwsortedfiles")
    vdwfile="%s.nb2.vdwsorted" %(lhs)
    individualfileobject=open(vdwfile)
    filenamestring="%s.hbplushbvdwsorted" %(lhs)
    os.chdir("/Users/curtisma/bioresearch/hbplushbvdwcombined")
    file = open(filenamestring,"a")
    for line in individualfileobject:
      file.write(line)
#note to self stip only works at the beginning and end of the string

def hbplustodssrnbstringcleaner(nbtobecleaned):
  nbclean=nbtobecleaned
  nbclean=nbclean.strip()
  # nbclean=nbclean.replace("B","")
  nbclean=nbclean.strip('-')
  return nbclean

def hbplushbvdwtodssrcomparer():
  i=0
  os.chdir('/Users/curtisma/bioresearch/hbplushbvdwcombined')
  listofprocessedhbplusfiles = os.listdir('.')
  # listofprocessedhbplusfiles = ["pdb1mnb.hbplushbvdwsorted"]
  for hbplusfile in listofprocessedhbplusfiles:
    os.chdir('/Users/curtisma/bioresearch/hbplushbvdwcombined')
    print hbplusfile
    hbplusfilestore = open(hbplusfile)
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
        for dssrline in dssrfilestore:
          dssrcomparer(hbline,dssrline,filenamestring)

# def dssr_parsed_file_opener(pdbfile):
#     os.chdir(os.path.expanduser('~/bioresearch/DSSRparsedfiles'))
#     pdbfile.replace("pdb","")
#     pdbfile.replace("hbplushbvdwsorted","")
#     pdbfile+=".dsr"
#     dssrfilestore=open(pdbfile)


def dssrcomparer(hbline,dssrline,filenamestring):
  hblinecompare=hbline.split(' ')
  dssrcompare=dssrline.strip().split(' ')
  # this dssrline is returning a really long black space after it
  # print "%s hi" %dssrline

  #determine if first it is a helix
  # then determine if it it's a helix check if it's A form
  # then check if it's a back backbone
  #since it was not a helix then go to check
  # if this was a backbone or not.
  # if dssrcompare[0] in ["hairpins","bulges","loops"]:
  #     non_helix_form_comparer():
  # elif dssrcompare[0] == "helices" :
  #     a_form_checker():

  # note match the 2d structure, check that it is the same chain, then check if it's the same residue name from hbPlus and DSSR
  if (dssrcompare[0] in ["hairpins","bulges","loops"] and str(hblinecompare[1].strip()) == str(dssrcompare[1].strip()) and str(hblinecompare[0].strip()) == str(dssrcompare[2].strip())):
    result = non_helix_form_comparer(hblinecompare[2].strip())
    bondcategorizedwriter(filenamestring,result,hbline,dssrcompare)
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
        result = helix_comparer(dssrcompare[5].strip(),hblinecompare[2].strip(),dssrcompare[7].strip(),residue_type)
    except Exception as e:
        result = "error: %s" %e
    bondcategorizedwriter(filenamestring,result,hbline,dssrcompare)

def helix_comparer(aformmarker,nbatom,dssr_canonical_pair_determinant,nucleotide_residue_type):
    # will need to build a separate one for hydrogen bond and one for vanderwaals
    # category_store=[0,0,0,0,0,0,0,0,0]
    if aformmarker == "A":
        hbplacer = backbonechecker(nbatom)

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
                a_form_mg_placeholder = a_form_major_or_minor_groove_checker(nbatom,nucleotide_residue_type)
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
        not_a_form_placeholder = backbonechecker(nbatom)
        # return "CAT_4"
        if not_a_form_placeholder == "backbone":
            # category_store[6]+=1
            return "CAT_7"
        elif not_a_form_placeholder == "base":
            not_a_form_base_placeholder = not_a_form_checker(dssr_canonical_pair_determinant)

            if not_a_form_base_placeholder == "canonical_base_pair":
                not_a_form_base_mg_placeholder = not_a_form_major_or_minor_groove_checker(nbatom,nucleotide_residue_type)

                if not_a_form_base_mg_placeholder == "major_groove":
                    return "CAT_4"
                elif not_a_form_base_mg_placeholder == "not_major_groove":
                    return "CAT_6"
            elif not_a_form_base_placeholder == "not_canonical_base_pair":
                not_a_form_base_mg_placeholder = not_a_form_major_or_minor_groove_checker(nbatom,nucleotide_residue_type)

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

def not_a_form_checker (dssr_canonical_pair_determinant):
    if dssr_canonical_pair_determinant == "cW-W":
        base_pair_placeholder = "canonical_base_pair"
    else:
        base_pair_placeholder = "not_canonical_base_pair"
    return base_pair_placeholder

def a_form_major_or_minor_groove_checker (nbatom,nucleotide_residue_type):
    # this is assuming it was a canonical_base_pair
    if nucleotide_residue_type == "A":
        # major groove atoms in A are:
        # H61,H62,N6,C6,C5,N7,C8,H8
        # minor groove atoms in A are:
        # N2,C2,C5,C4
        if nbatom in ["H61","H62","N6","C6","C5","N7","C8","H8"]:
            nb_mg_placeholder = "major_groove"
        else:
            nb_mg_placeholder = "not_major_groove"
    elif nucleotide_residue_type == "U":
        # major groove atoms in U are:
        # O4,C4,C5,H5,C6,H6
        # minor groove atoms in U are:
        # H3,N3,C2,O2,N1
        if nbatom in ["O4","C4","C5","H5","C6","H6"]:
            nb_mg_placeholder = "major_groove"
        else:
            nb_mg_placeholder = "not_major_groove"
    elif nucleotide_residue_type == "C":
        # major groove atoms in C are:
        # H6,C6,C5,H5,C4,N4,H42,H41,N3
        # minor groove atoms in C are:
        # N1,C2,O2
        if nbatom in ["H6","C6","C5","H5","C4","N4","H42","H41","N3"]:
            nb_mg_placeholder = "major_groove"
        else:
            nb_mg_placeholder = "not_major_groove"
    elif nucleotide_residue_type == "G":
        # major groove atoms in G are:
        # H1,N2,C6,O6,C5,N7,C8,H8
        # minor groove atoms in G are:
        # N2,C2,C5,C4
        if nbatom in ["H1", "N2", "C6", "O6", "C5", "N7", "C8", "H8"]:
            nb_mg_placeholder = "major_groove"
        else:
            nb_mg_placeholder = "not_major_groove"
    # need to return this inside each one instead.
    # else:
    #     nb_mg_placeholder = "not_major_groove"
    return nb_mg_placeholder
    print "checking if the atom is in the major groove or minor groove"

def not_a_form_major_or_minor_groove_checker (nbatom,nucleotide_residue_type):
    # we need to check the residue type.... If it's in A or C then it'll be these atoms
    if nucleotide_residue_type == "A":
        # major groove atoms in A are:
        # H61,H62,N6,C6,C5,N7,C8,H8
        # minor groove atoms in A are:
        # N2,C2,C5,C4
        if nbatom in ["H61","H62","N6","C6","C5","N7","C8","H8"]:
            nb_mg_placeholder = "major_groove"
        else:
            nb_mg_placeholder = "not_major_groove"
    elif nucleotide_residue_type == "U":
        # major groove atoms in U are:
        # O4,C4,C5,H5,C6,H6
        # minor groove atoms in U are:
        # H3,N3,C2,O2,N1
        if nbatom in ["O4","C4","C5","H5","C6","H6"]:
            nb_mg_placeholder = "major_groove"
        else:
            nb_mg_placeholder = "not_major_groove"
    elif nucleotide_residue_type == "C":
        # major groove atoms in C are:
        # H6,C6,C5,H5,C4,N4,H42,H41,N3
        # minor groove atoms in C are:
        # N1,C2,O2
        if nbatom in ["H6","C6","C5","H5","C4","N4","H42","H41","N3"]:
            nb_mg_placeholder = "major_groove"
        else:
            nb_mg_placeholder = "not_major_groove"
    elif nucleotide_residue_type == "G":
        # major groove atoms in G are:
        # H1,N2,C6,O6,C5,N7,C8,H8
        # minor groove atoms in G are:
        # N2,C2,C5,C4
        if nbatom in ["H1","N2","C6","O6","C5","N7","C8","H8"]:
            nb_mg_placeholder = "major_groove"
        else:
            nb_mg_placeholder = "not_major_groove"
    # need to return this inside each one instead.
    # else:
    #     nb_mg_placeholder = "not_major_groove"
    return nb_mg_placeholder
    print "checking if the atom is in the major groove or minor groove"

def non_helix_form_comparer(backboneatom):
    nhfplaceholder = backbonechecker(backboneatom)
    if nhfplaceholder == "backbone":
        return "CAT_9"
    elif nhfplaceholder == "base":
        return "CAT_8"

def backbonechecker(backboneatom):
    if backboneatom in ["C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3"]:
        nb_placeholder = "backbone"
    else:
        nb_placeholder = "base"
    return nb_placeholder

### DECISION TREE FOR CATEGORIES ###
#CAT 1 : Helix => A Form ("A") => major grove => major groove like => canonical basepair (U to A; C to G) => check atom type
#CAT 2 : Helix => A Form ("A") => not major grove => not a canonical base pair (U to G) => check atom type
#CAT 3 : Helix => A Form ("A") => backbone ("C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3")
#CAT 4 : Helix => Not A Form ("B","Z",".","x") => major groove like => canonical basepair (U to A; C to G)
#CAT 5 : Helix => Not A Form ("B","Z",".","x") => major groove like => not a canonical base pair (U to G)
#CAT 6 : Helix => Not A Form ("B","Z",".","x") => not major groove like
#CAT 7 : Helix => Not A Form ("B","Z",".","x") => backbone ("C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3")
#CAT 8 : Not Helix (Hairpin, Bulge, Loops) => base
#CAT 9 : Not Helix (Hairpin, Bulge, Loops) => backbone ("C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3")

# assumptions is that A Form Helix must likely only be a canonical base pair or wobble pair.
# thus if it's not a canonical basepair then check
# CAT 1 or CAT 2
# check atom type in A-form and which atom types are major groove in A-form and which ones are not; thus CAT 1 or CAT 2
# CAT 4 or CAT 5 or CAT 6
# CAT 4: canonical base pair is CAT 4 by cW-w
# CAT 5: non cW-w and check the atom if it's in the major groove or minor groove
# CAT 6: non cW-w and atom is not in major groove


### DECISION TREE TWO FOR CATEGORIES ###
#CAT 1 : Helix => A Form ("A") => canonical basepair ("cW-W" where U to A; C to G; wobble) => check atom type for major groove or minor groove
#CAT 2 : Helix => A Form ("A") => not a canonical base pair or  canonical basepair atom type is on minor groove
#CAT 3 : Helix => A Form ("A") => backbone ("C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3")
#CAT 4 : Helix => Not A Form ("B","Z",".","x") => canonical basepair ("cW-W" where U to A; C to G; wobble) => check atom type for major groove or minor groove
#CAT 5 : Helix => Not A Form ("B","Z",".","x") => not a canonical base pair => check atom type for major groove or minor groove
#CAT 6 : Helix => Not A Form ("B","Z",".","x") => not major groove like
#CAT 7 : Helix => Not A Form ("B","Z",".","x") => backbone ("C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3")
#CAT 8 : Not Helix (Hairpin, Bulge, Loops) => base
#CAT 9 : Not Helix (Hairpin, Bulge, Loops) => backbone ("C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3")

## major groove atoms in A are: H61,H62,N6,C6,C5,N7,C8,H8
## minor groove atoms in A are: N2,C2,C5,C4
## major groove atoms in U are: O4,C4,C5,H5,C6,H6
## minor groove atoms in U are: H3,N3,C2,O2,N1
## major groove atoms in C are: H6,C6,C5,H5,C4,N4,H42,H41,N3
## minor groove atoms in C are: N1,C2,O2
## major groove atoms in G are: H1,N2,C6,O6,C5,N7,C8,H8
## minor groove atoms in G are: N2,C2,C5,C4

def bondcategorizedwriter(filenamestring,category,hbline,dssrcompare):
    file = open(filenamestring,"a")
    hbline=hbline.strip('\n')
    dssrcompare=' '.join(dssrcompare).strip('\n').strip()
    file.write("%s %s %s\n" %(category,hbline,dssrcompare))

def helpoutput():
  print "\nWelcome to the Computatinal Biology RNA-Protein Interaction help section\n"
  print "\nThe following commands are available:\n"
  print "\ndssrcli \nfileretriever \nhbplushbcli \nhbplusvdwcli \ndssrprocessedreader \nhbplusprocessedreader\n"

# testprotein = structureRetriever('1mnb','pdbfiles/pdb1mnb.ent')
# neighborSearcher(testprotein,3.0)
# aminoacidmatcher()

proinput=str(sys.argv[1])
if proinput == "help":
  helpoutput()
elif proinput == "pdbfileretriever":
  pdbfileretriever()
elif proinput == "hbplushbcli":
  hbplusclihelper().hbplushbcli()
elif proinput == "hbplusvdwcli":
  hbplusclihelper().hbplusvdwcli()
elif proinput == "hbplushbparse":
  hbplusprocessedreader("hb")
elif proinput == "hbplusvdwparse":
  hbplusprocessedreader("vdw")
elif proinput == "hbplushbvdwcombine":
  hbplushbandvdwcombiner()
elif proinput == "dssrcli":
  dssrclihelper().dssrcli()
elif proinput == "dssrparse":
  DssrParser().dssrprocessedreader()
elif proinput == "hbcategorizedssr":
  hbplushbvdwtodssrcomparer()
elif proinput == "completerun":
  hbplusprocessedreader("hb")
  hbplusprocessedreader("vdw")
  hbplushbandvdwcombiner()
  DssrParser().dssrprocessedreader()
  hbplushbvdwtodssrcomparer()
