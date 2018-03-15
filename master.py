import os as os                                             #import the python os library
import csv as csv                                           #import the python csv library
import sys                                                  #import the sys to retrieve command line arguments
import re                                                   #provides regular expression matching
import json                                                 #import the biopython json library
from pprint import pprint                                   #pprint for json
from Bio import *                                           #use the BioPython Python library
from Bio.PDB import *                                       #more specifically import the BioPython PDB library
from helper.biopythonpdbretriever import FileRetriever      #from the helper folder biopythonpdbretriever.py import FileRetriever class
from helper.PDBRestService import PDBRestServicehelper
from library.hbplusclilib.hbpluscli import hbplusclihelper
from library.dssrclilib.dssrcli import dssrclihelper
from optparse import OptionParser # will need a option parsing library

def structureRetriever(structure_id,filename):
    # this function returns to the caller a structure object from the pdb file name by parsing the pdb file and returns a structure object
    parser = PDBParser()
    structure = parser.get_structure(structure_id,filename)
    return structure

def structureParser(structure_id,filename):
    # this function returns to the caller a structure object from the pdb file name by parsing the pdb file
    parser = PDBParser()
    structure = parser.get_structure(structure_id,filename)
    # The Structure object follows the so-called SMCRA (Structure/Model/Chain/Residue/Atom)
    for model in structure:
        print model
        for chain in model:
            print chain
            for residue in chain:
                print residue
                for atom in residue:
                    print atom

def neighborsingleSearcher(structure,distance):
    atom_list = Selection.unfold_entities(structure, 'A')
    ns = NeighborSearch(atom_list)
    # print dir(atom_list[0]) #print all possible attributes or methods of the object
    # print atom_list[0].__dict__ #print all info of this object
    print len(atom_list)
    center = atom_list[0].get_coord()
    neighbors = ns.search(center,distance)
    residue_list = Selection.unfold_entities(neighbors, 'R')
    print residue_list
    print "There are %s of nearest residues", len(residue_list)
    for i in residue_list:
        print dir(i)
        print i.get_resname()
        print i.child_list
        # the next thing we can do is get the individual atoms in the residue, if it's a hydrogen atom then it's hydrogen bonding
        # we can start classifyign the contacts here if this is then incrementer counter to. If this is then increment the counter here.

def neighborSearcher(structure,distance):
    atom_list = Selection.unfold_entities(structure, 'A')
    ns = NeighborSearch(atom_list)
    # print dir(atom_list[0]) #print all possible attributes or methods of the object
    # print atom_list[0].__dict__ #print all info of this object
    indexing = [0,0,0,0]
    # for atom in atom_list:
    for atom in atom_list[0:100]:
        center = atom.get_coord()
        neighbors = ns.search(center,distance)
        residue_list = Selection.unfold_entities(neighbors, 'R')
        print atom
        print residue_list
        print "There are %s nearest residues." % len(residue_list)
        for residue in residue_list:
            # index = 0
            # indexing = [0,0,0,0]
            print residue.get_resname()
            # if residue.get_resname().strip() == "  G":
            # we have to clean this string up using .strip
            for atom in residue:
                print atom.get_name()
                if residue.get_resname().strip() == "G" and atom.get_name() == "H1":
                #print atom
                  print "matching G"
                #   index += 1
                  indexing[0] +=1
                elif residue.get_resname() == "C" and atom.get_name() == "H1":
                  print "matching C"
                  indexing[1] +=1
                elif residue.get_resname() == "A" and atom.get_name() == "H1":
                  print "matching A"
                  indexing[2] +=1
                elif residue.get_resname() == "U" and atom.get_name() == "H1":
                  print "matching U"
                  indexing[3] +=1
            # print "There are %s G to H1, %s C to H1, %s A to H1, and %s U to H1 interactions." %(indexing[0],indexing[1],indexing[2],indexing[3])
    print "There are %s G to H1, %s C to H1, %s A to H1, and %s U to H1 interactions." %(indexing[0],indexing[1],indexing[2],indexing[3])
            # print residue.child_list
        # the next thing we can do is get the individual atoms in the residue, if it's a hydrogen atom then it's hydrogen bonding
        # we can start classifyign the contacts here if this is then incrementer counter to. If this is then increment the counter here.
# run hbplus in hydrogenbond mode


def pdbfileretriever():
    with open('structures_of_interest.csv','rb') as csvfile:
        csvstore = csv.reader(csvfile,delimiter = ',')
        for i in csvstore:
            print (i)
            # note: had to strip the string as it was returning ['<string>'] for []
            FileRetriever().fileretrieving(''.join(i))

def hbplusprocessedreader(hborvdw):
  os.chdir('/Users/curtisma/bioresearch/hbplusprocessed%sfiles' %hborvdw)
  # listofprocessedfiles = os.listdir('.')
  listofprocessedfiles=["pdb1mnb.hb2"]
  for hbplusprocessedfile in listofprocessedfiles:
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
              info6.append(store[9:12])
              info1.append(store[14:20])
              info2.append(store[20:23])
              info3.append(store[24:27])
            #appears that we need to account if the amino acid or the nucleotide base is either side. however we need to fix the files so hbplus sorted always have nucelotide base is always on first column.
            elif (store[6:9].strip() == nb and store[20:23].strip() == aa):
              info1.append(store[0:6])
              info2.append(store[6:9])
              info3.append(store[9:12])
              info4.append(store[14:20])
              info5.append(store[20:23])
              info6.append(store[24:27])
    print hbplusprocessedfile
  hbplusstorefilewriter(hbplusprocessedfile,info1,info2,info3,info4,info5,info6,hborvdw)
  for i in range(len(info1)):
    print "%s %s %s %s %s %s" %(info1[i],info2[i],info3[i],info4[i],info5[i],info6[i])

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
  # WE HAVE TO ASSUME FOR NOW THAT RNA IS THE B STRAND!
  for i in range(len(col1)):
    #we need to clean the col1 and col 2 to make hbplus compatibile to dssr i.e. B0018- A means B strand Argenine18"
    cleanednb=col2[i]+col1[i].strip()
    # strandnb=col1[i[0]]
    strandnb=col1[i]
    strandnb=strandnb[0]
    cleanednb=hbplustodssrnbstringcleaner(cleanednb)
    file.write("%s %s %s %s %s %s %s\n" %(cleanednb,strandnb,col3[i],col4[i],col5[i],col6[i],hborvdw))

def hbplushbandvdwcombiner():
  print "combining the hbplushb sorted and hbplusvdw sorted into one file"
  os.chdir("/Users/curtisma/bioresearch")
  os.system("mkdir hbplushbvdwcombined")
  os.chdir("/Users/curtisma/bioresearch/hbplushbsortedfiles")
  listofprocessedhbplusfiles = os.listdir('.')
  for i in listofprocessedhbplusfiles:
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
  nbclean=nbclean.replace("B","")
  nbclean=nbclean.strip('-')
  return nbclean

def hbplushbtodssrcomparer():
  # os.chdir('/Users/curtisma/bioresearch/hbplushbsortedfiles')
  os.chdir('/Users/curtisma/bioresearch/hbplushbvdwcombined')
  listofprocessedhbplusfiles = os.listdir('.')
  listofprocessedhbplusfiles = ["pdb1mnb.hb2.hbsorted"]
  listofprocessedhbplusfiles = ["pdb1mnb.hbplushbvdwsorted"]
  for hbplusfile in listofprocessedhbplusfiles:
    hbplusfilestore = open(hbplusfile)
    for hbline in hbplusfilestore:
      dssrcomparer(hbline)

# this is getting called everytime, which is doesn't make sense in regards to making directories. Proof of concept this is okay though.
def dssrcomparer(hbline):
  os.chdir('/Users/curtisma/bioresearch')
  os.system('mkdir bondcategorized')
  os.chdir('/Users/curtisma/bioresearch/DSSRparsedfiles')
  listofprocesseddssrfiles = os.listdir('.')
  listofprocesseddssrfiles = ["1mnb.dsr"]
  for dssrfile in listofprocesseddssrfiles:
    # print dssrfile
    os.chdir('/Users/curtisma/bioresearch/DSSRparsedfiles')
    dssrfilestore = open(dssrfile)
    lhs,rhs=dssrfile.split(".",1)
    filenamestring="%s.bondcategorized" %(lhs)
    # categorizedbond=[]
    os.chdir('/Users/curtisma/bioresearch/bondcategorized')
    # os.system("rm %s" %filenamestring)
    for dssrline in dssrfilestore:
      hblinecompare=hbline.split(' ')
      dssrcompare=dssrline.split(' ')
      # print type(dssrcompare[0])
    #   print dssrline
    #   print dssrcompare
    #   print hblinecompare[0]
    #   print dssrcompare[2] #something here is coming back weird where it's a whole line.
    #   print dssrcompare[1]
    #   print hblinecompare[1]
      # if (dssrcompare[0] == "hairpins" and str(hblinecompare[0]) == str(dssrcompare[2])):
      # note match the 2d structure, check that it is the same chain, then check if it's the same residue name
      if (dssrcompare[0] == "hairpins" and str(hblinecompare[1].strip()) == str(dssrcompare[1].strip()) and str(hblinecompare[0].strip()) == str(dssrcompare[2].strip())):
        # print "This is a bond on a HAIRPIN! %s %s" %(hbline,dssrcompare)
        # hbline=hbline.strip('\n')
        # dssrcompare=' '.join(dssrcompare).strip('\n').strip()
        # file.write("Cat_1_hairpin %s %s\n" %(hbline,dssrcompare))
        bondcategorizedwriter(filenamestring,"CAT_1_HAIRPIN",hbline,dssrcompare)
      elif (dssrcompare[0] == "stems" and str(hblinecompare[1].strip()) == str(dssrcompare[1].strip()) and str(hblinecompare[0].strip()) == str(dssrcompare[2].strip())):
        # print "This is a bond on a STEM! %s %s" %(hbline,dssrcompare)
        # hbline=hbline.strip('\n')
        # dssrcompare=' '.join(dssrcompare).strip('\n').strip()
        # file.write("Cat_2_stem %s %s\n" %(hbline,dssrcompare))
        bondcategorizedwriter(filenamestring,"CAT_2_STEM",hbline,dssrcompare)
      elif (dssrcompare[0] == "helice" and str(hblinecompare[1].strip()) == str(dssrcompare[1].strip()) and str(hblinecompare[0].strip()) == str(dssrcompare[2].strip())):
        # hbline=hbline.strip('\n')
        # dssrcompare=' '.join(dssrcompare).strip('\n').strip()
        # file.write("Cat_3_helix %s %s\n" %(hbline,dssrcompare))
        bondcategorizedwriter(filenamestring,"CAT_3_HELIX",hbline,dssrcompare)

def bondcategorizedwriter(filenamestring,category,hbline,dssrcompare):
    file = open(filenamestring,"a")
    hbline=hbline.strip('\n')
    dssrcompare=' '.join(dssrcompare).strip('\n').strip()
    file.write("%s %s %s\n" %(category,hbline,dssrcompare))

#     for line in categorizedbond:
#       file.write("%s\n" %line)
#       os.chdir('/Users/curtisma/bioresearch/bondcategorized')
#       file = open(filenamestring,"a")
      # if (dssrcompare[0] == "hairpins" and str(hblinecompare[1].strip()) == str(dssrcompare[1].strip()) and str(hblinecompare[0].strip()) == str(dssrcompare[2].strip())):
      #   print "This is a vdw bond on a HAIRPIN! %s %s" %(hbline,dssrcompare)
      # elif (dssrcompare[0] == "stems" and str(hblinecompare[1].strip()) == str(dssrcompare[1].strip()) and str(hblinecompare[0].strip()) == str(dssrcompare[2].strip())):
      #   print "This is a vdw bond on a STEM! %s %s" %(hbline,dssrcompare)
      # print dssrline

# we need to add 0's back in... from dssr to hbplus comparison
def dssrtohbplusstringcleaner(dssnbtobecleaned):
  rhs = dssnbtobecleaned[3:]
  lhs = dssnbtobecleaned[2]
  if len(rhs) == 1:
    rhs="000"+rhs
  elif len(rhs) == 2:
    rhs="00"+rhs
  elif len(rhs) == 3:
    rhs.join("0",rhs)
    rhs="0"+rhs
  # elif len(rhs) ==4:
  #   rhs.join("0",rhs)
  #   rhs="0"+rhs
  # consider the case where the residue name is 5 digits 99999
  dssncleaned=lhs+rhs
  # print dssncleaned
  return dssncleaned

def dssrprocessedreader():
    os.chdir('/Users/curtisma/bioresearch')
    listofdssrprocessedfiles = os.listdir('.')
    listofdssrprocessedfiles = ["1mnb.json"]
    for i in listofdssrprocessedfiles:
        dssrstore=[]
        lhs,rhs=i.split(".",1)
        print lhs
        data = json.load(open(i))
        dssroutputcategories = ["stems","junction","hairpins","torsions","stacks","splays","pairs","multiplets","helices","bulges","atom2bases"]
        for j in dssroutputcategories:
            print j
            try:
                # print(data[j])
                if j == "hairpins":
                  hairpinresult=hairpindssrparser(data,j)
                  for line in hairpinresult:
                    dssrstore.append(line)
                elif j == "stems":
                  stemresult=dssrstemParser(data,j)
                  for line in stemresult:
                    dssrstore.append(line)
                elif j == "helices":
                  helixresult=dssrhelixParser(data,j)
                  for line in helixresult:
                    dssrstore.append(line)
                  # reminder when parsing helix structures must take into the account strand 1 and strand 2 these residue structures are the helices
            except:
                print "There was an exception likely null"
        # print dssrstore
        dssrstorefilewriter(lhs,dssrstore)


def dssrhelixParser(data,j):
  print "This will take a dssr json output and parse out the helixes"
  helixtoappend=[]
  k=0
  while k < len(data[j]):
    h=0
    while h < len(data[j][k]["pairs"][0]):
       # while h might need to loop depending upon number of helices
      hindex=data[j][k]["index"]
      hiindex = data[j][k]["pairs"][h]["index"]
      hnt1 = data[j][k]["pairs"][h]["nt1"]
      hnt2 = data[j][k]["pairs"][h]["nt2"]
      # print "%s %s %s %s %s" %(j,hindex, hiindex,hnt1,hnt2)
      helixtoappend.append("%s %s %s %s %s" %(j,hindex, hiindex,hnt1,hnt2))
      # helixtoapp="%s %s %s %s %s" %(j,hindex, hiindex,hnt1,hnt2)
      # print helixtoappend
      h+=1
    k+=1
  return helixtoappend


def dssrstemParser(data,j):
  stemtoappend=[]
  # this is stepping out to the loop, may want to try switching it over to for loop instead
  for k in range(0,len(data[j])):
    for h in range(0,len(data[j][k]["pairs"])):
      sindex=data[j][k]["index"]
      siindex = data[j][k]["pairs"][h]["index"]
      snt1 = data[j][k]["pairs"][h]["nt1"]
      snt2 = data[j][k]["pairs"][h]["nt2"]
      # print "%s %s %s %s %s" %(j,sindex, siindex, snt1,snt2)
      # snt1=snt1[2:]
      # snt2=snt2[2:]
      chain1=snt1[0]
      chain2=snt2[0]
      snt1cleaned = dssrtohbplusstringcleaner(snt1)
      snt2cleaned = dssrtohbplusstringcleaner(snt2)
      print "%s %s %s %s %s" %(j,chain1,snt1cleaned,sindex,siindex)
      stemtoappend.append("%s %s %s %s %s" %(j,chain1,snt1cleaned,sindex,siindex))
      stemtoappend.append("%s %s %s %s %s" %(j,chain2,snt2cleaned,sindex,siindex))
  return stemtoappend


#I should start calling individual handlers here for DSSR and build them here
def hairpindssrparser(data,j):
  # print data[j][0]["nts_long"]
  hairpintoappend=[]
  hairpinnt=data[j][0]["nts_long"]
  print "%s %s" %(j,hairpinnt)
  #hairpinnt is a comma separated string #"B.U13,B.C14,B.A15,B.U16,B.U17,B.A18"
  for hairpinnb in hairpinnt.split(","):
    #we are stripping B. for now lets assume RNA is always the Bstrand
    dssrnbaddedzeros = dssrtohbplusstringcleaner(hairpinnb)
    # print testing
    chain = hairpinnb[0]
    hairpinnb=hairpinnb[2:]
    hairpintoappend.append("%s %s %s"%(j,chain,dssrnbaddedzeros))
  return hairpintoappend
    #we need to add the 0s back in for hbplus...

#def stemdssrparser(data):

def dssrstorefilewriter(proteinname,dssrstore):
  os.chdir("/Users/curtisma/bioresearch")
  os.system('mkdir DSSRprocessedfiles')
  os.chdir("/Users/curtisma/bioresearch/DSSRparsedfiles")
  # os.system('touch dssrparsed.dsr')
  filenamestring="%s.dsr" %proteinname
  os.system("rm %s" %filenamestring )
  file = open(filenamestring,"a")
  print "%s has been created" %filenamestring
  for line in dssrstore:
    file.write("%s\n" %line)
  os.chdir('/Users/curtisma/bioresearch')

# we need to write this output into another file for now
#lets call this file dssr_cleaned
#secondarystructuretype residue

#compare the hbplus hydrogen bond residues against this file
# if residue = residue in dssr_cleaned file, bin this as Category
#Category 1 = Helix
#Category 2 = Stem
#Category 3 = Junction

def aminoacidrnamatcher(aminoacid,nucleotidebase):
    aminoacidlist=["ARG","ALA","ARG","GLY","CYS,""ILE","LYS","MET","PHE,""PRO","SER","THR","TYR","VAL"]
    nucleotidebase=["U","G","C","A"]
    for aa in aminoacidlist:
        for nb in nucleotidebase:
            print aa
            print nb

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
elif proinput == "dssrcli":
  dssrclihelper().dssrcli()
elif proinput == "dssrparse":
  dssrprocessedreader()
elif proinput == "hbplushbparse":
  hbplusprocessedreader("hb")
elif proinput == "hbplusvdwparse":
  hbplusprocessedreader("vdw")
elif proinput == "hbcategorizedssr":
  hbplushbtodssrcomparer()
elif proinput == "hbplushbvdwcombine":
  hbplushbandvdwcombiner()
