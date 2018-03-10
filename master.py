import os as os #import the python os library
import csv as csv #import the python csv library
import sys # import the sys to retrieve command line arguments
import re #provides regular expression matching
import json #import the biopython json library
from pprint import pprint #pprint for json
from Bio import *      #use the BioPython Python library
from Bio.PDB import *  #more specifically import the BioPython PDB library
from helper.biopythonpdbretriever import FileRetriever
from helper.PDBRestService import PDBRestServicehelper
from helper.hbpluscli import hbplusclihelper

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


def fileretriever():
    with open('structures_of_interest.csv','rb') as csvfile:
        csvstore = csv.reader(csvfile,delimiter = ',')
        for i in csvstore:
            print (i)
            # note: had to strip the string as it was returning ['<string>'] for []
            FileRetriever().fileretrieving(''.join(i))


def hbplushbcli():
  os.chdir('/Users/curtisma/bioresearch/compbio/pdbfiles')
  listoffiles = os.listdir('.')
  os.chdir('/Users/curtisma/bioresearch/')
  os.system('mkdir hbplusprocessedhbfiles')
  os.chdir('/Users/curtisma/bioresearch/hbplus')
  for f in listoffiles:
    stringprep=('~/bioresearch/compbio/pdbfiles/%s' %f)
    os.system('./hbplus %s' %stringprep)
    lhs,rhs=f.split(".",1)
    print lhs
    #output file is an .hb2
    os.system('mv %s.hb2 ~/bioresearch/hbplusprocessedhbfiles' %lhs)
    print f

# run hbplus in van der Waal mode
def hbplusvdwcli():
  os.chdir('/Users/curtisma/bioresearch/compbio/pdbfiles')
  listoffiles = os.listdir('.')
  os.chdir('/Users/curtisma/bioresearch/')
  os.system('mkdir hbplusprocessedvdwfiles')
  os.chdir('/Users/curtisma/bioresearch/hbplus')
  for f in listoffiles:
    stringprep=('~/bioresearch/compbio/pdbfiles/%s' %f)
    os.system('./hbplus %s -N' %stringprep)
    lhs,rhs=f.split(".",1)
    print lhs
    #output file is an .nb2
    os.system('mv %s.nb2 ~/bioresearch/hbplusprocessedvdwfiles' %lhs)
    print f

def dssrcli():
  os.chdir('/Users/curtisma/bioresearch/compbio/pdbfiles')
  listoffiles = os.listdir('.')
  os.chdir('/Users/curtisma/bioresearch/')
  os.system('mkdir DSSRprocessedfiles')
  for f in listoffiles:
    stringprep=('~/bioresearch/compbio/pdbfiles/%s' %f)
    os.chdir('/Users/curtisma/bioresearch/')
    os.system('./x3dna-dssr input=%s' %stringprep)
    pdblhs,pdbrhs=f.split(".",1)
    # print pdblhs
    os.system('mv dssr* /Users/curtisma/bioresearch/DSSRprocessedfiles')
    os.chdir('/Users/curtisma/bioresearch/DSSRprocessedfiles')
    dssrlistoffiles = os.listdir('.')
    newdssrlist = filter(re.compile("dssr-").match,dssrlistoffiles)
    for g in newdssrlist:
      dssrlhs,dssrrhs=g.split("-",1)
      dssrstringprep=('%s%s-%s' %(dssrlhs,pdblhs,dssrrhs))
      os.system('cp %s %s' %(g,dssrstringprep))
    #   print f

def hbplusprocessedhbreader():
  os.chdir('/Users/curtisma/bioresearch/hbplusprocessedhbfiles')
  listofprocessedfiles = os.listdir('.')
  listofprocessedfiles=["pdb1mnb.hb2"]
  totalHH=0
  totalMS=0
  totalMH=0
  totalSS=0
  totalSM=0
  totalaanb=0
  totalCGHH=0
  info1=[]
  info2=[]
  info3=[]
  info4=[]
  info5=[]
  info6=[]
  info7=[]
  info8=[]
  aminoacid_proteinstore=[]
  for hbplusprocessedfile in listofprocessedfiles:
    individualfileobject=open(hbplusprocessedfile)
    # print individualfileobject
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
            if (store[6:9].strip() == aa and store[20:23].strip() == nb) or (store[6:9].strip() == nb and store[20:23].strip() == aa):
              totalaanb += 1
              info1.append(store[0:6])
              info2.append(store[6:9])
              info3.append(store[9:12])
              info4.append(store[14:20])
              info5.append(store[20:23])
              info6.append(store[24:27])
    print hbplusprocessedfile
    print totalaanb
  for i in range(len(info1)):
    print "%s %s %s %s %s %s" %(info1[i],info2[i],info3[i],info4[i],info5[i],info6[i])

    #    contactcategorizer(value)
    #    if store[6:9].strip() == "A" and store[20:23].strip() == "U" and store[33:35] == "HH":
    #      totalAUHH += 1
    #    elif store[6:9].strip() == "C" and store[20:23].strip() == "G" and store[33:35] == "HH":
    #      totalCGHH += 1
    #    if store[33:35] == "HH":
    #      totalHH += 1
    #    elif store[33:35] == "MH":
    #      totalMH += 1
    #    elif store[33:35] == "MS":
    #      totalMS += 1
    #    elif store[33:35] == "SS":
    #      totalSS += 1
    #    elif store[33:35] == "SM":
    #      totalSM += 1
    #  fileprocesser(value)
    #  hbondcontactmatcher(info1,info2)
  # print info1
  # print info2
  # print info3
  # print info4
  # print info5
  # print info6
  # print info7
  # print info8
  # print "In these files %s the following was determined" %listofprocessedfiles
  # print "Total of %s A-U HHbonds" %totalAUHH
  # print "Total of %s C-G HHbonds" %totalCGHH
  # print "Total of %s HH bonds" %totalHH
  # print "Total of %s MH bonds" %totalMH
  # print "Total of %s MS bonds" %totalMS
  # print "Total of %s SS bonds" %totalSS
  # print "Total of %s SM bonds" %totalSM
# def hbondcontactmatcher(info1,info2):
#     if infor1 == info2:

def fileprocesser(value):
    info.append(value)

def dssrprocessedreader():
    os.chdir('/Users/curtisma/bioresearch')
    listofdssrprocessedfiles = os.listdir('.')
    listofdssrprocessedfiles = ["1mnb.json"]
    for i in listofdssrprocessedfiles:
        data = json.load(open(i))
        dssroutputcategories = ["stems","hairpins","torsions","stacks","splays","pairs","multiplets","helices","bulges","atom2bases"]
        for j in dssroutputcategories:
            print j
            try:
                # print(data[j])
                if j == "hairpins":
                  # print data[j][0]["nts_long"]
                  hairpinnt=data[j][0]["nts_long"]
                  print "%s %s" %(j,hairpinnt)
                elif j == "stems":
                  k=0 
                  while k < len(data[j]):
                    h=0
                    while h < len(data[j][k]["pairs"][0]):
                      sindex=data[j][0]["index"]
                      siindex = data[j][k]["pairs"][h]["index"]
                      snt1 = data[j][k]["pairs"][h]["nt1"]
                      snt2 = data[j][k]["pairs"][h]["nt2"]
                      print "%s %s %s %s %s" %(j,sindex, siindex, snt1,snt2)
                      h+=1
                    k+=1
                elif j == "helices":
                  # reminder when parsing helix structures must take into the account strand 1 and strand 2 these residue structures are the helices
                  k=0
                  while k < len(data[j]):
                    h=0
                    while h < len(data[j][k]["pairs"][0]):
                    # while h might need to loop depending upon number of helices
                      hindex=data[j][k]["index"]
                      hiindex = data[j][k]["pairs"][h]["index"]
                      hnt1 = data[j][k]["pairs"][h]["nt1"]
                      hnt2 = data[j][k]["pairs"][h]["nt2"]
                      print "%s %s %s %s %s" %(j,hindex, hiindex,hnt1,hnt2)
                      h+=1
                    k+=1
            except:
                print "There was an exception likely null"


def aminoacidrnamatcher(aminoacid,nucleotidebase):
    aminoacidlist=["ARG","ALA","ARG","GLY","CYS,""ILE","LYS","MET","PHE,""PRO","SER","THR","TYR","VAL"]
    nucleotidebase=["U","G","C","A"]
    for aa in aminoacidlist:
        for nb in nucleotidebase:
            print aa
            print nb

def helpoutput():
  print"\nWelcome to the compbio project help section\n \nThe following commands are available: \ndssrcli \nfileretriever \nhbplushbcli \nhbplusvdwcli \ndssrprocessedreader \nhbplusprocessedreader"
# testprotein = structureRetriever('1mnb','pdbfiles/pdb1mnb.ent')
# neighborSearcher(testprotein,3.0)
# aminoacidmatcher()

proinput=str(sys.argv[1])
if proinput == "help":
  helpoutput()
elif proinput == "dssrcli":
  dssrcli()
elif proinput == "fileretriever":
  fileretriever()
elif proinput == "hbplushbcli":
  hbplusclihelper().hbplushbcli()
elif proinput == "hbplusvdwcli":
  hbplusvdwcli()
elif proinput == "dssrparse":
  dssrprocessedreader()
elif proinput == "hbplushbparse":
  hbplusprocessedhbreader()


