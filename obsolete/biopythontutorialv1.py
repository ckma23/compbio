import os as os

os.chdir('/Users/curtisma/bioresearch/practicecompbio')

from Bio import SeqIO

for seq_record in SeqIO.parse("ls_orchid.fasta","fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))


def hbplusprocessedhbreader():
  os.chdir('/Users/curtisma/bioresearch/hbplusprocessedhbfiles')
  # listofprocessedfiles = os.listdir('.')
  listofprocessedfiles=["pdb1mnb.hb2"]
  totalHH=0
  totalMS=0
  totalMH=0
  totalSS=0
  totalSM=0
  totalaanb=0
  totalCGHH=0
  # info1=[]
  # info2=[]
  # info3=[]
  # info4=[]
  # info5=[]
  # info6=[]
  # info7=[]
  # info8=[]
  aminoacid_proteinstore=[]
  for hbplushbprocessedfile in listofprocessedfiles:
    #bring this inside to scale it across all files lets write this into a for loop
    # for k in range(0:8):
    #   info=info+k
    #   info=[]
    info1=[]
    info2=[]
    info3=[]
    info4=[]
    info5=[]
    info6=[]
    info7=[]
    info8=[]
    individualfileobject=open(hbplushbprocessedfile)
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
    print hbplushbprocessedfile
    print totalaanb
  hbplushbstorefilewriter(hbplushbprocessedfile,info1,info2,info3,info4,info5,info6) 
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

# def hbplushbtodssrcomparer():
#   os.chdir('/Users/curtisma/bioresearch/hbplushbsortedfiles')
#   listofprocessedhbplusfiles = os.listdir('.')
#   listofprocessedhbplusfiles = ["pdb1mnb.hb2.hbsorted"]
#   for hbplusfile in listofprocessedhbplusfiles:
#     hbplusfilestore = open(hbplusfile)
#     # os.chdir('/Users/curtisma/bioresearch/DSSRparsedfiles')
#     # listofprocesseddssrfiles = ["1mnb.dsr"]
#     # dssrfilestore = open(dssrfile)
#     for hbline in hbplusfilestore:
#       print hbline
#       # for dssrfile in listofprocesseddssrfiles:
#       #   dssrcomparer(hbline,dssrfile)

# def hbplushbtodssrcomparer():
#   os.chdir('/Users/curtisma/bioresearch/hbplushbsortedfiles')
#   listofprocessedhbplusfiles = os.listdir('.')
#   listofprocessedhbplusfiles = ["pdb1mnb.hb2.hbsorted"]
#   for hbplusfile in listofprocessedhbplusfiles:
#     hbplusfilestore = open(hbplusfile)
#     for hbline in hbplusfilestore:
#       dssrcomparer(hbline)