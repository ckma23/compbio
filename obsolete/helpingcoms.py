import os as os
import csv as csv

from helper.biopythonpdbretriever import FileRetriever
from helper.PDBRestService import PDBRestServicehelper

def fileretriever ():
    with open('structures_of_interest.csv','rb') as csvfile:
        csvstore = csv.reader(csvfile,delimiter = ',')
        for i in csvstore:
            print (i)
            # note: had to strip the string as it was returning ['<string>'] for []
            FileRetriever().fileretrieving(''.join(i))

def payloaddescriberetriever ():
    with open('structures_of_interest.csv','rb') as csvfile:
        csvstore = csv.reader(csvfile,delimiter = ',')
        for i in csvstore:
            print (i)
            store = PDBRestServicehelper().describeMolRetriever(i)

def payloadstatusretriever ():
    with open('structures_of_interest.csv','rb') as csvfile:
        csvstore = csv.reader(csvfile,delimiter = ',')
        for i in csvstore:
            print (i)
            store = PDBRestServicehelper().statusSequenceRetriever(i)
            
#this grabs all the pdb files required
fileretriever()
#this uses the pdb rest service
payloaddescriberetriever()
#this uses the pdb rest service
payloadstatusretriever()
# this returns the list of structures interested
structure_array=FileRetriever().structure_of_interest_return()
print structure_array
