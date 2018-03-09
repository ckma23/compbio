import os as os
import csv as csv           #use Python's built in csv handler
import sys
import requests as request  #use the requests Python library
from Bio import *           #use the BioPython Python library
from Bio.PDB import *

#this service will be used to support Restful API requests to the PDB server.
#https://www.rcsb.org/pdb/software/rest.do
#linking files together

#current PDB Ids are here https://www.rcsb.org/pdb/rest/getCurrent

# os.chdir('/Users/curtisma/bioresearch/compbio')

class PDBRestServicehelper(object):
    def describeMolRetriever (self,structure):
        baseurl='https://www.rcsb.org/pdb/rest/describeMol'
        payload={'structureId':structure}
        apiRequest = request.get(baseurl,payload)
        if apiRequest.status_code == 200:
            print (apiRequest.status_code)
            print (apiRequest.content)
            return apiRequest
        else:
            print "200 status was not returned"

    def statusSequenceRetriever (self,structure):
        baseurl='https://www.rcsb.org/pdb/rest/getStatusSequence'
        payload={'structureId':structure}
        apiRequest = request.get(baseurl,payload)
        if apiRequest.status_code == 200:
            print (apiRequest.status_code)
            print (apiRequest.content)
            return apiRequest
        else:
            print "200 status was not returned"
