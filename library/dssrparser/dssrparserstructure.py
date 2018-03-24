import os as os #import the python os library
import csv as csv #import the python csv library
import re #provides regular expression matching
import json

class DssrParserjson(object):
    # we need to add 0's back in... from dssr to hbplus comparison
    def dssrtohbplusstringcleaner(self,dssnbtobecleaned):
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

    def dssrhelixParser(self,data,j):
      print "This will take a dssr json output and parse out the helixes"
      helixtoappend=[]
      k=0
      while k < len(data[j]):
        h=0
        # hf=0
        while h < len(data[j][k]["pairs"]):
        #   print len(data[j][k]["pairs"][1])
           # while h might need to loop depending upon number of helices because it's evaluated in stacks of two.
        #   left, right = str(round(hf,0)).split(".")
          hindex = data[j][k]["index"]
          helixform = data[j][k]["helix_form"]
          hiindex = data[j][k]["pairs"][h]["index"]
          hnt1 = data[j][k]["pairs"][h]["nt1"]
          hnt2 = data[j][k]["pairs"][h]["nt2"]
          strand1=hnt1[0]
          strand2=hnt2[0]
          hnt1cleaned = self.dssrtohbplusstringcleaner(hnt1)
          hnt2cleaned = self.dssrtohbplusstringcleaner(hnt2)
          dssr = data[j][k]["pairs"][h]["DSSR"]
          try:
              helixf=helixform[h]
          except:
              helixf="end"
          helixtoappend.append("%s %s %s %s %s %s %s %s" %(j,strand1,hnt1cleaned,hindex,hiindex,helixf, "nt1", dssr))
          helixtoappend.append("%s %s %s %s %s %s %s %s" %(j,strand2,hnt2cleaned,hindex,hiindex,helixf, "nt2", dssr))

        #   hf=hf+0.5
          h+=1
        k+=1
      return helixtoappend

    # this can be deprecated, will only work with the helices output. Stem is a subset of Helices
    # def dssrstemParser(self,data,j):
    #   stemtoappend=[]
    #   # this is stepping out to the loop, may want to try switching it over to for loop instead
    #   for k in range(0,len(data[j])):
    #     for h in range(0,len(data[j][k]["pairs"])):
    #       sindex=data[j][k]["index"]
    #       siindex = data[j][k]["pairs"][h]["index"]
    #       snt1 = data[j][k]["pairs"][h]["nt1"]
    #       snt2 = data[j][k]["pairs"][h]["nt2"]
    #       # print "%s %s %s %s %s" %(j,sindex, siindex, snt1,snt2)
    #       # snt1=snt1[2:]
    #       # snt2=snt2[2:]
    #       chain1=snt1[0]
    #       chain2=snt2[0]
    #       snt1cleaned = self.dssrtohbplusstringcleaner(snt1)
    #       snt2cleaned = self.dssrtohbplusstringcleaner(snt2)
    #       print "%s %s %s %s %s" %(j,chain1,snt1cleaned,sindex,siindex)
    #       stemtoappend.append("%s %s %s %s %s" %(j,chain1,snt1cleaned,sindex,siindex))
    #       stemtoappend.append("%s %s %s %s %s" %(j,chain2,snt2cleaned,sindex,siindex))
    #   return stemtoappend

    #I should start calling individual handlers here for DSSR and build them here
    def hairpindssrparser(self,data,j):
      # print data[j][0]["nts_long"]
      hairpintoappend=[]
      hairpinnt=data[j][0]["nts_long"]
      print "%s %s" %(j,hairpinnt)
      #hairpinnt is a comma separated string #"B.U13,B.C14,B.A15,B.U16,B.U17,B.A18"
      for hairpinnb in hairpinnt.split(","):
        #we are stripping B. for now lets assume RNA is always the Bstrand
        dssrnbaddedzeros = self.dssrtohbplusstringcleaner(hairpinnb)
        # print testing
        chain = hairpinnb[0]
        hairpinnb=hairpinnb[2:]
        hairpintoappend.append("%s %s %s"%(j,chain,dssrnbaddedzeros))
      return hairpintoappend
        #we need to add the 0s back in for hbplus...
