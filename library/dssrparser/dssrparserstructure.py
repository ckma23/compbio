import os as os #import the python os library
import csv as csv #import the python csv library
import re #provides regular expression matching
import json

class DssrParserjson(object):
    # we need to add 0's back in... from dssr to hbplus comparison
    def dssrtohbplusstringcleaner(self,dssnbtobecleaned):
    #   print dssnbtobecleaned
      rhs = dssnbtobecleaned[3:]
    #   print rhs
    #   print len(rhs)
      lhs = dssnbtobecleaned[2]
      if len(rhs) == 1:
        rhs="000"+rhs
      elif len(rhs) == 2:
        rhs="00"+rhs
      elif len(rhs) == 3:
        # rhs.join("0",rhs)
        rhs="0"+rhs
        # print rhs
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
          dssr_base_pair_type = data[j][k]["pairs"][h]["DSSR"]
          dssr_base_pair_name = data[j][k]["pairs"][h]["name"]
          # NEED TO ADD THE NAME for wobble, ETC!
          try:
              helixf=helixform[h]
          except:
              helixf="end"
          helixtoappend.append("%s,%s,%s,%s,%s,%s,%s,%s,%s" %(j,strand1,hnt1cleaned,hindex,hiindex,helixf, "nt1", dssr_base_pair_type, dssr_base_pair_name))
          helixtoappend.append("%s,%s,%s,%s,%s,%s,%s,%s,%s" %(j,strand2,hnt2cleaned,hindex,hiindex,helixf, "nt2", dssr_base_pair_type, dssr_base_pair_name))

        #   hf=hf+0.5
          h+=1
        k+=1
      return helixtoappend

    def dssrbulgeParser(self,data,j):
        bulgetoappend=[]
        k=0
        while k < len(data[j]):
            bulges_data=data[j][k]["nts_long"]
            for bulges_residue in bulges_data.split(","):
                bulge_nb_cleaned = self.dssrtohbplusstringcleaner(bulges_residue)
                chain = bulges_residue[0]
                bulges_residue=bulges_residue[2:]
                bulgetoappend.append("%s,%s,%s" %(j,chain,bulge_nb_cleaned))
            k+=1
        return bulgetoappend


    def dssrhairpinParser(self,data,j):
      # need to check how many hairpins there are... it's only going through one right now..
      # print data[j][0]["nts_long"]
        hairpintoappend=[]
        k = 0
        while k <len(data[j]):
            hairpinnt=data[j][k]["nts_long"]
      #hairpinnt is a comma separated string #"B.U13,B.C14,B.A15,B.U16,B.U17,B.A18"
            for hairpinnb in hairpinnt.split(","):
        #we are stripping B. for now lets assume RNA is always the Bstrand
                dssrnbaddedzeros = self.dssrtohbplusstringcleaner(hairpinnb)
        # print testing
                chain = hairpinnb[0]
                hairpinnb=hairpinnb[2:]
                hairpintoappend.append("%s,%s,%s"%(j,chain,dssrnbaddedzeros))
            k+=1
        return hairpintoappend
        #we need to add the 0s back in for hbplus...

    def dssrjunctionsParser(self,data,j):
        junctiontoappend=[]
        k=0
        while k <len(data[j]):
            junctionnt=data[j][k]["nts_long"]
            for junctionnb in junctionnt.split(","):
                dssjunctionsaddedzeros = self.dssrtohbplusstringcleaner(junctionnb)
                chain = junctionnb[0]
                junctionnb = junctionnb[2:]
                junctiontoappend.append("%s,%s,%s"%(j,chain,dssjunctionsaddedzeros))
            k+=1
        return junctiontoappend

    def dssriloopsParser(self,data,j):
        ilooptoappend=[]
        k=0
        while k <len(data[j]):
            iloopnt=data[j][k]["nts_long"]
      #hairpinnt is a comma separated string #"B.U13,B.C14,B.A15,B.U16,B.U17,B.A18"
            for iloopnb in iloopnt.split(","):
        #we are stripping B. for now lets assume RNA is always the Bstrand
                dssrnbaddedzeros = self.dssrtohbplusstringcleaner(iloopnb)
        # print testing
                chain = iloopnb[0]
                iloopnb=iloopnb[2:]
                ilooptoappend.append("%s,%s,%s"%(j,chain,dssrnbaddedzeros))
            k+=1
        return ilooptoappend
