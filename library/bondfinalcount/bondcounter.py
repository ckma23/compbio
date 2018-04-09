import os

class Bondcounter(object):

    def bondcounter(self):
        os.chdir(os.path.expanduser('~/bioresearch/bondcategorized'))
        listofcategorizedfiles = os.listdir('.')
        category_storehb=[0,0,0,0,0,0,0,0,0]
        category_storevdw=[0,0,0,0,0,0,0,0,0]
        for categorized_file in listofcategorizedfiles:
            # with open(categorized_file) as csvfile:
                # csvstore = csv.reader(csvfile, delimiter = '')
                # for line in csvstore:
            categorized_file_open = open(categorized_file,"r")
            for line in categorized_file_open:
                #check each line what CAT the CAT is line[0]
                #check each line is hb or vdw is line[8]
                cat = line[0:5]
                # this is returning "hb " the space is not matching up
                hborvdw = line[34:37].strip(' ')
                # bondcounterfinal(cat,hborvdw)
                    # print line[8]
                    # print line[0]
                if hborvdw == "hb":
                    if cat == "CAT_1":
                        category_storehb[0]+=1
                    elif cat == "CAT_2":
                        category_storehb[1]+=1
                    elif cat == "CAT_3":
                        category_storehb[2]+=1
                    elif cat == "CAT_4":
                        category_storehb[3]+=1
                    elif cat == "CAT_5":
                        category_storehb[4]+=1
                    elif cat == "CAT_6":
                        category_storehb[5]+=1
                    elif cat == "CAT_7":
                        category_storehb[6]+=1
                    elif cat == "CAT_8":
                        category_storehb[7]+=1
                    elif cat == "CAT_9":
                        category_storehb[8]+=1
                elif hborvdw  == "vdw":
                    if cat == "CAT_1":
                        category_storevdw[0]+=1
                    elif cat == "CAT_2":
                        category_storevdw[1]+=1
                    elif cat == "CAT_3":
                        category_storevdw[2]+=1
                    elif cat == "CAT_4":
                        category_storevdw[3]+=1
                    elif cat == "CAT_5":
                        category_storevdw[4]+=1
                    elif cat == "CAT_6":
                        category_storevdw[5]+=1
                    elif cat == "CAT_7":
                        category_storevdw[6]+=1
                    elif cat == "CAT_8":
                        category_storevdw[7]+=1
                    elif cat == "CAT_9":
                        category_storevdw[8]+=1
        print categorized_file
        print "HB  CAT_1:%s,CAT_2:%s,CAT_3:%s,CAT_4:%s,CAT_5:%s,CAT_6:%s,CAT_7:%s,CAT_8:%s,CAT_9:%s" %(category_storehb[0],category_storehb[1],category_storehb[2],category_storehb[3],category_storehb[4],category_storehb[5],category_storehb[6],category_storehb[7],category_storehb[8])
        print "VDW CAT_1:%s,CAT_2:%s,CAT_3:%s,CAT_4:%s,CAT_5:%s,CAT_6:%s,CAT_7:%s,CAT_8:%s,CAT_9:%s" %(category_storevdw[0],category_storevdw[1],category_storevdw[2],category_storevdw[3],category_storevdw[4],category_storevdw[5],category_storevdw[6],category_storevdw[7],category_storevdw[8])
