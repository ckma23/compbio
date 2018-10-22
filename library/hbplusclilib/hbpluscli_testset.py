import os as os #import the python os library
import csv as csv #import the python csv library
import re #provides regular expression matching
import time

class HbplusCliTestset(object):
    def hbpluscli(self,folder_name,vdw_or_hb):
        # lets move into the files_wip directory
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/'))
        os.system('mkdir ~/bioresearch/compbio/files_wip/%s' %folder_name)
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes'))

        list_of_protein_complexes = os.listdir('.')
        list_of_protein_complexes = ["1e7k"]
        if vdw_or_hb == "hb":
            hbpluscommand = "./hbplus"
            hb_vdw_file_format = "hb2"
        elif vdw_or_hb == "vdw":
            hbpluscommand = "./hbplus -N"
            hb_vdw_file_format = "nb2"
        # This becomes exponential O(N)^2 not ideal

        for protein_complex in list_of_protein_complexes:
            os.system("mkdir ~/bioresearch/compbio/files_wip/%s/%s" %(folder_name,protein_complex))
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s' %protein_complex))
            #make this into a set..... maybe? to not make it O(n) and make it (O1)
            list_of_pose_complexes = os.listdir('.')
            os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/hbplus'))
            os.system('rm *%s'%hb_vdw_file_format)
            print protein_complex
            for pose in set(list_of_pose_complexes):
                lhs_pose,rhs=pose.split(".",1)
                stringprep=('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s/%s' %(protein_complex,pose))
                os.system('%s %s > ~/bioresearch/compbio/logs/hbplus%s.txt &' %(hbpluscommand,stringprep,vdw_or_hb))
                # want to implement an await solution here same as ftdock, it's memory leaking or thread leaking. Jumping to the mv command
                # this was having file movement issues as it was flying through the directory
                time.sleep(.100)
                os.system("mv %s.%s ~/bioresearch/compbio/files_wip/%s/%s" %(lhs_pose,hb_vdw_file_format,folder_name,protein_complex))
                # this was having file movement issues as it was flying through the directory
                time.sleep(.100)
                os.system("mv Complex* ~/bioresearch/compbio/files_wip/%s/%s" %(folder_name,protein_complex))
