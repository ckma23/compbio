import os as os #import the python os library
import csv as csv #import the python csv library
import re #provides regular expression matching
import time

class HbplusCliTestSet(object):

    def hbplushbcli(self,folder):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip'))
        os.system('mkdir %s' %folder)
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes'))
        list_of_protein_complexes = os.listdir('.')
        # This becomes exponential O(N)^2 not ideal
        for protein_complex in list_of_protein_complexes:
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s' %protein_complex))
            #make this into a set..... maybe? to not make it O(n) and make it (O1)
            list_of_pose_complexes = os.listdir('.')
            os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/hbplus'))
            os.system('rm output')
            # for pose in set(list_of_pose_complexes):
            #     stringprep=('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s/%s' %(protein_complex,pose))
            #     os.system('./hbplus %s >> output &' %stringprep)
            #     time.sleep(.050)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s' %folder))
            # os.system("mkdir %s" %protein_complex)
            # for numb in range(1,55):
            #     os.system("mv Complex_%i* ~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s" %(numb,testset_ftdocked_pdbfile[0:4]))
            # os.system("mv Complex*.pdb ~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s" %testset_ftdocked_pdbfile[0:4])
            # os.system('mv %s.hb2 ~/bioresearch/compbio/files_wip/%s/%s' %(pose,folder,protein_complex))

    def hbpluscli(self,folder,vdw_or_hb):
        # lets move into the files_wip directory
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip'))
        os.system('mkdir %s' %folder)
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes'))
        list_of_protein_complexes = os.listdir('.')
        # lets set this hbplus command here so we don't bring it into the for loop increasing big O complexity
        if vdw_or_hb == "hb":
            hbpluscommand = "./hbplus"
        elif vdw_or_hb == "vdw":
            hbpluscommand = "./hbplus -N"
        # This becomes exponential O(N)^2 not ideal
        for protein_complex in list_of_protein_complexes:
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s' %protein_complex))
            #make this into a set..... maybe? to not make it O(n) and make it (O1)
            list_of_pose_complexes = os.listdir('.')
            os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/hbplus'))
            os.system('rm output')
            for pose in set(list_of_pose_complexes):
                stringprep=('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s/%s' %(protein_complex,pose))
                os.system('./hbplus %s >> output &' %stringprep)
                time.sleep(.050)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/%s' %folder))
            # os.system("mkdir %s" %protein_complex)
            # for numb in range(1,55):
            #     os.system("mv Complex_%i* ~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s" %(numb,testset_ftdocked_pdbfile[0:4]))
            # os.system("mv Complex*.pdb ~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s" %testset_ftdocked_pdbfile[0:4])
            # os.system('mv %s.hb2 ~/bioresearch/compbio/files_wip/%s/%s' %(pose,folder,protein_complex))

# HbplusCliTestSet().hbplushbcli("hbplus_processed_hb_files_testset")
# HbplusCliTestSet().hbplushbcli("hbplus_processed_hb_files_testset","hb")
# HbplusCliTestSet().hbplushbcli("hbplus_processed_hb_files_testset","vdw")
# HbplusCliTestSet().hbplusvdwcli("hbplus_processed_hb_files_testset")
