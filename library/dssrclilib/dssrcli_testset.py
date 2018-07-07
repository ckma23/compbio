import os
import time

class DssrCliHelperTestset(object):
    def dssrcli(self):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/'))
        os.system('mkdir DSSRprocessed_testset')
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes'))
        list_of_testcomplexes_protein = os.listdir('.')
        for testprotein in list_of_testcomplexes_protein:
            os.system('mkdir ~/bioresearch/compbio/files_wip/DSSRprocessed_testset/%s' %testprotein)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s' %testprotein))
            list_of_testcomplexes_pose = os.listdir('.')
            for testcomplex_pose in list_of_testcomplexes_pose:
                stringprep=('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s/%s' %(testprotein,testcomplex_pose))
                os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/'))
                lhs,rhs=testcomplex_pose.split(".",1)
                os.system('./x3dna-dssr input=%s -o=%s.json -json > ~/bioresearch/compbio/logs/DSSR.log' %(stringprep,lhs))
                # can remove all the auxiliary dssr-complimentary files since we are using the json one anyways.
                os.system('rm dssr-*')
                # this was having file movement issues as it was flying through the directory
                time.sleep(.050)
                os.system('mv *.json ~/bioresearch/compbio/files_wip/DSSRprocessed_testset/%s' %(testprotein))
                # this was having file movement issues as it was flying through the directory
                time.sleep(.050)
