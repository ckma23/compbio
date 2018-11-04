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
                time.sleep(.050)
                os.system('mv *.json ~/bioresearch/compbio/files_wip/DSSRprocessed_testset/%s' %(testprotein))
                time.sleep(.050)

    def dssrcli_revamped(self):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/'))
        os.system('mkdir DSSRprocessed_testset')
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes'))
        list_of_testcomplexes_protein = os.listdir('.')
        # list_of_testcomplexes_protein = ["1ec6"]
        for testprotein in list_of_testcomplexes_protein:
            os.system('mkdir ~/bioresearch/compbio/files_wip/DSSRprocessed_testset/%s' %testprotein)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s' %testprotein))
            list_of_testcomplexes_pose = os.listdir('.')
            for testcomplex_pose in list_of_testcomplexes_pose:
                stringprep=('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s/%s' %(testprotein,testcomplex_pose))
                os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/'))
                lhs,rhs=testcomplex_pose.split(".",1)
                output_file_string_prep = ('~/bioresearch/compbio/files_wip/DSSRprocessed_testset/%s/%s.json' %(testprotein,lhs))
                os.system('./x3dna-dssr input=%s --auxfile=no -o=%s -json > ~/bioresearch/compbio/logs/DSSR.log' %(stringprep,output_file_string_prep))
