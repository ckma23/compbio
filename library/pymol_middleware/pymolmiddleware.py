import sys
import os
import pymol

import __main__

__main__.pymol_argv = ['pymol', '-c']
pymol.finish_launching()

os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset'))
protein_files=os.listdir('.')
first_file_name="1e7k_A.pdb"
second_file_name="1ec6_A.pdb"
#load the files that are to be compared
def pymol_middleware ():
    for protein_file in protein_files:
        pymol.cmd.do("load %s,pose" %first_file_name)
        pymol.cmd.do("load %s,complex" %protein_file)
        # pymol.cmd.do("align pose ,complex ")
        pymol.cmd.do("align pose and name CA,complex and name CA")
    # pymol.cmd.do('select ou, /complex//R//')
    # pymol.cmd.do('select co, /pose//R//')
    # pymol.cmd.do('rms_cur ou,co')
        pymol.cmd.do("delete %s" %first_file_name)
        pymol.cmd.do("delete %s" %protein_file)
        pymol.cmd.do("delete pose")
        pymol.cmd.do("delete complex")
    pymol.cmd.quit()

pymol_middleware()

def pymol_middleware_test ():
    os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset'))
    protein_files=os.listdir('.')
    #expect some string cleaning
    for protein_file in protein_files:
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset'))
        pose_directory = os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset')
        pose_files = os.listdir('.')
        # only grab the list of files that match that protein for those poses.
        for pose_file in pose_files:
            pymol.cmd.do("load %s/%s,pose" %pose_directory,pose_file)
            pymol.cmd.do("load %s,complex" %protein_file)
            pymol.cmd.do("align pose and name CA,complex and name CA")
    # pymol.cmd.do('select ou, /complex//R//')
    # pymol.cmd.do('select co, /pose//R//')
    # pymol.cmd.do('rms_cur ou,co')
            pymol.cmd.do("delete %s" %pose_file)
            pymol.cmd.do("delete %s" %protein_file)
            pymol.cmd.do("delete pose")
            pymol.cmd.do("delete complex")


# import os, fnmatch
#
# listOfFiles = os.listdir('.')
# pattern = "*.py"
# for entry in listOfFiles:
#     if fnmatch.fnmatch(entry, pattern):
#             print (entry)



# import pymol
#
# class PymolManipulator(object):
#     def execute_pymol(self):
#         print "hi"
#
# PymolManipulator.execute_pymol()


# import __main__
# __main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
#
# import sys, time, os
# import pymol
#
# pymol.finish_launching()


#!/usr/bin/python2.6 -i

# import sys, os
#
# # autocompletion
# import readline
# import rlcompleter
# readline.parse_and_bind('tab: complete')
#
# # pymol environment
# moddir='/opt/pymol-svn/modules'
# sys.path.insert(0, moddir)
# os.environ['PYMOL_PATH'] = os.path.join(moddir, 'pymol/pymol_path')
#
# # pymol launching: quiet (-q), without GUI (-c) and with arguments from command line
# import pymol
# pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
# pymol.finish_launching()
# cmd = pymol.cmd
#
#
# os.system("pymol -cqr")
