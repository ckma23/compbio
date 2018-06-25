import sys
import os
import pymol
import __main__

__main__.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()

# os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset'))
# protein_files=os.listdir('.')
# first_file_name="1e7k_A.pdb"
# second_file_name="1ec6_A.pdb"
# #load the files that are to be compared
# def pymol_middleware ():
#     for protein_file in protein_files:
#         pymol.cmd.do("load %s,pose" %first_file_name)
#         pymol.cmd.do("load %s,complex" %protein_file)
#         # pymol.cmd.do("align pose ,complex ")
#         pymol.cmd.do("align pose and name CA,complex and name CA")
#     # pymol.cmd.do('select ou, /complex//R//')
#     # pymol.cmd.do('select co, /pose//R//')
#     # pymol.cmd.do('rms_cur ou,co')
#         pymol.cmd.do("delete %s" %first_file_name)
#         pymol.cmd.do("delete %s" %protein_file)
#         pymol.cmd.do("delete pose")
#         pymol.cmd.do("delete complex")
#     pymol.cmd.quit()

# pymol_middleware()
class PymolMiddleware(object):
    def pymol_middleware_test(self):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip'/))
        os.system("mkdir native_poses")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset'))
        protein_files=os.listdir('.')
        protein_files=["1e7k_A.pdb"]
        #expect some string cleaning
        for protein_file in protein_files:
            protein_file_name_directory = protein_file[0:4]
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset'))
            pose_directory = os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s' %protein_file_name_directory)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s' %protein_file_name_directory))
            #sort the pose files after the move.. since they are out of order
            pose_files =sorted(os.listdir('.'))
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset'))
            pose_result_file_name="%s.pose_calculated" %protein_file[0:6]
            os.system("rm *pose_calculated")
            file = open(pose_result_file_name,"a")
            # only grab the list of files that match that protein for those poses.
            for pose_file in pose_files:
                pymol.cmd.do("load %s/%s,pose" %(pose_directory,pose_file))
                pymol.cmd.do("load %s,complex" %protein_file)
                # rms = pymol.cmd.do("align pose and name CA,complex and name CA")
                rms = pymol.cmd.align("pose////CA","complex////CA")
                # typical rms response (0.0004964197287335992, 967, 1, 0.0004964197287335992, 967, 627.0, 125)
                native_or_nonnative_value = native_nonnative_checker(rms[0])
                # if rms[0] < 10.0:
                #     native_or_nonnative = "native"
                # elif rms[0] > 10.0:
                #     native_or_nonnative = "nonnative"
                print "%s %s %s %s" %(protein_file,pose_file,rms[0],native_or_nonnative_value)
                file.write("%s %s %s %s\n" %(protein_file[0:6],pose_file,rms[0],native_or_nonnative_value))
                pymol.cmd.do("delete %s" %pose_file)
                pymol.cmd.do("delete %s" %protein_file)
                pymol.cmd.do("delete pose")
                pymol.cmd.do("delete complex")
            os.system("mv %s ~/bioresearch/compbio/files_wip/native_poses") pose_result_file_name
            #eventually need to just make a directory maker....
            pymol.cmd.quit()

    def native_nonnative_checker(rms_value):
        if rms_value < 3.0:
            native_or_nonnative = "native"
        elif rms_value > 3.0:
            native_or_nonnative = "nonnative"
        return native_or_nonnative


PymolMiddleware().pymol_middleware_test()

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
