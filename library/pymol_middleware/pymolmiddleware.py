import sys
import os
import pymol
import __main__


__main__.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()

class TestsetHashstore:
    def __init__(self):
        self.hashy = {
        "1WSU_A_E":{"protein":"1LVA_A","rna":"1MFK_A"},
        "2PJP_A_B":{"protein":"2PJP_A","rna":"1MFK_A"},
        "1LNG_A_B":{"protein":"1LNG_A","rna":"1Z43_A"},
        "1E7K_A_C":{"protein":"2JNB_A","rna":"1E7K_C"},
        "1WPU_A_C":{"protein":"1WPV_A","rna":"1WPU_C"},
        "2QUX_A_C":{"protein":"2QUD_A","rna":"2QUX_C"},
        "2JEA_A_C":{"protein":"2JE6_A","rna":"2JEA_C"},
        "2FMT_A_C":{"protein":"1FMT_A","rna":"3CW5_A"}, 
        "1MFQ_C_A":{"protein":"1QB2_B","rna":"1L9A_B"},
        "1U0B_B_A":{"protein":"1L17_A","rna":"1B23_R"},
        "1EC6_A_D":{"protein":"1DTJ_A","rna":"1EC6_D"},
        "1HC8_A_C":{"protein":"1FOY_A","rna":"1HC8_C"},
        "1JBR_B_D":{"protein":"1AQZ_A","rna":"1JBR_D"},
        "1KOG_A_I":{"protein":"1EVL_A","rna":"1KOG_I"},
        "1M8W_A_C":{"protein":"1M8Z_A","rna":"1M8W_C"},
        "1F7U_A_B":{"protein":"1BS2_A","rna":"1F7U_B"},
        "1K8W_A_B":{"protein":"1R3F_A","rna":"1K8W_B"},
        "1N78_A_C":{"protein":"1J09_A","rna":"1N78_C"},
        "1U63_A_B":{"protein":"1I2A_A","rna":"1U63_B"},
        "2BTE_A_B":{"protein":"1H3N_A","rna":"2BTE_B"},
        "2HW8_A_B":{"protein":"1AD2_A","rna":"2HW8_B"},
        }
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

    def which_chains(self,identifier):
        testset_retrieve_hash = TestsetHashstore().hashy
        identifier_uppercase = identifier.upper()
        for native in testset_retrieve_hash.keys():
            if identifier_uppercase == native[0:4]:
                native_protein_chain = native[5]
                native_rna_chain = native[7]
                pose_protein_chain = testset_retrieve_hash[native]["protein"][5]
                pose_rna_chain = testset_retrieve_hash[native]["rna"][5]
        return native_rna_chain, native_protein_chain, pose_rna_chain, pose_protein_chain

    def pymol_middleware_test(self):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip'))
        os.system("mkdir native_poses_testset")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/cleaned_pdbfiles_testset_bound'))
        # os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/combined_pdbfiles_testset'))
        # os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset'))
        protein_files = os.listdir('.')
        # protein_files = ["1e7k_A.pdb"]

        #expect some string cleaning
        for protein_file in protein_files:
            protein_file_name_directory = protein_file[0:4]
            native_rna_chain_use, native_protein_chain_use, pose_rna_chain_use, pose_protein_chain_use = self.which_chains(protein_file_name_directory)

            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/cleaned_pdbfiles_testset_bound'))
            # os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset'))
            pose_directory = os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s' %protein_file_name_directory)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s' %protein_file_name_directory))
            #sort the pose files after the move.. since they are out of order
            pose_files = sorted(os.listdir('.'))
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/cleaned_pdbfiles_testset_bound'))
            # os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_testset'))
            pose_result_file_name="%s" %protein_file[0:4]
            #remove the file previously created since we are opening in append mode
            os.system("rm *pose_calculated")
            file = open(pose_result_file_name,"a")
            # from analysis on 10/09/2018
            # align the native structure and the complex by alpha carbons

            # only grab the list of files that match that protein for those poses.
            for pose_file in pose_files:
                pymol.cmd.do("load %s/%s,pose" %(pose_directory,pose_file))
                pymol.cmd.do("load %s,complexes" %protein_file)
                # pymol.cmd.do("align pose and name P,complexes and name P")

                # pymol.cmd.do("align pose and name CA and chain %s,complexes and name CA and chain %s") %(pose_protein_chain_use,native_protein_chain_use)
                # as of 10-09-18 since protein is static and rna is mobile
                # align the pose and complex by alpha carbon
                # note that poses are the 50373 and complexes are the "native" need to change syntax
                # pymol.cmd.do("align pose and name CA and chain %s,complexes and name CA and chain %s" %(pose_protein_chain_use,native_protein_chain_use))
                # no need to align by chain, lets jsut align by the alpha carbons
                pymol.cmd.do("align pose and name CA,complexes and name CA")
                # grab all the RNA atoms could also align by Phosphate probably not
                # pymol.cmd.do('select complex_atoms, /complex//%s//P' %native_rna_chain_use)
                pymol.cmd.do('select complex_atoms, /complex//%s//' %native_rna_chain_use)
                # now as of 10/10/2018 found out that if the chain identifiers don't match pymol will throw an error ExecutiveRMS-Error: No atoms selected.
                # https://pymolwiki.org/index.php/Fit
                # therefore we actually need the rna chains to have the same identifier
                # pymol.cmd.do('select pose_atoms, /pose//%s//P' %pose_rna_chain_use)
                pymol.cmd.do('select pose_atoms, /pose//%s//' %pose_rna_chain_use)
                pymol.cmd.do('alter pose_atoms, chain="%s"' %native_rna_chain_use)
                rms = pymol.cmd.rms_cur("pose_atoms","complex_atoms")
                # rms = pymol.cmd.do('rms_cur ou,co')
                # rms = pymol.cmd.rms_cur(ou,co)
                # rms = pymol.cmd.do("rms_cur pose////CA,complexes////CA")
                # print rms
                #https://pymol.org/dokuwiki/doku.php?id=command:rms_cur
                #“rms_cur” computes the RMS difference between two atom selections without performing any fitting.
                # rms = pymol.cmd.rms_cur("pose////CA","complexes////CA")
                print rms
                # rms = pymol.cmd.align("pose////CA","complexes////CA")
                # typical rms response (0.0004964197287335992, 967, 1, 0.0004964197287335992, 967, 627.0, 125)
                native_or_nonnative_value = self.native_nonnative_checker(rms)
                # if rms[0] < 10.0:
                #     native_or_nonnative = "native"
                # elif rms[0] > 10.0:
                #     native_or_nonnative = "nonnative"
                print "%s %s %s %s" %(protein_file,pose_file,rms,native_or_nonnative_value)
                pose_file_to_store,throwaway = pose_file.split(".",1)
                pose_file_to_store = pose_file_to_store.strip("g")
                file.write("%s,%s,%s,%s\n" %(protein_file[0:4],pose_file_to_store,rms,native_or_nonnative_value))
                pymol.cmd.do("delete %s" %pose_file)
                pymol.cmd.do("delete %s" %protein_file)
                pymol.cmd.do("delete pose")
                pymol.cmd.do("delete complexes")
                pymol.cmd.do("delete pose_atoms")
                pymol.cmd.do("delete complex_atoms")
            #move the file to where it will eventually be used
            os.system("mv %s ~/bioresearch/compbio/files_wip/native_poses_testset" %pose_result_file_name)
            #eventually need to just make a directory maker....
        pymol.cmd.quit()

    def native_nonnative_checker(self,rms_value):
        if rms_value < 10.0:
            native_or_nonnative = "native"
        elif rms_value > 10.0:
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
