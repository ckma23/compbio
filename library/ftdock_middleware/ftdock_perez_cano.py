import os
import time

from library.energy_formation.energy_formation_calculator import energyCalculator as energyCalculator

class FtdockPerezCano(object):
    def perez_cano_rankings(self):
        os.system("mkdir ~/bioresearch/compbio/files_wip/perez_cano_rankings")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockresults'))
        pere_cano_matrix = ""
        ftdocked_pdb_files = os.listdir('.')
        for ftdockedpdbfile in ftdocked_pdb_files:
            os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/ftdock-2-dev2/progs-2.0.3'))
            os.system("./rpscore -matrix -in ~/bioresearch/compbio/files_wip/ftdockresults/%s -out ~/bioresearch/compbio/files_wip/perez_cano_rankings/%s.perez_cano_scored > ~/bioresearch/compbio/logs/rpscore &" %(ftdockedpdbfile,ftdockedpdbfile))
            # os.system("./rpscore -matrix %s -in ~/bioresearch/compbio/files_wip/ftdockresults/%s -out ~/bioresearch/compbio/files_wip/perez_cano_rankings/%s.perez_cano_scored > ~/bioresearch/compbio/logs/rpscore &" %(pere_cano_matrix,ftdockedpdbfile,ftdockedpdbfile))

    def matching_complex_to_rpscore(self):
        native_perez_cano_pose_hash, native_pose_file = energyCalculator.native_checker_pose_hash_retriever(testset_protein)
        for line in perez_cano_protein:
            #Complex numbers by the
            # Data
            #Type       ID    prvID    SCscore        RPscore         Coordinates            Angles
            #G_DATA      1    25187         65          0.000         1   -1  -30        76  99 108

            #ID is the Complex Number
            #RP Score is what we need to rank on

            native_pose_hash[key]=line[0]
            #line[0] should be the protein name
            #line[1] is the Complex Number
            #line[2] is the RMS
            #line[3] is the native or nonnative
            native_pose_hash[native_pose_file][line[1]] = {}
            #originally looked this
            #native_pose_hash[native_pose_file][line[1]] = {"RMS":,"nativeness":line[3]}
            native_pose_hash[native_pose_file][line[1]]["perez_cano_score"] = line[4]
