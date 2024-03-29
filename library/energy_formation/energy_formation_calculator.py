import os as os                                             #import the python os library
import csv as csv                                           #import the python csv library
import sys                                                  #import the sys to retrieve command line arguments
import re                                                   #provides regular expression matching
import json                                                 #import the biopython json library
import ConfigParser                                         #import the configParser library
import glob                                                 #import the glob library used for regression grabbing
import fnmatch
#need the Statistical Potential Hash
from library.bondfinalcount.statistical_potential_calculator import StatisticalPotential as StatisticalPotential
from library.energy_formation.energy_balancer import Energyknockdown as Energyknockdown

class EnergyCalculator(object):
    #given the protein name, load in the complex # and the RMS and the native or nonnative column
    def native_checker_pose_hash_retriever(self,protein_name_string):
        # move into the native_poses directory which has all the native poses
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/native_poses_testset'))
        # native_pose_files = os.listdir('.')
        # for native_pose_file in native_pose_files:
        native_pose_file = protein_name_string
        # fileregex_match = protein_name_string + "*"
        # if fnmatch.fnmatch(native_pose_file,fileregex_match):
        #     native_posefile = native_pose_file
        native_pose_hash = {}
        native_pose_hash[native_pose_file] = {}
        with open( native_pose_file,'rb') as posefile:
            pose_file  = csv.reader(posefile)
            for line in pose_file:
                    #line[0] should be the protein name
                    #line[1] is the Complex Number
                    #line[2] is the RMS
                    #line[3] is the native or nonnative
                native_pose_hash[native_pose_file][line[1]] = {}
                native_pose_hash[native_pose_file][line[1]] = {"RMS":line[2],"nativeness":line[3]}
        return native_pose_hash , native_pose_file

        print "Checking if the structure is native"

    def retrieve_statistical_potential(self):
        #create the file_path to the statistical potential json hash
        stat_potential_file_path = os.path.join(os.path.expanduser('~/bioresearch/compbio/metadata_folder/'),"statistical_potential.json")
        #open the statistical_potentials which are stored as json
        stat_pot_file = open(stat_potential_file_path)
        #read method on the file converts the file object into the string
        stat_pot_file_string = stat_pot_file.read()
        #since it is a string now can use json loads method to convert this into a Python hash
        stat_pot_hash = json.loads(stat_pot_file_string)
        print stat_pot_hash
        #return the statistical potential hash now to the caller.
        return stat_pot_hash

    def energy_calculator(self,hb_only_or_hb_and_vdw):
        os.system("mkdir ~/bioresearch/compbio/files_wip/energy_calculations_testset%s" %hb_only_or_hb_and_vdw)
        # statistical_potential_hash = StatisticalPotential().statistical_potential()
        #grab the statistical_potential_hash
        statistical_potential_hash = self.retrieve_statistical_potential()
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/bondcategorized_testset'))
        # for each of the proteins in the testset
        testset_proteins = os.listdir('.')
        for testset_protein in testset_proteins:
            energy_knockdown_hash = Energyknockdown().energyknockdownretriever()
            native_pose_hash, native_pose_file = self.native_checker_pose_hash_retriever(testset_protein)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/energy_calculations_testset%s' %hb_only_or_hb_and_vdw))
            os.system("rm %s" %testset_protein)
            #create the energy_file where this is needed to be populated
            energy_file = open(testset_protein,"a")
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/bondcategorized_testset/%s' %testset_protein))
            testset_protein_poses = os.listdir('.')
            for testset_protein_pose_number in testset_protein_poses:
                hb_energy = 0
                vdw_energy = 0
                total_energy = 0
                # testset_protein_categorized_file = open(testset_protein_pose_number)
                with open( testset_protein_pose_number,'rb') as testset_protein_categorized_file_meta:
                    testset_protein_categorized_file = csv.reader(testset_protein_categorized_file_meta)

                    for line in testset_protein_categorized_file:
                        # print line
                        # use strip() as not all the columns have the same amount of spaces.
                        # category = line[0:5].strip(' ')
                        # hborvdw = line[34:37].strip(' ')
                        # nucleotide_base = line[6].strip(' ')
                        # amino_acid = line[25:29].strip(' ')
                        category = line[0]
                        hborvdw = line[7]
                        nucleotide_base = line[1][0]
                        amino_acid = line[5]

                        # print category
                        # print hborvdw
                        # print nucleotide_base
                        # print amino_acid

                        bond_energy = 0
                        try:
                            bond_energy = statistical_potential_hash[hborvdw][category][nucleotide_base][amino_acid]["statistical_potential"]
                            print bond_energy
                            bond_energy = energy_knockdown_hash[hborvdw][category]*bond_energy
                            print energy_knockdown_hash[hborvdw][category]
                            print bond_energy
                            # print energy_knockdown_hash[hborvdw][category]
                            # print type(energy_knockdown_hash)
                            # print energy_knockdown_hash
                        except Exception as e:
                            print "There was an error: %s" %e
                        if hborvdw == "hb":
                            hb_energy = hb_energy + bond_energy
                        elif hborvdw == "vdw":
                            vdw_energy = vdw_energy + bond_energy
                        print hb_energy
                        print vdw_energy
                    if hb_only_or_hb_and_vdw == "hb_only":
                        vdw_energy = 0;
                    elif hb_only_or_hb_and_vdw == "vdw_only":
                        hb_energy = 0;
                    total_energy = hb_energy + vdw_energy

                    testset_protein_pose_number_to_store,throwaway = testset_protein_pose_number.split('.',1)
                    testset_protein_pose_number_to_store = testset_protein_pose_number_to_store.strip('g')
                    try:
                        pose_RMS = native_pose_hash[native_pose_file][testset_protein_pose_number_to_store]['RMS']
                        pose_nativeness = native_pose_hash[native_pose_file][testset_protein_pose_number_to_store]['nativeness']
                    except Exception as e:
                        print "There was an error in retrieving the pose RMS: %s" %e
                    line = "%s,%s,%s,%s,%s,%s,%s,%s\n" %(testset_protein,testset_protein_pose_number_to_store,hb_energy,vdw_energy,total_energy,pose_RMS,pose_nativeness,native_pose_file)
                    print line
                    energy_file.write(line)
                energy_file.close

    def pose_rankings(self,hb_only_or_hb_and_vdw):
        os.system("mkdir ~/bioresearch/compbio/files_wip/best_rankings")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/best_rankings'))
        os.system("rm %srankings" %hb_only_or_hb_and_vdw)
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/energy_calculations_testset%s' %hb_only_or_hb_and_vdw))
        testset_proteins = os.listdir('.')
        for testset_protein in testset_proteins:
            print testset_protein
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/energy_calculations_testset%s' %hb_only_or_hb_and_vdw))
            native_pose_ranking_hash = {}
            native_pose_ranking_hash[testset_protein] = {}
            with open(testset_protein,'rb') as energy_file:
                energy_file_to_be_ranked  = csv.reader(energy_file)
                for line in energy_file_to_be_ranked:
                    native_pose_ranking_hash[testset_protein][line[1]] = {"energy":line[4],"RMS":line[5],"nativeness":line[6]}
            # for complex_pose,values in sorted(native_pose_ranking_hash[testset_protein].items()):
            #     for energy in sorted(values["energy"]):
            # print native_pose_ranking_hash[testset_protein].keys
            sorted_complex_list = sorted(native_pose_ranking_hash[testset_protein].keys(), key=lambda x: (float(native_pose_ranking_hash[testset_protein][x]["energy"])))
            print sorted_complex_list
            native_count = 0
            for complex_num in sorted_complex_list:
                print complex_num
                print native_pose_ranking_hash[testset_protein][complex_num]["energy"]
                value_index = sorted_complex_list.index(complex_num)
                print value_index
                if native_pose_ranking_hash[testset_protein][complex_num]["nativeness"] == "native":
                    print "These structures are native"
                    native_count+=1
                    print native_count
                    # the first native structure should either be the maximum or minumum energy.
            for complex_num in sorted_complex_list:
                if native_pose_ranking_hash[testset_protein][complex_num]["nativeness"] == "native":
                    # add 1 because the array starts at 0
                    best_rank = float(sorted_complex_list.index(complex_num)+1)/len(sorted_complex_list)
                    best_rank = best_rank * 100
                    break
            if native_count == 0:
                best_rank = 0
            print "%s,%s,%s" %(testset_protein,native_count,best_rank)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/best_rankings'))
            file_name_string_prepper = "%srankings" %hb_only_or_hb_and_vdw
            rank_file = open(file_name_string_prepper,"a")
            rank_file.write("%s,%s,%s\n" %(testset_protein,native_count,best_rank))
            rank_file.close
