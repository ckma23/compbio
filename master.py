import os as os                                             #import the python os library
import csv as csv                                           #import the python csv library
import sys                                                  #import the sys to retrieve command line arguments
import re                                                   #provides regular expression matching
import json                                                 #import the biopython json library
import ConfigParser                                         #import the configParser library
import glob                                                 #import the glob library used for regression grabbing
import fnmatch                                              #import the fnmatch library

from pprint import pprint                                   #pprint for json
from Bio import *                                           #use the BioPython Python library
from Bio.PDB import *                                       #more specifically import the BioPython PDB library
from optparse import OptionParser as optionparser           # will need a option parsing library

from helper.biopythonpdbretriever import FileRetriever      #from the helper folder biopythonpdbretriever.py import FileRetriever class
from helper.PDBRestService import PDBRestServicehelper      #from the helper folder PDBRestService import FileRetriever class
from metadata_folder.testset_hashstore import TestsetHashstore #from the metadata_folder folder testset_hashstore import FileRetriever class

from library.culling_helper.cull_helper import CullHelper
from library.hbplusclilib.hbpluscli import hbplusclihelper
from library.hbplusclilib.hbpluscli_testset import HbplusCliTestset as HbplusCliTestset
from library.hbplusparser.hbplusparser import HbPlusProcesser as HbPlusProcesser
from library.hbplusparser.hbplusparser_testset import HbPlusProcesserTestSet as HbPlusProcesserTestSet
from library.dssrclilib.dssrcli import dssrclihelper
from library.dssrclilib.dssrcli_testset import DssrCliHelperTestset as DssrCliHelperTestset
from library.dssrparser.dssrparser import DssrParser as DssrParser
from library.dssrparser.dssrparser_testset import DssrParserTestSet
from library.hbplus_to_dssr_comparer.hbplus_to_dssr_comparer import HbPlusToDssrComparer as HbPlusToDssrComparer
from library.hbplus_to_dssr_comparer.hbplus_to_dssr_comparer_testset import HbPlusToDssrComparerTestset as HbPlusToDssrComparerTestset
from library.bondfinalcount.statistical_potential_calculator import StatisticalPotential as StatisticalPotential
from library.energy_formation.energy_formation_calculator import energyCalculator as energyCalculator
from library.bondfinalcount.bondcount import Bondcounter
from library.ftdockprepper.chain_stripping import FTDockChainStripping
from library.ftdock_middleware.ftdock_middleware import FtdockMiddleware

def directory_builder():
    os.chdir(os.path.expanduser('~/bioresearch/compbio'))
    #make the directory logs
    os.system("mkdir logs")
    #make the work in progress folder directories used
    os.system("mkdir files_wip")

def pdb_file_retriever(configuration_file,folder_path_name):
    with open(configuration_file,'rb') as csvfile:
        csvstore = csv.reader(csvfile,delimiter = ',')
        for structure in csvstore:
            print (structure)
            # note: had to strip the string as it was returning ['<string>'] for []
            FileRetriever().fileretrieving(''.join(structure),folder_path_name)

def help_output():
  print "\nWelcome to the Computatinal Biology RNA-Protein Interaction help section\n"
  print "\nThe following commands are available:\n"
  print "\ndssrcli \nfileretriever \nhbplushbcli \nhbplusvdwcli \ndssrprocessedreader \nhbplusprocessedreader\n"

#from the system take the second argument can switch this out with optparse
proinput=str(sys.argv[1])

#default help menu
if proinput == "help":
    help_output()

#build out the needed directories
elif proinput == "initialize":
    directory_builder()
elif proinput == "cullhelp":
    CullHelper().cullhelper("structures_of_interest.csv")

#function to split the pdb files into rna only and protein only since we'll need to grab the unbound/bound protein
elif proinput == "parse_testset_list":
    TestsetHashstore().testset_list_prepare("protein","test_complexes_pdb_protein_unbound.csv")
    TestsetHashstore().testset_list_prepare("rna","test_complexes_pdb_rna_unbound.csv")

#note that the way the code is written. only run once, it does not clear the directory
elif proinput == "chainstrip_testset_pdb":
    TestsetHashstore().file_cleaner_based_on_chain("files_wip/test_complexes_pdb_protein_unbound","protein")
    TestsetHashstore().file_cleaner_based_on_chain("files_wip/test_complexes_pdb_rna_unbound","rna")
    TestsetHashstore().pdbfile_native_fileprepper_rna_protein("files_wip/test_complexes_pdb","protein")
    TestsetHashstore().pdbfile_native_fileprepper_rna_protein("files_wip/test_complexes_pdb","rna")

# grab the actual pdb files from the PDB database
elif proinput == "pdbfileretriever":
    pdb_file_retriever("structures_of_interest.csv","files_wip/base_complexes_pdb")
    pdb_file_retriever("test_complexes_pdb_complex_bound.csv","files_wip/test_complexes_pdb")
    pdb_file_retriever("test_complexes_pdb_protein_unbound.csv","files_wip/test_complexes_pdb_protein_unbound")
    pdb_filer_etriever("test_complexes_pdb_rna_unbound.csv","files_wip/test_complexes_pdb_rna_unbound")

#run hbplus initially for the baseset in hydrogen bonding mode
elif proinput == "hbplushbcli":
    hbplusclihelper().hbplushbcli()

#run hbplus initially for the baseset in van der Waals mode
elif proinput == "hbplusvdwcli":
    hbplusclihelper().hbplusvdwcli()

#process the hbplus baseset hydrogne bonding and van der Waals results into the needed format
elif proinput == "hbplus_parse_baseset":
    HbPlusProcesser().hbplusprocessedreader("hb")
    HbPlusProcesser().hbplusprocessedreader("vdw")

#combine the hbplus baseset hydrogen bonding and van der Waals results into a single file
elif proinput == "hbplushbvdwcombine":
    HbPlusProcesser().hbplushbandvdwcombiner()

elif proinput == "hbplus_baseset":
    hbplusclihelper().hbplushbcli()
    hbplusclihelper().hbplusvdwcli()
    HbPlusProcesser().hbplusprocessedreader("hb")
    HbPlusProcesser().hbplusprocessedreader("vdw")
    HbPlusProcesser().hbplushbandvdwcombiner()

elif proinput == "dssrcli":
    dssrclihelper().dssrcli_test()
elif proinput == "dssrparse":
    DssrParser().dssrprocessedreader()

elif proinput == "dssr_cli_parse_combine":
    dssrclihelper().dssrcli_test()
    DssrParser().dssrprocessedreader()

elif proinput == "hbcategorizedssr":
    HbPlusToDssrComparer().hbplushbvdwtodssrcomparer()

#method to only count the bounds
elif proinput == "countbonds":
    Bondcounter().file_opener()
    # Bondcounter().total_atom_counter()
elif proinput == "statistical_potential":
    StatisticalPotential().statistical_potential()

#Running it together does not seem to be working
elif proinput == "bound_count_stat_pot_baseset":
    Bondcounter().file_opener()
    StatisticalPotential().statistical_potential()

elif proinput == "pdbsplit":
    FTDockChainStripping().pdbfile_splitter_rna_protein("protein_seperated_pdbfiles_baseset","rna_seperated_pdbfiles_baseset","base_complexes_pdb")
    FTDockChainStripping().pdbfile_splitter_rna_protein("protein_seperated_pdbfiles_testset","rna_seperated_pdbfiles_testset","test_complexes_pdb")
elif proinput == "preprocessftdock":
    FTDockChainStripping().preprocess_ftdock("rna")
    FTDockChainStripping().preprocess_ftdock("protein")
# this is deprectaed
elif proinput == "testprotein_pdb_combine":
    FTDockChainStripping().pdb_file_combine_rms_calc()

elif proinput == "copyfilestosjsucluster":
    FTDockChainStripping().file_copier_to_sjsu_cluster_testset("base")
    # FTDockChainStripping().file_copier_to_sjsu_cluster_testset("test")
    # FTDockChainStripping().file_copier_to_sjsu_cluster("rna")
    # FTDockChainStripping().file_copier_to_sjsu_cluster("protein")
elif proinput == "ftdockgen":
    FtdockMiddleware().ftdock_kicker()
elif proinput == "ftdockbuild":
    # FtdockMiddleware().ftdock_directory_cleaner("home","cma","Users","curtisma")
    # FtdockMiddleware().ftdock_directory_cleaner("Users","curtisma","home","cma")
    FtdockMiddleware().ftdock_builder(1073)

elif proinput == "hbpluscli_testset":
    HbplusCliTestset().hbpluscli("hbplus_processed_hb_files_testset","hb")
    HbplusCliTestset().hbpluscli("hbplus_processed_vdw_files_testset","vdw")

elif proinput == "hbplusprocess_testset":
    HbPlusProcesserTestSet().hbplusprocessed_file_prepper_reader("hbplus_processed_hb_files_testset","hb")
    HbPlusProcesserTestSet().hbplusprocessed_file_prepper_reader("hbplus_processed_vdw_files_testset","vdw")

elif proinput == "hbplusprocess_testset_combine":
    HbPlusProcesserTestSet().hbplushbandvdwcombiner()

elif proinput == "dssrcli_testset":
    DssrCliHelperTestset().dssrcli_revamped()
elif proinput == "dssrparse_testset":
    DssrParserTestSet().dssrprocessedreader()
elif proinput == "bondcategorizer_testset":
    HbPlusToDssrComparerTestset().hbplushbvdwtodssrcomparer()

#run pymol wrapper for the native poses here....

#native checker

elif proinput == "native_checker":
    energyCalculator().native_checker_pose_hash_retriever("<protein>")

elif proinput == "energy_calculate_testset":
    energyCalculator().energy_calculator("")
#    energyCalculator().energy_calculator("hb_only")
#    energyCalculator().energy_calculator("vdw_only")
elif proinput == "energy_calculate_testset_ranking":
    energyCalculator().pose_rankings("")
#    energyCalculator().pose_rankings("hb_only")
#    energyCalculator().pose_rankings("vdw_only")

elif proinput == "completerun":
    hbplusprocessedreader("hb")
    hbplusprocessedreader("vdw")
    hbplushbandvdwcombiner()
    DssrParser().dssrprocessedreader()
    hbplushbvdwtodssrcomparer()
    Bondcounter().bondcounter()
    Bondcounter().total_atom_counter()
    preprocess_ftdock("rna")
    preprocess_ftdock("protein")
