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

from helper.biopythonpdbretriever import FileRetriever         #from the helper folder biopythonpdbretriever.py import FileRetriever class
from helper.PDBRestService import PDBRestServicehelper         #from the helper folder PDBRestService import FileRetriever class
from metadata_folder.testset_hashstore import TestsetHashstore #from the metadata_folder folder testset_hashstore import FileRetriever class

from library.culling_helper.cull_helper import CullHelper
from library.hbplusclilib.hbpluscli import HbplusCliHelper
from library.hbplusclilib.hbpluscli_testset import HbplusCliTestset as HbplusCliTestset
from library.hbplusparser.hbplusparser import HbPlusProcesser as HbPlusProcesser
from library.hbplusparser.hbplusparser_testset import HbPlusProcesserTestSet as HbPlusProcesserTestSet
from library.dssrclilib.dssrcli import DssrCliHelper
from library.dssrclilib.dssrcli_testset import DssrCliHelperTestset as DssrCliHelperTestset
from library.dssrparser.dssrparser import DssrParser as DssrParser
from library.dssrparser.dssrparser_testset import DssrParserTestSet
from library.hbplus_to_dssr_comparer.hbplus_to_dssr_comparer import HbPlusToDssrComparer as HbPlusToDssrComparer
from library.hbplus_to_dssr_comparer.hbplus_to_dssr_comparer_testset import HbPlusToDssrComparerTestset as HbPlusToDssrComparerTestset
from library.bondfinalcount.statistical_potential_calculator import StatisticalPotential as StatisticalPotential
from library.energy_formation.energy_formation_calculator import EnergyCalculator as EnergyCalculator
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

#from the system take the second argument from command line. Can switch this out with optparse
prog_input=str(sys.argv[1])

#default help menu
if prog_input == "help":
    help_output()

#build out the needed directories
elif prog_input == "initialize":
    directory_builder()

# the cullhelper takes the PIECES file output and converts it into the "structures_of_interest.csv" learning set
elif prog_input == "cullhelp":
    CullHelper().cull_helper("structures_of_interest.csv")

#function to split the pdb files into rna only and protein only since we'll need to grab the unbound/bound protein
elif prog_input == "parse_testset_list":
    TestsetHashstore().testset_list_prepare("protein","test_complexes_pdb_protein_unbound.csv")
    TestsetHashstore().testset_list_prepare("rna","test_complexes_pdb_rna_unbound.csv")

#note that the way the code is written. only run once, it does not clear the directory
elif prog_input == "chainstrip_testset_pdb":
    TestsetHashstore().file_cleaner_based_on_chain("files_wip/test_complexes_pdb_protein_unbound","protein")
    TestsetHashstore().file_cleaner_based_on_chain("files_wip/test_complexes_pdb_rna_unbound","rna")
    TestsetHashstore().pdbfile_native_fileprepper_rna_protein("files_wip/test_complexes_pdb","protein")
    TestsetHashstore().pdbfile_native_fileprepper_rna_protein("files_wip/test_complexes_pdb","rna")

# grab the actual pdb files from the PDB database
elif prog_input == "pdbfileretriever":
    pdb_file_retriever("structures_of_interest.csv","files_wip/base_complexes_pdb")
    pdb_file_retriever("test_complexes_pdb_complex_bound.csv","files_wip/test_complexes_pdb")
    pdb_file_retriever("test_complexes_pdb_protein_unbound.csv","files_wip/test_complexes_pdb_protein_unbound")
    pdb_filer_etriever("test_complexes_pdb_rna_unbound.csv","files_wip/test_complexes_pdb_rna_unbound")

#run hbplus initially for the baseset in hydrogen bonding mode
elif prog_input == "hbplus_hb_cli_baseset":
    HbPlusCliHelper().hbplus_hb_cli()

#run hbplus initially for the baseset in van der Waals mode
elif prog_input == "hbplus_vdw_cli_baseset":
    HbPlusCliHelper().hbplus_vdw_cli()

#process the hbplus baseset hydrogne bonding and van der Waals results into the needed format
elif prog_input == "hbplus_parse_baseset":
    HbPlusProcesser().hbplusprocessedreader("hb")
    HbPlusProcesser().hbplusprocessedreader("vdw")

#combine the hbplus baseset hydrogen bonding and van der Waals results into a single file
elif prog_input == "hbplushbvdwcombine":
    HbPlusProcesser().hbplushbandvdwcombiner()

elif prog_input == "hbplus_baseset":
    HbPlusCliHelper().hbplus_hb_cli()
    HbPlusCliHelper().hbplus_vdw_cli()
    HbPlusProcesser().hbplusprocessedreader("hb")
    HbPlusProcesser().hbplusprocessedreader("vdw")
    HbPlusProcesser().hbplushbandvdwcombiner()

elif prog_input == "dssrcli_baseset":
    DssrCliHelper().dssrcli()
elif prog_input == "dssrparse_baseset":
    DssrParser().dssrprocessedreader()

elif prog_input == "dssr_cli_parse_baseset_combine":
    DssrCliHelper().dssrcli()
    DssrParser().dssrprocessedreader()

elif prog_input == "hbcategorizedssr_baseset":
    HbPlusToDssrComparer().hbplushbvdwtodssrcomparer()

#method to only count the bounds
elif prog_input == "countbonds":
    Bondcounter().file_opener()
    # Bondcounter().total_atom_counter()
elif prog_input == "statistical_potential":
    StatisticalPotential().statistical_potential()

#Running it together does not seem to be working
elif prog_input == "bound_count_stat_pot_baseset":
    Bondcounter().file_opener()
    StatisticalPotential().statistical_potential()

elif prog_input == "pdbsplit":
    FTDockChainStripping().pdbfile_splitter_rna_protein("protein_seperated_pdbfiles_baseset","rna_seperated_pdbfiles_baseset","base_complexes_pdb")
    FTDockChainStripping().pdbfile_splitter_rna_protein("protein_seperated_pdbfiles_testset","rna_seperated_pdbfiles_testset","test_complexes_pdb")
elif prog_input == "preprocessftdock":
    FTDockChainStripping().preprocess_ftdock("rna")
    FTDockChainStripping().preprocess_ftdock("protein")
# this is deprectaed
elif prog_input == "testprotein_pdb_combine":
    FTDockChainStripping().pdb_file_combine_rms_calc()

elif prog_input == "copyfilestosjsucluster":
    FTDockChainStripping().file_copier_to_sjsu_cluster_testset("base")
    # FTDockChainStripping().file_copier_to_sjsu_cluster_testset("test")
    # FTDockChainStripping().file_copier_to_sjsu_cluster("rna")
    # FTDockChainStripping().file_copier_to_sjsu_cluster("protein")
elif prog_input == "ftdockgen":
    FtdockMiddleware().ftdock_kicker()
elif prog_input == "ftdockbuild":
    # FtdockMiddleware().ftdock_directory_cleaner("home","cma","Users","curtisma")
    # FtdockMiddleware().ftdock_directory_cleaner("Users","curtisma","home","cma")
    FtdockMiddleware().ftdock_builder(50373)

elif prog_input == "hbpluscli_testset":
    HbplusCliTestset().hbpluscli("hbplus_processed_hb_files_testset","hb")
    HbplusCliTestset().hbpluscli("hbplus_processed_vdw_files_testset","vdw")

elif prog_input == "hbplusprocess_testset":
    HbPlusProcesserTestSet().hbplusprocessed_file_prepper_reader("hbplus_processed_hb_files_testset","hb")
    HbPlusProcesserTestSet().hbplusprocessed_file_prepper_reader("hbplus_processed_vdw_files_testset","vdw")

elif prog_input == "hbplusprocess_testset_combine":
    HbPlusProcesserTestSet().hbplushbandvdwcombiner()

elif prog_input == "dssrcli_testset":
    DssrCliHelperTestset().dssrcli_revamped()
elif prog_input == "dssrparse_testset":
    DssrParserTestSet().dssrprocessedreader()
elif prog_input == "bondcategorizer_testset":
    HbPlusToDssrComparerTestset().hbplushbvdwtodssrcomparer()

#run pymol wrapper for the native poses here....
#this is decprecated
elif prog_input == "native_checker":
    EnergyCalculator().native_checker_pose_hash_retriever("<protein>")

elif prog_input == "energy_calculate_testset":
    EnergyCalculator().energy_calculator("")
#    energyCalculator().energy_calculator("hb_only")
#    energyCalculator().energy_calculator("vdw_only")
elif prog_input == "energy_calculate_testset_ranking":
    EnergyCalculator().pose_rankings("")
#    energyCalculator().pose_rankings("hb_only")
#    energyCalculator().pose_rankings("vdw_only")

elif prog_input == "completerun":
    hbplusprocessedreader("hb")
    hbplusprocessedreader("vdw")
    hbplushbandvdwcombiner()
    DssrParser().dssrprocessedreader()
    hbplushbvdwtodssrcomparer()
    Bondcounter().bondcounter()
    Bondcounter().total_atom_counter()
    preprocess_ftdock("rna")
    preprocess_ftdock("protein")
