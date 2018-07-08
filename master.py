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
from helper.biopythonpdbretriever import FileRetriever      #from the helper folder biopythonpdbretriever.py import FileRetriever class
from helper.PDBRestService import PDBRestServicehelper
from library.hbplusclilib.hbpluscli import hbplusclihelper
from library.hbplusclilib.hbpluscli_testset import HbplusCliTestset as HbplusCliTestset  #This is just placeholder for now.
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

from library.bondfinalcount.bondcounter import Bondcounter
from library.ftdockprepper.chain_stripping import FTDockChainStripping
from library.ftdock_middleware.ftdock_middleware import FtdockMiddleware
from optparse import OptionParser as optionparser # will need a option parsing library

def directory_builder():
    os.chdir(os.path.expanduser('~/bioresearch/compbio'))
    #make the directory logs
    os.system("mkdir logs")
    #make the work in progress folder directories used
    os.system("mkdir files_wip")

def pdbfileretriever(configuration_file,folder_path_name):
    with open(configuration_file,'rb') as csvfile:
        csvstore = csv.reader(csvfile,delimiter = ',')
        for structure in csvstore:
            print (structure)
            # note: had to strip the string as it was returning ['<string>'] for []
            FileRetriever().fileretrieving(''.join(structure),folder_path_name)

### DECISION TREE FOR CATEGORIES ###
#CAT 1 : Helix => A Form ("A") => major grove => major groove like => canonical basepair (U to A; C to G) => check atom type
#CAT 2 : Helix => A Form ("A") => not major grove => not a canonical base pair (U to G) => check atom type
#CAT 3 : Helix => A Form ("A") => backbone ("C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3")
#CAT 4 : Helix => Not A Form ("B","Z",".","x") => major groove like => canonical basepair (U to A; C to G)
#CAT 5 : Helix => Not A Form ("B","Z",".","x") => major groove like => not a canonical base pair (U to G)
#CAT 6 : Helix => Not A Form ("B","Z",".","x") => not major groove like
#CAT 7 : Helix => Not A Form ("B","Z",".","x") => backbone ("C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3")
#CAT 8 : Not Helix (Hairpin, Bulge, Loops) => base
#CAT 9 : Not Helix (Hairpin, Bulge, Loops) => backbone ("C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3")

# assumptions is that A Form Helix must likely only be a canonical base pair or wobble pair.
# thus if it's not a canonical basepair then check
# CAT 1 or CAT 2
# check atom type in A-form and which atom types are major groove in A-form and which ones are not; thus CAT 1 or CAT 2
# CAT 4 or CAT 5 or CAT 6
# CAT 4: canonical base pair is CAT 4 by cW-w
# CAT 5: non cW-w and check the atom if it's in the major groove or minor groove
# CAT 6: non cW-w and atom is not in major groove


### DECISION TREE TWO FOR CATEGORIES ###
#CAT 1 : Helix => A Form ("A") => canonical basepair ("cW-W" where U to A; C to G; wobble) => check atom type for major groove or minor groove
#CAT 2 : Helix => A Form ("A") => not a canonical base pair or  canonical basepair atom type is on minor groove
#CAT 3 : Helix => A Form ("A") => backbone ("C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3")
#CAT 4 : Helix => Not A Form ("B","Z",".","x") => canonical basepair ("cW-W" where U to A; C to G; wobble) => check atom type for major groove or minor groove
#CAT 5 : Helix => Not A Form ("B","Z",".","x") => not a canonical base pair => check atom type for major groove or minor groove
#CAT 6 : Helix => Not A Form ("B","Z",".","x") => not major groove like
#CAT 7 : Helix => Not A Form ("B","Z",".","x") => backbone ("C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3")
#CAT 8 : Not Helix (Hairpin, Bulge, Loops) => base
#CAT 9 : Not Helix (Hairpin, Bulge, Loops) => backbone ("C1\'","C2\'","C3\'","C4\'","C5\'","O3\'","O4\'","O5\'","H1\'","H2\'","H3\'","H4\'","H5\'","H5\'\'","P","OP1","OP2","OP3")

## major groove atoms in A are: H61,H62,N6,C6,C5,N7,C8,H8
## minor groove atoms in A are: N2,C2,C5,C4
## major groove atoms in U are: O4,C4,C5,H5,C6,H6
## minor groove atoms in U are: H3,N3,C2,O2,N1
## major groove atoms in C are: H6,C6,C5,H5,C4,N4,H42,H41,N3
## minor groove atoms in C are: N1,C2,O2
## major groove atoms in G are: H1,N2,C6,O6,C5,N7,C8,H8
## minor groove atoms in G are: N2,C2,C5,C4

def helpoutput():
  print "\nWelcome to the Computatinal Biology RNA-Protein Interaction help section\n"
  print "\nThe following commands are available:\n"
  print "\ndssrcli \nfileretriever \nhbplushbcli \nhbplusvdwcli \ndssrprocessedreader \nhbplusprocessedreader\n"

# testprotein = structureRetriever('1mnb','pdbfiles/pdb1mnb.ent')
# neighborSearcher(testprotein,3.0)
# aminoacidmatcher()

#from the system take the second argument
proinput=str(sys.argv[1])
#default help menu
if proinput == "help":
    helpoutput()
elif proinput == "initialize":
    directory_builder()
#function to split the pdb files into rna only and protein only
elif proinput == "pdbfileretriever":
    pdbfileretriever("structures_of_interest.csv","files_wip/base_complexes_pdb")
    pdbfileretriever("test_complexes_pdb.csv","files_wip/test_complexes_pdb")
#need to build hbplus initially for the baseset
elif proinput == "hbplushbcli":
    hbplusclihelper().hbplushbcli()
#need to build hbplus for the testset and complexes
elif proinput == "hbplusvdwcli":
    hbplusclihelper().hbplusvdwcli()
elif proinput == "hbplus_parse_baseset":
    HbPlusProcesser().hbplusprocessedreader("hb")
    HbPlusProcesser().hbplusprocessedreader("vdw")
# elif proinput == "hbplushbparse":
#     hbplusprocessedreader("hb")
# elif proinput == "hbplusvdwparse":
#     hbplusprocessedreader("vdw")
elif proinput == "hbplushbvdwcombine":
    HbPlusProcesser().hbplushbandvdwcombiner()
elif proinput == "dssrcli":
    dssrclihelper().dssrcli()
elif proinput == "dssrparse":
    DssrParser().dssrprocessedreader()
elif proinput == "hbcategorizedssr":
    HbPlusToDssrComparer().hbplushbvdwtodssrcomparer()
    # hbplushbvdwtodssrcomparer()
elif proinput == "countbonds":
    Bondcounter().bondcounter()
    Bondcounter().total_atom_counter()
elif proinput == "statistical_potential":
    StatisticalPotential().statistical_potential()


elif proinput == "pdbsplit":
    FTDockChainStripping().pdbfile_splitter_rna_protein("protein_seperated_pdbfiles_baseset","rna_seperated_pdbfiles_baseset","base_complexes_pdb")
    FTDockChainStripping().pdbfile_splitter_rna_protein("protein_seperated_pdbfiles_testset","rna_seperated_pdbfiles_testset","test_complexes_pdb")
elif proinput == "preprocessftdock":
    FTDockChainStripping().preprocess_ftdock("rna")
    FTDockChainStripping().preprocess_ftdock("protein")
elif proinput == "copyfilestosjsucluster":
    FTDockChainStripping().file_copier_to_sjsu_cluster_testset("base")
    FTDockChainStripping().file_copier_to_sjsu_cluster_testset("test")
    FTDockChainStripping().file_copier_to_sjsu_cluster("rna")
    FTDockChainStripping().file_copier_to_sjsu_cluster("protein")
elif proinput == "ftdockgen":
    FtdockMiddleware().ftdock_kicker()
elif proinput == "ftdockbuild":
    FtdockMiddleware().ftdock_directory_cleaner("home","cma","Users","curtisma")
    FtdockMiddleware().ftdock_builder(10)



elif proinput == "hbpluscli_testset":
    HbplusCliTestset().hbpluscli("hbplus_processed_hb_files_testset","hb")
    HbplusCliTestset().hbpluscli("hbplus_processed_vdw_files_testset","vdw")
elif proinput == "hbplusprocess_testset":
    HbPlusProcesserTestSet().hbplusprocessed_file_prepper_reader("hbplus_processed_hb_files_testset","hb")
    HbPlusProcesserTestSet().hbplusprocessed_file_prepper_reader("hbplus_processed_vdw_files_testset","vdw")
    HbPlusProcesserTestSet().hbplushbandvdwcombiner()
elif proinput == "dssrcli_testset":
    DssrCliHelperTestset().dssrcli()
elif proinput == "dssrparse_testset":
    DssrParserTestSet().dssrprocessedreader()
elif proinput == "bondcategorizer_testset":
    HbPlusToDssrComparerTestset().hbplushbvdwtodssrcomparer()
elif proinput == "energy_calculate_testset":
    energyCalculator().energy_calculator("")
    energyCalculator().energy_calculator("hb_only")
    energyCalculator().energy_calculator("vdw_only")
elif proinput == "energy_calculate_testset_ranking":
    energyCalculator().pose_rankings("")
    energyCalculator().pose_rankings("hb_only")
    energyCalculator().pose_rankings("vdw_only")
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
