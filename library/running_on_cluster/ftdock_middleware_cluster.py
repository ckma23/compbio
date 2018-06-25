import os
import time

class FtdockMiddleware(object):
    def ftdock_kicker (self):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip'))
        os.system("mkdir ftdockresults")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockresults'))
        os.system("rm *.dat")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/test_complexes_pdb'))
        testset_pdb_files = os.listdir('.')
	    testset_pdb_files =["pdb1e7k.ent","pdb1f7u.ent"]
        os.chdir(os.path.expanduser('~/curtisma/bioresearch/compbio/bin/ftdock-2-dev2/progs-2.0.3'))
        for pdbfile in testset_pdb_files:
            pdbfile = pdbfile.strip("pdb")
            pdbfile = pdbfile.strip(".ent")
            print pdbfile
            os.system("./ftdock -noelec -static ~/bioresearch/compbio/rna_seperated_pdbfiles_preprocessperl_testset/%s_*.parsed -mobile ~/bioresearch/compbio/bioresearch/protein_seperated_pdbfiles_preprocessperl_testset/%s_*.parsed -out ~/bioresearch/compbio/ftdockresults/%sftdock.out > ~/bioresearch/compbio/logs/ftdock_output &" %(pdbfile,pdbfile,pdbfile))
            # sleep FTdock for 5 mins
            time.sleep(60)

FtdockMiddleware().ftdock_kicker()
