import os
import time

class FtdockMiddleware(object):
    def ftdock_kicker (self):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip'))
        os.system("mkdir ftdockresults")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockresults'))
        ftdocked_files = os.listdir('.')
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/test_complexes_pdb'))
        testset_pdb_files = set(os.listdir('.'))
        # essentially here is no need to re-run FTdock for files we've created already.
        for ftdocked_file_already in ftdocked_files:
            # copy the set because if not the stacktrace is RuntimeError: Set changed size during iteration
            for test_complexes_pdb in testset_pdb_files.copy():
                if ftdocked_file_already [0:3] == test_complexes_pdb[3:6]:
                    # remove works for set.
                    testset_pdb_files.remove(test_complexes_pdb)
        print ftdocked_files
        print testset_pdb_files
        # test_pdb_files = ["",""]
        for pdbfile in testset_pdb_files:
            pdbfile = pdbfile.strip("pdb")
            pdbfile = pdbfile[0:4]
            # add a line to clean up the logs each time
            print pdbfile
            os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/ftdock-2-dev2/progs-2.0.3'))
            os.system("./ftdock -noelec -static ~/bioresearch/compbio/files_wip/rna_seperated_pdbfiles_preprocessperl_testset/%s_*.parsed -mobile ~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_preprocessperl_testset/%s_*.parsed -out ~/bioresearch/compbio/files_wip/ftdockresults/%sftdock.out > ~/bioresearch/compbio/logs/ftdock_output &" %(pdbfile,pdbfile,pdbfile))
            # sleep FTdock for 5 mins
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockresults'))
            boolean = True
            while boolean:
                ftdockresults = os.listdir('.')
                ftdockresults = set(ftdockresults)
                print ftdockresults
                if "%sftdock.out" %pdbfile in ftdockresults:
                    boolean = False
                    break
                else:
                    print "Not in there yet! Sleeping!"
                    boolean = True
                    time.sleep(60)

FtdockMiddleware().ftdock_kicker()
