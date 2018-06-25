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
        os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/ftdock-2-dev2/progs-2.0.3'))
        for pdbfile in testset_pdb_files:
            pdbfile = pdbfile.strip("pdb")
            pdbfile = pdbfile.strip(".ent")
            # add a line to clean up the logs each time
            print pdbfile
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
