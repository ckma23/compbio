import os
import time

class FtdockMiddleware(object):
    def ftdock_kicker (self):
        # os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/rna_seperated_pdbfiles_preprocessperl_testset'))
        # rna_seperated_pdbfiles = os.listdir()
        # os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_preprocessperl_testset'))
        # protein_seperated_pdbfiles = os.listdir()
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockresults'))
        # os.system("rm *.dat")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/test_complexes_pdb'))
        testset_pdb_files = os.listdir('.')
        # check that if this file EXISTS then kick off next ftdock, if not hold!
        for pdbfile in testset_pdb_files:
            print pdbfile
            pdbfile = pdbfile.strip("pdb")
            pdbfile = pdbfile[0:4]
            print pdbfile
            os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/ftdock-2-dev2/progs-2.0.3'))
            os.system("./ftdock -noelec -mobile ~/bioresearch/compbio/files_wip/rna_seperated_pdbfiles_preprocessperl_testset/%s_*.parsed -static ~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_preprocessperl_testset/%s_*.parsed -out ~/bioresearch/compbio/files_wip/ftdockresults/%sftdock.out > ~/bioresearch/compbio/library/ftdock_middleware/output &" %(pdbfile,pdbfile,pdbfile))
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
    #it turns out that the .out files have a callout where to find the static and mobile strings it came from.. for example...
    # Curtiss-MacBook-Pro:ftdockresults curtisma$ cat 1ec6ftdock.out | head -8
    # FTDOCK data file
    #
    # Global Scan
    #
    # Command line controllable values
    # Static molecule                    :: /home/cma/bioresearch/compbio/files_wip/rna_seperated_pdbfiles_preprocessperl_testset/1ec6_C.parsed
    # Mobile molecule                    :: /home/cma/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_preprocessperl_testset/1ec6_A.parsed
    #
    # Curtiss-MacBook-Pro:ftdockresults curtisma$
    def ftdock_directory_cleaner(self,root_path_cluster,user_path_cluster,root_path_local,user_path_local):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockresults'))
        ftdock_file_cleans = os.listdir('.')
        # ftdock_file_cleans = ["1e7kftdock.out"]
        for ftdock_file in ftdock_file_cleans:
            print ftdock_file
            print os.getcwd()
            print root_path_cluster
            print user_path_cluster
            print root_path_local
            print user_path_local
            os.system("sed -i.bak 's/%s/%s/g'  %s" %(root_path_cluster,root_path_local,ftdock_file))
            os.system("sed -i.bak 's/%s/%s/g'  %s" %(user_path_cluster,user_path_local,ftdock_file))
            os.system("rm *.bak")

    #the build_number deteremines the amount of poses that will be built into pdbs
    def ftdock_builder (self,build_number):
        os.system("mkdir ~/bioresearch/compbio/files_wip/ftdockbuiltposes")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockresults'))
        testset_ftdocked_pdbfiles = os.listdir('.')

        for testset_ftdocked_pdbfile in testset_ftdocked_pdbfiles:
            #make the directory for this to be built into
            os.system("mkdir ~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s" %testset_ftdocked_pdbfile[0:4])
            os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/ftdock-2-dev2/progs-2.0.3'))
            os.system("./build -in ~/bioresearch/compbio/files_wip/ftdockresults/%s -b1 1 -b2 %s > ~/bioresearch/compbio/logs/ftdock_middleware.txt" %(testset_ftdocked_pdbfile,build_number))
            # os.system("mv Complex_* ~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s" %testset_ftdocked_pdbfile[0:4])
            # os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes'))
            #
            # os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s' %testset_ftdocked_pdbfile[0:4]))
            # # os.system("rm Complex*.pdb")
            # os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/ftdock-2-dev2/progs-2.0.3'))
            for numb in range(1,55):
                os.system("mv Complex_%i* ~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s" %(numb,testset_ftdocked_pdbfile[0:4]))

            os.system("mv Complex*.pdb ~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s" %testset_ftdocked_pdbfile[0:4])

        time.sleep(5)
