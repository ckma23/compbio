import os
import time

class FtdockMiddleware(object):
    def ftdock_kicker (self):
        # os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/rna_seperated_pdbfiles_preprocessperl_testset'))
        # rna_seperated_pdbfiles = os.listdir()
        # os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_preprocessperl_testset'))
        # protein_seperated_pdbfiles = os.listdir()
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockresults'))
        os.system("rm *.dat")
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/test_complexes_pdb'))
        testset_pdb_files = os.listdir('.')
        # check that if this file EXISTS then kick off next ftdock, if not hold!
        for pdbfile in testset_pdb_files:
            pdbfile = pdbfile.strip("pdb")
            pdbfile = pdbfile.strip(".ent")
            print pdbfile
            os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/ftdock-2-dev2/progs-2.0.3'))
            os.system("./ftdock -noelec -static ~/bioresearch/compbio/files_wip/rna_seperated_pdbfiles_preprocessperl_testset/%s_*.parsed -mobile ~/bioresearch/compbio/files_wip/protein_seperated_pdbfiles_preprocessperl_testset/%s_*.parsed -out ~/bioresearch/compbio/files_wip/ftdockresults/%sftdock.out > ~/bioresearch/compbio/library/ftdock_middleware/output &" %(pdbfile,pdbfile,pdbfile))
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


    def ftdock_builder (self):
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockresults'))
        testset_ftdocked_pdbfiles = os.listdir('.')
        os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip'))
        os.system("mkdir ftdockbuiltposes")
        for testset_ftdocked_pdbfile in testset_ftdocked_pdbfiles:
            os.chdir(os.path.expanduser('~/bioresearch/ftdock-2-dev2/progs-2.0.3'))
            # os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/ftdock-2-dev2/progs-2.0.3'))
            os.system("./build -in ~/bioresearch/compbio/files_wip/ftdockresults/%s -b1 1 -b2 50373 > ~/bioresearch/compbio/library/ftdock_middleware/cd ..output &" %testset_ftdocked_pdbfile)
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes'))
            os.system("mkdir %s" %testset_ftdocked_pdbfile[0:4])
            os.chdir(os.path.expanduser('~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s' %testset_ftdocked_pdbfile[0:4]))
            os.system("rm Complex*.pdb")
            os.chdir(os.path.expanduser('~/bioresearch/ftdock-2-dev2/progs-2.0.3'))
            # os.chdir(os.path.expanduser('~/bioresearch/compbio/bin/ftdock-2-dev2/progs-2.0.3'))
            for numb in range(1,55):
                os.system("mv Complex_%i* ~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s" %(numb,testset_ftdocked_pdbfile[0:4]))
            os.system("mv Complex*.pdb ~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s" %testset_ftdocked_pdbfile[0:4])
            # built_pose_filenames = os.listdir('.')

            # for now, lets not rename the files
            # for file_to_be_renamed in built_pose_filenames:
            #     print file_to_be_renamed
            #     print testset_ftdocked_pdbfile[0:4]
            #     os.system("mkdir %s" %test_ftdocked_pdbfile)
            #     os.system("mv Complex*.pdb ~/bioresearch/compbio/files_wip/ftdockbuiltposes/%s" %s)
                # os.system("mv ~/bioresearch/compbio/files_wip/ftdockresults/ftdockbuiltposes")

        time.sleep(5)
