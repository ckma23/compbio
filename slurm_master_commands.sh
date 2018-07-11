#!/bin/bash
#SBATCH --job-name=dssrcurtis_ma_test
#SBATCH --cpus-per-task=8
#SBATCH --output=dssr.txt

python master.py dssrcli_testset
