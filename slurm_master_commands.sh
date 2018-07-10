#!/bin/bash
#SBATCH --job-name=curtis_ma_test
#SBATCH --cpus-per-task=8
#SBATCH --output=ftdock_build.txt

python master.py ftdockbuild
