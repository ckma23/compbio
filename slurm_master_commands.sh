#!/bin/bash
#SBATCH --job-name=stat_pot_curtis_ma_test
#SBATCH --cpus-per-task=8
#SBATCH --output=stat_pot.txt

python master.py statistical_poential
