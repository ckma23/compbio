#!/bin/sh
#
#SBATCH --job-name=curtis_ma_test
#SBATCH --output=~/bioresearch/compbio/logs/ftdock_log.txt
sbatch python ftdock_middleware.py
