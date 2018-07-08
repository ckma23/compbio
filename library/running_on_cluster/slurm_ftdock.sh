#!/bin/bash
#SBATCH --job-name=curtis_ma_test
#SBATCH --cpus-per-task=4
#SBATCH --output=ftdock_log.txt

python ftdock_middleware_cluster.py
