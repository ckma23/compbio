#!/bin/bash
#
#SBATCH --job-name=curtis_ma_test
#SBATCH --output=ftdock_log.txt
python ftdock_middleware_cluster.py
