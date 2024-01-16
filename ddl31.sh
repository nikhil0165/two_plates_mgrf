#!/bin/bash

#SBATCH --job-name=ocr31
#
# Account:
#SBATCH --account=fc_ionactivity
#
# Partition:
#SBATCH --partition=savio2
#
# Request one node:
#SBATCH --nodes=1
#
# Specify number of tasks for use case (example):
#SBATCH --ntasks-per-node=20

# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=10:30:30
#SBATCH --output=ocr31%j.out

#command
echo "this job starts at"

date

source activate phd

python3 -u simulator.py physical_param.py

echo "this job ends at"

date
##

