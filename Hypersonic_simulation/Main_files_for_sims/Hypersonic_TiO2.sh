#!/bin/bash

#SBATCH --job-name=TiO2_300_1
#SBATCH --time=24:00:00
#SBATCH --mem=4096
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lkm41@cam.ac.uk

python3 Hypersonic_TiO2_Python.py