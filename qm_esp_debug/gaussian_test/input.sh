#!/bin/bash

#SBATCH -A DCCADD
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL


module load gaussian/09
module load ~/mambaforge/condabin/conda

conda activate /mnt/nfs/home/nca121/mambaforge/envs/openff
#chargecraft activate openff

python ase_test.py -n 4

