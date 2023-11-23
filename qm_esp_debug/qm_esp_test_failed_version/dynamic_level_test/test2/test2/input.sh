#!/bin/bash

#SBATCH -A DCCADD
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL


module load Anaconda3
chargecraft activate p4env

psi4 input.dat out.dat -n 4
