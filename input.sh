#!/bin/bash

#SBATCH -A DCCADD
#SBATCH --mail-type=ALL
#SBATCH --export=ALL
#SBATCH -c 4


chargecraft activate /mnt/nfs/home/nca121/mambaforge/envs/openff

export PYTHONPATH=/mnt/storage/nobackup/nca121/test_jobs/QM_ESP_Psi4/QM_ESP_Psi4/chargecraft:$PYTHONPATH

python /mnt/storage/nobackup/nca121/test_jobs/QM_ESP_Psi4/QM_ESP_Psi4/run_example.py
