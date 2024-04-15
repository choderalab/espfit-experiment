#!/bin/bash
#BSUB -P "peptide"
#BSUB -J "peptide"
#BSUB -n 1
#BSUB -R rusage[mem=64]
#BSUB -R span[hosts=1]
#BSUB -q cpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -W 23:59
#BSUB -o peptide_%J_%I.stdout
#BSUB -eo peptide_%J_%I.stderr
#BSUB -L /bin/bash

source ~/.bashrc
export OE_LICENSE=~/.openeye/oe_license.txt   # Open eye license activation/env

# change dir
echo "changing directory to ${LS_SUBCWD}"
cd $LS_SUBCWD

# Report node in use
echo "======================"
hostname
env | sort | grep 'CUDA'
nvidia-smi
echo "======================"

# run
conda activate test_env
echo "# Peptide"
python prep.py --category peptide --dataset pepconf-dlc --dataset protein-torsion --dataset spice-dipeptide
