#!/bin/bash
#BSUB -P "rna"
#BSUB -J "rna"
#BSUB -n 1
#BSUB -R rusage[mem=256]
#BSUB -R span[hosts=1]
#BSUB -q cpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -W 23:59
#BSUB -o rna_%J_%I.stdout
#BSUB -eo rna_%J_%I.stderr
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

# test
#python prep.py --category small --dataset gen2 --dataset gen2-torsion

# production
echo "# RNA"
python prep.py --category rna-nucleoside --dataset rna-nucleoside
python prep.py --category rna-diverse --dataset rna-diverse
python prep.py --category rna-trinucleotide --dataset rna-trinucleotide

