#!/bin/bash
#BSUB -P "small"
#BSUB -J "small"
#BSUB -n 1
#BSUB -R rusage[mem=256]
#BSUB -R span[hosts=1]
#BSUB -q cpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -W 47:59
#BSUB -o small_%J_%I.stdout
#BSUB -eo small_%J_%I.stderr
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
echo "# Small"
python prep.py --category small --dataset gen2 --dataset gen2-torsion --dataset spice-des-monomers --dataset spice-pubchem
