#!/bin/bash
#BSUB -P "eval"
#BSUB -J "eval-[1-100]"
#BSUB -n 1
#BSUB -R rusage[mem=256]
#BSUB -R span[hosts=1]
#BSUB -q cpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -W 2:59
#BSUB -o stdout/out_%J_%I.stdout
#BSUB -eo stderr/out_%J_%I.stderr
#BSUB -L /bin/bash

source ~/.bashrc
OPENMM_CPU_THREADS=1
#export OE_LICENSE=~/.openeye/oe_license.txt   # Open eye license activation/env

# change dir
echo "changing directory to ${LS_SUBCWD}"
cd $LS_SUBCWD


# Report node in use
echo "======================"
hostname
env | sort | grep 'CUDA'
nvidia-smi
echo "======================"

# Make output directories
mkdir -p stdout stderr

# conda
conda activate espfit
#conda activate test_env
python eval.py --path_to_data /home/takabak/data/espfit-experiment/experiments/spice-openff-default/data --index $(( $LSB_JOBINDEX ))

