#!/bin/bash
#BSUB -P "esp"
#BSUB -J "metric"
#BSUB -n 1
#BSUB -R rusage[mem=512]
#BSUB -R span[hosts=1]
#BSUB -q cpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -W 12:00
#BSUB -o out_%J_%I.stdout
#BSUB -eo out_%J_%I.stderr
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


# paramters
train_ratio=0.8
ckptfile="../eval/ckpt.pt"
configfile="../config.toml"
path_to_data="/home/takabak/data/espfit-experiment/experiments/spice-openff-default/data"

# run jobs
conda activate espfit

echo "# Peptide"
python metric.py --train_ratio ${train_ratio} --ckptfile ${ckptfile} --configfile ${configfile} --path_to_data ${path_to_data} --category peptide

echo "# Small"
python metric.py --train_ratio ${train_ratio} --ckptfile ${ckptfile} --configfile ${configfile} --path_to_data ${path_to_data} --category small
