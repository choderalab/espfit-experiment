#!/bin/bash
#BSUB -P "download"
#BSUB -J "download"
#BSUB -n 1
#BSUB -R rusage[mem=8]
#BSUB -R span[hosts=1]
#BSUB -q cpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -W 12:00
#BSUB -o out_%J_%I.stdout
#BSUB -eo out_%J_%I.stderr
#BSUB -L /bin/bash

source ~/.bashrc
OPENMM_CPU_THREADS=1
export OE_LICENSE=~/.openeye/oe_license.txt   # Open eye license activation/env


# change dir
echo "changing directory to ${LS_SUBCWD}"
cd $LS_SUBCWD


# Report node in use
echo "======================"
hostname
#env | sort | grep 'CUDA'
#nvidia-smi
echo "======================"


# run job
wget https://zenodo.org/records/10835749/files/SPICE-2.0.0.hdf5