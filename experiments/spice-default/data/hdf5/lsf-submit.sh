#!/bin/bash
#BSUB -P "esp"
####BSUB -J "esp-[1-15694]"
####BSUB -J "esp-[1-5000]"
####BSUB -J "esp-[5001-10000]"
#BSUB -J "esp-[10001-15694]"
#BSUB -n 1
#BSUB -R rusage[mem=8]
#BSUB -R span[hosts=1]
#BSUB -q cpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -W 1:00
#BSUB -o ./stdout/out_%J_%I.stdout
#BSUB -eo ./stderr/out_%J_%I.stderr
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
conda activate espfit
python hdf2graph.py --filename SPICE-1.1.4.hdf5 --array_index $(( $LSB_JOBINDEX - 1 ))