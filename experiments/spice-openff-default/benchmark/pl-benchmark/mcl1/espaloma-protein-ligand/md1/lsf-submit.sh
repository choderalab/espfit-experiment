#!/bin/bash
#BSUB -P "mcl1-benchmark"
###BSUB -J "mcl1-[1-24]"
#BSUB -J "mcl1-[1]"
#BSUB -n 1
###BSUB -R rusage[mem=32]
#BSUB -R rusage[mem=8]
#BSUB -R span[hosts=1]
#BSUB -q gpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -gpu num=1:j_exclusive=yes:mode=shared
###BSUB -W  2:00
#BSUB -W  24:00
#BSUB -o stdout/out_%J_%I.stdout
#BSUB -eo stderr/out_%J_%I.stderr
#BSUB -L /bin/bash

source ~/.bashrc
OPENMM_CPU_THREADS=1

echo "changing directory to ${LS_SUBCWD}"
cd $LS_SUBCWD
conda activate perses-espaloma-0.3.0-v3

# Report node in use
hostname

# Report CUDA info
env | sort | grep 'CUDA'

mkdir -p stdout stderr
script_path="/home/takabak/data/espfit-experiment/scripts/pl-benchmark"
# launching a benchmark pair (target, edge) per job (0-based thus substract 1)
python ${script_path}/run_benchmarks.py --target mcl1 --edge $(( $LSB_JOBINDEX - 1 ))
