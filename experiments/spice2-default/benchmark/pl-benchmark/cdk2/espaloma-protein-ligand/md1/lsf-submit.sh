#!/bin/bash
#BSUB -P "cdk2-benchmark"
#BSUB -J "cdk2-[1-9]"
#BSUB -n 1
###BSUB -R rusage[mem=32]
#BSUB -R rusage[mem=8]
#BSUB -R span[hosts=1]
#BSUB -q gpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -gpu num=1:j_exclusive=yes:mode=shared
###BSUB -W  3:00
#BSUB -W  35:00
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
python ${script_path}/run_benchmarks.py --target cdk2 --edge $(( $LSB_JOBINDEX - 1 ))
