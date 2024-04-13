#!/bin/bash
#BSUB -P "p38-benchmark"
#BSUB -J "p38"
#BSUB -n 1
#BSUB -R rusage[mem=8]
#BSUB -R span[hosts=1]
#BSUB -q gpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -gpu num=1:j_exclusive=yes:mode=shared
#BSUB -W  3:00
#BSUB -o analysis_%J_%I.stdout
#BSUB -eo analysis_%J_%I.stderr
#BSUB -L /bin/bash

source ~/.bashrc
OPENMM_CPU_THREADS=1

echo "changing directory to ${LS_SUBCWD}"
cd $LS_SUBCWD
conda activate perses-espaloma-0.3.0

script_path="/home/takabak/data/espfit-experiment/scripts/pl-benchmark"
python ${script_path}/benchmark_analysis.py --target p38
