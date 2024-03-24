#!/bin/bash
#BSUB -P "rna"
#BSUB -J "c"
#BSUB -n 1
#BSUB -R rusage[mem=8]
#BSUB -R span[hosts=1]
#BSUB -q gpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -gpu num=1:j_exclusive=yes:mode=shared
###BSUB -W 0:10
#BSUB -W 8:59
#BSUB -o out_%J_%I.stdout
#BSUB -eo out_%J_%I.stderr
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

# Fixed
script_path='/home/takabak/data/espfit-experiment/scripts/rna-nucleoside'
configfile='config.toml'

# Settings
small_molecule_forcefield='/home/takabak/.espaloma/espaloma-0.3.2.pt'
override_with_espaloma='True'  # False: amber, True: espaloma
suffix='c_tip3p'

# Run
conda activate test_env
python ${script_path}/sim.py --configfile ${configfile} --small_molecule_forcefield ${small_molecule_forcefield} --override_with_espaloma ${override_with_espaloma}
python ${script_path}/analysis.py --suffix ${suffix}

