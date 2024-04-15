#!/bin/bash

# Dataset
python symlink.py --dataset dataset/rna-diverse
python symlink.py --dataset dataset/rna-nucleoside
python symlink.py --dataset dataset/rna-trinucleotide
python symlink.py --dataset dataset/spice-des-monomers
python symlink.py --dataset dataset/spice-dipeptide
python symlink.py --dataset dataset/spice-pubchem

# OptimizationDataset
python symlink.py --dataset dataset/gen2
python symlink.py --dataset dataset/pepconf-dlc

# TorsionDriveDataset
python symlink.py --dataset dataset/gen2-torsion
python symlink.py --dataset dataset/protein-torsion
