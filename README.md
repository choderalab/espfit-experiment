# Workspace for espfit
This is a workspace to test and experiment [espfit](https://github.com/choderalab/espfit) which is a refactored version of [refit-espaloma](https://github.com/choderalab/refit-espaloma) to train and validate [espaloma](https://github.com/choderalab/espaloma). 
Espaloma is an end-to-end differentibale graph neural network framework to parameterize classical molecular mechanics force fields.


## Description
The motivation and final goal of this repository to refit espaloma using QC data and NMR experimental observables of RNA nucleosides to simulate RNA systems.


## Manifest
- `experiments/` - workspace for each espaloma model trained on different dataset and/or protocol
    - `espaloma-0.3.2/` - experiments using `espaloma-0.3.2` which is equivalent to the model created in Ref1.
    - `openff-default/` - experiments using espaloma model trained with the same dataset as `espaloma-0.3.2` but with [espfit](https://github.com/choderalab/espfit)
    - `spice-default/` - experiments using espaloma model trained with the [SPICE-1.1.4](https://zenodo.org/records/8222043) dataset
    - `spice-openff-default/` - experiments using espaloma model trained with the SPICE dataset included in the `espaloma-0.3.2` training dataset
    - `spice2-default/` - experiments using espaloma model trained with the [SPICE-2.0.0](https://zenodo.org/records/10835749) dataset
- `scripts/` - common scripts to run benchmark experiments
    - `pl-benchmark/` - alchemical protein-ligand binding free energy calculations
    - `rna-nucleoside/` - RNA nucleoside simulations


Note that `espaloma-0.3.2`, `openff-default`, and `spice-openff-default` uses QC data generated using ***B3LYP-D3BJ/DZVP*** level of theory whereas `spice-default` and `spice2-default` uses ***ωB97M-D3BJ/def2-TZVPPD***. See Ref2 for the details about the impact of QM level of theory.

## Reference
[1] [Takaba, K et al., Machine-learned molecular mechanics force field for the simulation of protein-ligand systems and beyond, 2023, arXiv:2307.07085 ](https://arxiv.org/abs/2307.07085)  
[2] [Behara, P. K. et al., Benchmarking QM theory for drug-like molecules to train force fields, 2022, OpenEye CUP XII, Santa Fe, NM. Zenodo](https://doi.org/10.5281/zenodo.7548777)