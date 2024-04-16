## Description
Each directory is a workspace for a different espaloma model (see main README.md for more details). 


## Manifest
The basic constuction of the directories are as follows:
- `spice-default/`
    - `data/` - download QC data and preprocess DGL graphs
        - `hdf5/` - download QC data as HDF5 file
    - `train/` - train and evaluate espaloma
        - `eval/` - find the best model
        - `metric/` - compute RMSE metric of QC energetic properties
    - `benchmark/` - perform benchmark experiments
        - `nucleoside/` - RNA nucleoside simulation
        - `pl-benchmark` - Protein-ligand alechmical binding free energy calculation
    - `baseline/` - (optional) compute RMSE metric for other force fields (e.g. openff-2.1.0 and amber14SB)

Note that depending on the espaloma model, the QC data is downloaded slightly differently. 
For example, for `openff-default`, the QC data in `data/` is gathered by soft linking to existing data sources within the hpc storage.


## Basic usage (e.g. `spice-default/`)
- Download QC data and convert them to DGL graphs (`data/hdf5/`)    
    >\# download SPICE-1.1.4 from Zenodo  
    > wget.sh  

    >\# run hdf2.graph.py to convert hdf5 to dgl graphs  
    >\# meta data of jobs that failed to compute the partial charges are exported to `partial_charge_failures.txt`  
    > bsub < lsf-submit.sh  

    (Optional)  
    `check_status.py`: check if there are any failed jobs by checking if output files exists  
    `getkeys.py`: get the meta data from the HDF5 file

    | HDF_Index | Array_Index | Key                                    | Subset                                      |
    |-----------|-------------|----------------------------------------|---------------------------------------------|
    | 0         | 0           | 103147721                              | SPICE PubChem Set 6 Single Points Dataset v1.2 |
    | 1         | 1           | 103147756                              | SPICE PubChem Set 2 Single Points Dataset v1.2 |


- Preprocess DGL graphs
    >\# run prep.py for small molecule species
    > bsub < lsf-submit-small.sh

- Train espaloma (`train/`)
    >\# run training using `train.py` and `config.toml`  
    > bsub < lsf-submit.sh

- Find best model (`train/eval/`)
    >\# use validation dataset to find the best model  
    >\# `rmse.csv` is created and energy and force losses at a given epoch will be appended to this file  
    >\# note that the epochs in .csv are not in acsending order but this should be changed in the future  
    >bsub < lsf-submit.sh  

    >\# plot validation loss and find the best model using early stopping  
    >\# saves the best model as `ckpt.pt`  
    >python plot.py > plot.log

- Compute RMSE metrics of QC energies and forces (`train/metric/`)
    >\# summary csv files (e.g. `summary_small.csv`) will be created along with other related files  
    >bsub < lsf-submit.sh  

- Save espaloma model for further validation and simulations  
    Move to `train/` but this could be any directory.  

    ```python
    from espfit.app.train import EspalomaModel
    model = EspalomaModel.from_toml('config.toml')
    # save model as `espaloma.pt` to current directory (default)
    model.save_model(checkpoint_file='./eval/ckpt.pt')
    ```