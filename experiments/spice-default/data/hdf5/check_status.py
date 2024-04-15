#!/usr/bin/env python
import os
import glob


def check_file_exists():
    try:
        # key, mapped_smiles, hdf_index, array_index
        failures = open('partial_charge_failures.txt', 'r')
        charge_failure_keys = [ l.split()[0] for l in failures.readlines() ]
    except:
        charge_failure_keys = []

    import pandas as pd
    df = pd.read_csv('keys.txt', sep='\t')
    for _, col in df.iterrows():
        key = col['Key']
        hdf_index = col['HDF_Index']
        array_index = col['Array_Index']
        subset = col['Subset']

        # Check if file exists
        if key not in charge_failure_keys:
            path = glob.glob('spice*/' + str(hdf_index))
            if len(path) != 1:
                raise ValueError(f"Found {len(path)} files")
            else:
                path = path[0]
                molfile = os.path.join(path, "mol.json")
                heterograph = os.path.join(path, "heterograph.bin")
                homograph = os.path.join(path, "homograph.bin")
                if os.path.exists(molfile) and os.path.exists(heterograph) and os.path.exists(homograph):
                    pass
                else:
                    print(f'Subset: {subset}\tHDF Index: {hdf_index}\tArray Index: {array_index}')


if __name__ == "__main__":
    check_file_exists()
