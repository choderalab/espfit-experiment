#!/usr/bin/env python
import os
import click
import numpy as np
import h5py
import torch
import espaloma as esp
from espaloma.units import *
from openff.toolkit.topology import Molecule
from simtk import unit
from simtk.unit import Quantity


def get_graph(record, key, hdf_index, array_index):
    """
    Convert HDF5 entries into dgl graphs.
    """

    # Convert mapped smiles to openff molecules
    try:
        mapped_smiles = record["smiles"][0].decode('UTF-8')
    except:
        # This was necessarily for GEN2 Datasets from openff default spec
        # May not be needed for Spice default spec but just in case
        mapped_smiles = record["smiles"][0][0].decode('UTF-8')
    offmol = Molecule.from_mapped_smiles(mapped_smiles, allow_undefined_stereo=True)

    # Compute AM1-BCC ELF10 using openeye-toolkit
    try:
        offmol.assign_partial_charges(partial_charge_method="am1bccelf10")
    except:
        failures = open('partial_charge_failures.txt', 'a')
        failures.write("{}\t{}\t{}\t{}\n".format(key, mapped_smiles, hdf_index, array_index))
        failures.close()
        msg = 'could not assign partial charge'
        raise ValueError(msg)
    charges = offmol.partial_charges.to_openmm().value_in_unit(esp.units.CHARGE_UNIT)   # openmm.unit.unit.Unit
    g = esp.Graph(offmol)

    energy = record["dft_total_energy"]
    grad = record["dft_total_gradient"]
    conformations = record["conformations"]

    g.nodes["g"].data["u_ref"] = torch.tensor(
        [
            Quantity(
                _energy,
                esp.units.HARTREE_PER_PARTICLE,
            ).value_in_unit(esp.units.ENERGY_UNIT)
            for _energy in energy
        ],
        dtype=torch.float64,
    )[None, :]

    g.nodes["n1"].data["xyz"] = torch.tensor(
        np.stack(
            [
                Quantity(
                    xyz,
                    unit.bohr,
                ).value_in_unit(esp.units.DISTANCE_UNIT)
                for xyz in conformations
            ],
            axis=1,
        ),
        requires_grad=True,
        dtype=torch.float32,
    )

    g.nodes["n1"].data["u_ref_prime"] = torch.stack(
        [
            torch.tensor(
                Quantity(
                    _grad,
                    esp.units.HARTREE_PER_PARTICLE / unit.bohr,
                ).value_in_unit(esp.units.FORCE_UNIT),
                dtype=torch.float32,
            )
            for _grad in grad
        ],
        dim=1,
    )

    g.nodes['n1'].data['q_ref'] = torch.tensor(charges, dtype=torch.float32,).unsqueeze(-1)
    
    return g


def load_from_hdf5(kwargs):

    filename = kwargs["filename"]
    array_index = int(kwargs["array_index"])
    hdf = h5py.File(filename)

    import pandas as pd
    df = pd.read_csv('keys.txt', sep='\t')
    key = str(df.iloc[array_index]['Key'])
    hdf_index = int(df.iloc[array_index]['HDF_Index'])
    subset = str(df.iloc[array_index]['Subset'])
    assert array_index == df.iloc[array_index]['Array_Index']

    SELECTED_SUBSETS = ['SPICE PubChem Set', 'SPICE Dipeptides', 'SPICE DES Monomers']
    if subset.startswith(SELECTED_SUBSETS[0]):
        output_prefix = os.path.join("spice-pubchem", f"{hdf_index}")
    elif subset.startswith(SELECTED_SUBSETS[1]):
        output_prefix = os.path.join("spice-dipeptide", f"{hdf_index}")
    elif subset.startswith(SELECTED_SUBSETS[2]):
        output_prefix = os.path.join("spice-des-monomers", f"{hdf_index}")

    try:
        record = hdf[key]
    #except KeyError:
    #    print(f"Invalid key ({key}). Get key using HDF entry index ({hdf_index}).")
    #    key = list(hdf.keys())[int(hdf_index)]
    #    record = hdf[key]
    except Exception as e:
        print(f'Failed to load record for key {key}')
        with open('failed_keys.txt', 'a') as f:
            f.write(f'HDF_Index: {hdf_index}\tArray_Index: {array_index}\tKey: {key}\n')

    g = get_graph(record, key, hdf_index, array_index)
    g.save(output_prefix)


@click.command()
@click.option("--filename", required=True, help='hdf5 filename')
@click.option("--array_index", required=True, help='index of array to load from hdf5 file')
def cli(**kwargs):
    print(esp.__version__)
    load_from_hdf5(kwargs)


if __name__ == "__main__":
    cli()
