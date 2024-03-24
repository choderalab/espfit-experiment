#!/usr/bin/env python
import click
import h5py


SELECTED_SUBSETS = ['SPICE PubChem', 'SPICE Dipeptides', 'SPICE DES Monomers']


def load_from_hdf5(kwargs):
    filename = kwargs["hdf5"]
    hdf = h5py.File(filename)

    count = 0
    with open('keys.txt', 'w') as wf:
        wf.write(f'HDF_Index\tArray_Index\tKey\tSubset\n')
        for i, key in enumerate(list(hdf.keys())):
            subset = hdf[key]['subset'][0].decode('UTF-8')
            if subset.startswith(SELECTED_SUBSETS[0]) or subset.startswith(SELECTED_SUBSETS[1]) or subset.startswith(SELECTED_SUBSETS[2]):
                wf.write(f'{i}\t{count}\t{key}\t{subset}\n')
                count += 1
                

@click.command()
@click.option("--hdf5", default="SPICE-2.0.0.hdf5", help='hdf5 filename')
def cli(**kwargs):
    load_from_hdf5(kwargs)


if __name__ == "__main__":
    cli()