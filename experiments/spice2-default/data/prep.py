import click
import espfit
from espfit.utils.graphs import CustomGraphDataset  


def run(kwargs):
    dataset = kwargs['dataset']
    category = kwargs['category']
    
    # Load the dataset
    ds = CustomGraphDataset.load(f'hdf5/{dataset[0]}')
    if len(dataset) > 1:
        for _dataset in dataset[1:]:
            ds += CustomGraphDataset.load(f'hdf5/{_dataset}')

    print(f'(Default) reference force field: {ds.reference_forcefield}')
    ds.reference_forcefield = 'openff-2.1.0'
    print(f'(Updated) reference force field: {ds.reference_forcefield}')

    # Drop and merge duplicate molecules. Save merged dataset as a new dataset
    # If `output_directory_path` is None, then the current working directory is used
    ds.drop_duplicates(isomeric=False, keep=True, save_merged_dataset=True, dataset_name=category, output_directory_path='misc')
    # Subtract nonbonded energies and forces from QC reference (e.g. subtract all valence and ele interactions)
    # u_ref (QM reference) will be copied to u_qm/u_qm_prime and update u_ref/u_ref_prime in-place
    ds.subtract_nonbonded_interactions(subtract_vdw=False, subtract_ele=True)
    # Filter high energy conformers (u_qm: QM reference before nonbonded interations are subtracted)
    ds.filter_high_energy_conformers(relative_energy_threshold=0.1, node_feature='u_qm')
    # Filter high energy conformers (u_ref: QM reference after nonbonded interactions are subtracted)
    ds.filter_high_energy_conformers(relative_energy_threshold=0.1, node_feature='u_ref')
    # Filter conformers below certain number
    ds.filter_minimum_conformers(n_conformer_threshold=5)
    # Compute energies and forces using other force fields
    ds.compute_baseline_energy_force(forcefield_list=['openff-2.1.0', 'amber14-all.xml'])
    # Save the dataset
    ds.save(category)


@click.command()
@click.option("--category", required=True, type=click.Choice(['debug', 'small', 'peptide', 'rna']), help="category of the dataset")
@click.option("--dataset", required=True, multiple=True, help="name of the dataset")
def cli(**kwargs):
    run(kwargs)


if __name__ == "__main__":
    print(espfit.__version__)
    cli()