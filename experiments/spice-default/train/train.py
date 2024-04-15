import click
import espfit
from espaloma.graphs.utils.regenerate_impropers import regenerate_impropers
from espfit.utils.graphs import CustomGraphDataset  
from espfit.app.train import EspalomaModel


def run(kwargs):
    """Train Espaloma model on a dataset."""
    # Set parameters
    train_ratio = kwargs['train_ratio']
    path_to_data = kwargs['path_to_data']

    # Load dataset
    print(f'Get training dataset')
    dataset_small = CustomGraphDataset.load(f'{path_to_data}/small').split([train_ratio, 1-train_ratio])[0]
    dataset_peptide = CustomGraphDataset.load(f'{path_to_data}/peptide').split([train_ratio, 1-train_ratio])[0]
    print(f'Small molecule: {len(dataset_small)}')
    print(f'Peptide: {len(dataset_peptide)}')

    # Split dataset
    ds = dataset_small
    ds += dataset_peptide
    print(f'Total: {len(ds)}')

    # Prepare dataset
    print(f'Preprocess dataset...')
    ds.apply(regenerate_impropers, in_place=True)
    ds.compute_relative_energy()
    ds.reshape_conformation_size(n_confs=50, include_min_energy_conf=True)
    print(f'Final dataset size: {len(ds)}')

    # Create esplama model
    print(f'Create model')
    filename = 'config.toml'
    model = EspalomaModel.from_toml(filename)
    model.dataset_train = ds
    model.output_directory_path = 'checkpoints'

    # Change default training settings
    print(f'Train espaloma...')
    model.train()


@click.command()
@click.option('--train_ratio', default='0.8', help='Train ratio', type=float)
@click.option('--path_to_data', required=True, help='Path to QC datasets')
def cli(**kwargs):
    print(espfit.__version__)
    run(kwargs)


if __name__ == '__main__':
    cli()