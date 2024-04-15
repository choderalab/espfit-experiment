import os
import click
import glob


def run(kwargs):
    dataset = kwargs['dataset']
    basepath = '/home/takabak/data/exploring-rna/refit-espaloma/openff-default/01-create-dataset'

    # Create dataset directory
    os.makedirs(dataset, exist_ok=True)

    # Create softlink
    files = glob.glob(f'{basepath}/{dataset}/data/*/mydata')
    for file in files:        
        src = file
        
        idx = file.split('/')[-2]
        dst = os.path.join(dataset, idx)
        
        #print(src, dst)
        os.symlink(src, dst)


@click.command()
@click.option('--dataset', required=True, help='name of dataset in the form of [QCWorkflow]/[DATASET]')
def cli(**kwargs):
    run(kwargs)

if __name__ == '__main__':
    cli()