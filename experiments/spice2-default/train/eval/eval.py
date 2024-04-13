import os
import glob
import click
import torch
import espfit
import numpy as np
from espaloma.graphs.utils.regenerate_impropers import regenerate_impropers
from espfit.utils.graphs import CustomGraphDataset  
from espfit.app.train import EspalomaModel
from espfit.utils.units import HARTREE_TO_KCALPERMOL, BOHR_TO_ANGSTROMS


def run(kwargs):
    train_ratio = kwargs['train_ratio']
    path_to_data = kwargs['path_to_data']
    path_to_config = kwargs['path_to_config']
    path_to_checkpoints = kwargs['path_to_checkpoints']
    index = kwargs['index']

    # Load dataset
    ds_vl, ds_te = load_dataset(train_ratio, path_to_data)

    # Create model
    filename = os.path.join(path_to_config, "config.toml")
    model = EspalomaModel.from_toml(filename)
    net = model.net

    # Compute RMSE
    if not os.path.exists('rmse.csv'):
        with open('rmse.csv', 'w') as wf:
            wf.write('Epoch\tValidation_Total\tValidation_Energy\tValidation_Force\tTest_Total\tTest_Energy\tTest_Force\n')
    
    epoch = index * model.checkpoint_frequency
    file = os.path.join(path_to_checkpoints, f"ckpt{epoch}.pt")
    state_dict = torch.load(file, map_location=torch.device('cpu'))
    net.load_state_dict(state_dict)
    eval_loss(epoch, net, ds_vl, ds_te)


def load_dataset(train_ratio, path_to_data):
    """Train Espaloma model on a dataset."""
    # Load dataset
    print(f'Split dataset')
    
    _, small_vl_te = CustomGraphDataset.load(f'{path_to_data}/small').split([train_ratio, 1-train_ratio])
    small_vl_te.apply(regenerate_impropers, in_place=True)
    small_vl, small_te = small_vl_te.split([0.5, 0.5])
    print(f'Small molecule: {len(small_vl)}, {len(small_te)}')
    
    _, peptide_vl_te = CustomGraphDataset.load(f'{path_to_data}/peptide').split([train_ratio, 1-train_ratio])
    peptide_vl_te.apply(regenerate_impropers, in_place=True)
    peptide_vl, peptide_te = peptide_vl_te.split([0.5, 0.5])
    print(f'Peptide: {len(peptide_vl)}, {len(peptide_te)}')
    
    # Merge dataset
    ds_vl = small_vl + peptide_vl
    ds_te = small_te + peptide_te

    # Compute relative
    ds_vl.compute_relative_energy()
    ds_te.compute_relative_energy()

    return ds_vl, ds_te


def eval_loss(epoch, net, ds_vl, ds_te):
    net.eval()

    e_vl, f_vl = [], []
    e_te, f_te = [], []

    for g in ds_vl:
        g.nodes["n1"].data["xyz"].requires_grad = True
        e_loss, _ = net(g.heterograph)
        e_vl.append(HARTREE_TO_KCALPERMOL * e_loss.pow(0.5).item())
        du_dx_hat = torch.autograd.grad(
            g.nodes['g'].data['u'].sum(),
            g.nodes['n1'].data['xyz'],
            create_graph=True,
            retain_graph=True,
            allow_unused=True,
        )[0]
        du_dx = g.nodes["n1"].data["u_ref_prime"]
        f_loss = torch.nn.MSELoss()(du_dx, du_dx_hat)
        f_vl.append((HARTREE_TO_KCALPERMOL/BOHR_TO_ANGSTROMS) * f_loss.pow(0.5).item())

    for g in ds_te:
        g.nodes["n1"].data["xyz"].requires_grad = True
        e_loss, _ = net(g.heterograph)
        e_te.append(HARTREE_TO_KCALPERMOL * e_loss.pow(0.5).item())
        
        du_dx_hat = torch.autograd.grad(
            g.nodes['g'].data['u'].sum(),
            g.nodes['n1'].data['xyz'],
            create_graph=True,
            retain_graph=True,
            allow_unused=True,
        )[0]
        du_dx = g.nodes["n1"].data["u_ref_prime"]
        f_loss = torch.nn.MSELoss()(du_dx, du_dx_hat)
        f_te.append((HARTREE_TO_KCALPERMOL/BOHR_TO_ANGSTROMS) * f_loss.pow(0.5).item())

    # Energy and Force
    e_vl = np.array(e_vl).mean()
    e_te = np.array(e_te).mean()
    f_vl = np.array(f_vl).mean()
    f_te = np.array(f_te).mean()
    vl = e_vl + f_vl
    te = e_te + f_te

    with open('rmse.csv', 'a') as wf:
        wf.write(f'{epoch}\t{vl:.4f}\t{e_vl:.4f}\t{f_vl:.4f}\t{te:.4f}\t{e_te:.4f}\t{f_te:.4f}\n')


@click.command()
@click.option('--train_ratio', default='0.8', help='Train ratio', type=float)
@click.option('--path_to_data', required=True, help='Path to QC datasets')
@click.option('--index', required=True, help='Index number to define the epoch (<index> * model.checkpoint_frequency)', type=int)
@click.option('--path_to_checkpoints', default='../checkpoints', help='Path to checkpoints')
@click.option('--path_to_config', default='../', help='Path to config')
def cli(**kwargs):
    print(espfit.__version__)
    run(kwargs)


if __name__ == '__main__':
    cli()