#!/usr/bin/env python
import os
import click
import torch
import pandas as pd
import numpy as np
import espaloma as esp
import espfit
from espaloma.graphs.utils.regenerate_impropers import regenerate_impropers
from espfit.utils.graphs import CustomGraphDataset
from espfit.app.train import EspalomaModel
from espfit.utils.units import HARTREE_TO_KCALPERMOL, BOHR_TO_ANGSTROMS 
from rdkit.Chem import PandasTools



class GetForce(torch.nn.Module):
    def forward(self, g):
        g.nodes['n1'].data['u_prime'] = torch.autograd.grad(
            g.nodes['g'].data['u'].sum(),
            g.nodes['n1'].data['xyz'],
            create_graph=True,
            retain_graph=True,
            allow_unused=True,
        )[0]

 
def load_dataset(train_ratio, path_to_data, dataset):
    """Load dgl graphs"""

    # Load dataset
    print(f'Get training dataset')

    if dataset == 'rna-diverse':
        ds = CustomGraphDataset.load(f'{path_to_data}/{dataset}')
        ds.apply(regenerate_impropers, in_place=True)
        ds.apply(add_grad, in_place=True)
        ds.compute_relative_energy()
        ds_tr, ds_vl_te = ds.split([train_ratio, 1-train_ratio])
        ds_vl, ds_te = ds_vl_te.split([0.5, 0.5])    
        print("Train:Validate:Test = {}:{}:{}".format(len(ds_tr),len(ds_vl),len(ds_te)))

    elif dataset == 'rna-nucleoside':
        ds = CustomGraphDataset.load(f'{path_to_data}/{dataset}')
        ds.apply(regenerate_impropers, in_place=True)
        ds.apply(add_grad, in_place=True)
        ds.compute_relative_energy()
        ds_tr, ds_vl = [], []
        ds_te = ds
        print("Train:Validate:Test = {}:{}:{}".format(len(ds_tr),len(ds_vl),len(ds_te)))

    elif dataset == 'rna-trinucleotide':
        ds = CustomGraphDataset.load(f'{path_to_data}/{dataset}')
        ds.apply(regenerate_impropers, in_place=True)
        ds.apply(add_grad, in_place=True)
        ds.compute_relative_energy()
        ds_tr = ds
        ds_vl, ds_te = [], []
        print("Train:Validate:Test = {}:{}:{}".format(len(ds_tr),len(ds_vl),len(ds_te)))

    return ds_tr, ds_vl, ds_te


def add_grad(g):
    g.nodes["n1"].data["xyz"].requires_grad = True
    return g


def _bootstrap_mol(x, y, n_samples=1000, ci=0.95):
    z = []
    for _x, _y in zip(x, y):
        mse = torch.nn.functional.mse_loss(_x, _y).item()
        z.append(np.sqrt(mse))
    z = np.array(z)

    results = []
    for _ in range(n_samples):
        _z = np.random.choice(z, z.size, replace=True)
        results.append(_z.mean())

    results = np.array(results)
    low = np.percentile(results, 100.0 * 0.5 * (1 - ci))
    high = np.percentile(results, (1 - ((1 - ci) * 0.5)) * 100.0)
    mean = z.mean()

    return mean, low, high


def bootstrap_mol(u_ref, u, u_ref_prime, u_prime):
    """
    Bootstrap over molecules
    """
    mean, low, high = _bootstrap_mol(u_ref, u)
    ci_e = esp.metrics.latex_format_ci(
        mean * HARTREE_TO_KCALPERMOL, 
        low * HARTREE_TO_KCALPERMOL, 
        high * HARTREE_TO_KCALPERMOL
        )
    mean, low, high = _bootstrap_mol(u_ref_prime, u_prime)
    ci_f = esp.metrics.latex_format_ci(
        mean * (HARTREE_TO_KCALPERMOL/BOHR_TO_ANGSTROMS), 
        low * (HARTREE_TO_KCALPERMOL/BOHR_TO_ANGSTROMS), 
        high * (HARTREE_TO_KCALPERMOL/BOHR_TO_ANGSTROMS)
        )
    return ci_e, ci_f


def calc_metric(net, ds, dataset, suffix):
    """
    """
    u, u_ref = [], []
    u_prime, u_ref_prime = [], []
    df = pd.DataFrame(columns=["SMILES", "RMSE_ENERGY", "RMSE_FORCE", "n_snapshots"])

    for g in ds:
        # Compute energy and force
        #g.heterograph = g.heterograph.to("cuda:0")
        net(g.heterograph)
        # Energy
        _u = (g.nodes['g'].data['u'] - g.nodes['g'].data['u'].mean(dim=-1, keepdims=True)).detach().cpu().flatten()
        _u_ref = g.nodes['g'].data['u_ref_relative'].detach().cpu().flatten()
        u.append(_u)
        u_ref.append(_u_ref)
        # Force
        _u_prime = g.nodes['n1'].data['u_prime'].detach().cpu().flatten()
        _u_ref_prime = g.nodes['n1'].data['u_ref_prime'].detach().cpu().flatten()
        u_prime.append(_u_prime)
        u_ref_prime.append(_u_ref_prime)
        # Smiles
        smi = g.mol.to_smiles()
        # Pandas
        df.loc[len(df)] = {
            'SMILES': smi,
            'RMSE_ENERGY': esp.metrics.rmse(_u_ref, _u).item() * HARTREE_TO_KCALPERMOL,
            'RMSE_FORCE': esp.metrics.rmse(_u_ref_prime, _u_prime).item() * (HARTREE_TO_KCALPERMOL/BOHR_TO_ANGSTROMS),
            'n_snapshots': g.nodes['n1'].data['xyz'].shape[1]
        }

    # Export csv and html
    df = df.sort_values(by="RMSE_FORCE", ascending=False)
    df.to_csv(f'rmse_{suffix}_{dataset}.csv', sep='\t')
    PandasTools.AddMoleculeColumnToFrame(df, "SMILES", "MOL")
    open(f"rmse_{suffix}_{dataset}.html", "w").write(df.to_html())
    # Bootstrap over molecule
    ci_e, ci_f = bootstrap_mol(u_ref, u, u_ref_prime, u_prime)
    # Report summary
    ofile = f'summary_{dataset}.csv'
    if os.path.exists(ofile):
        wf = open(ofile, 'a')
    else:
        wf = open(ofile, 'w')
    wf.write(f">{dataset} ({suffix})\n")
    wf.write("----------\n")
    wf.write(f"energy: {ci_e}\n")
    wf.write(f"force: {ci_f}\n")
    wf.write("\n")
    wf.close()


def update_net(net):
    """Update net model"""
    modules = []
    for _net in net[:-1]:
        modules.append(_net)
    modules.append(GetForce())
    return torch.nn.Sequential(*modules)


def run(kwargs):
    """
    """
    train_ratio = kwargs['train_ratio']
    dataset = kwargs['dataset']
    path_to_data = kwargs['path_to_data']
    ckptfile = kwargs['ckptfile']
    configfile = kwargs['configfile']

    # Load dataset
    ds_tr, ds_vl, ds_te = load_dataset(train_ratio, path_to_data, dataset)

    # Load model
    model = EspalomaModel.from_toml(configfile)
    net = update_net(model.net)
    model.net = net
    state_dict = torch.load(ckptfile, map_location=torch.device('cpu'))
    net.load_state_dict(state_dict)
    net.eval()

    for suffix in ['tr', 'vl', 'te']:
        if suffix == 'tr' and ds_tr:
            calc_metric(net, ds_tr, dataset, suffix)
            del ds_tr
        elif suffix == 'vl' and ds_vl:
            calc_metric(net, ds_vl, dataset, suffix)
            del ds_vl
        elif suffix == 'te' and ds_te:
            calc_metric(net, ds_te, dataset, suffix)
            del ds_te
            

@click.command()
@click.option("--train_ratio", required=True, type=float, help="train ratio")
@click.option("--ckptfile", required=True, help="best model (e.g. ckpt.pt)", type=str)
@click.option("--configfile", required=True, help="config file (e.g. config.toml)", type=str)
@click.option('--path_to_data', required=True, help='Path to QC datasets')
@click.option("--dataset", required=True, type=click.Choice(['rna-diverse', 'rna-trinucleotide', 'rna-nucleoside']), help="RNA dataset")
def cli(**kwargs):
    print(espfit.__version__)
    run(kwargs)


if __name__ == "__main__":
    cli()
