#!/usr/bin/env python
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


class GetForce(torch.nn.Module):
    def forward(self, g):
        g.nodes['n1'].data['u_prime'] = torch.autograd.grad(
            g.nodes['g'].data['u'].sum(),
            g.nodes['n1'].data['xyz'],
            create_graph=True,
            retain_graph=True,
            allow_unused=True,
        )[0]

 
def load_dataset(train_ratio, path_to_data, category):
    """Load dgl graphs"""

    # Load dataset
    print(f'Get training dataset')
    ds = CustomGraphDataset.load(f'{path_to_data}/{category}')
    ds.apply(regenerate_impropers, in_place=True)
    ds.apply(add_grad, in_place=True)
    ds.compute_relative_energy()

    ds_tr, ds_vl_te = ds.split([train_ratio, 1-train_ratio])
    ds_vl, ds_te = ds_vl_te.split([0.5, 0.5])    
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


def inspect_rmse(ds, forcefields, suffix):
    """
    """
    # initialize
    mydict = {"u_qm": [], "u_qm_prime": []}
    for forcefield in forcefields:
        mydict["u_%s" % forcefield] = []
        mydict["u_%s_prime" % forcefield] = []
    # dataframe
    df = pd.DataFrame(columns=["SMILES"] + [forcefield + "_ENERGY_RMSE" for forcefield in forcefields] + [forcefield + "_FORCE_RMSE" for forcefield in forcefields])
    # loop over molecule
    for g in ds:
        row = {}
        row["SMILES"] = g.mol.to_smiles()
        # center mean
        u_qm = (g.nodes['g'].data['u_qm'] - g.nodes['g'].data['u_qm'].mean(dim=-1, keepdims=True)).detach().cpu().flatten()
        u_qm_prime = g.nodes['n1'].data['u_qm_prime'].detach().cpu().flatten()
        # append        
        mydict["u_qm"].append(u_qm)
        mydict["u_qm_prime"].append(u_qm_prime)
        for forcefield in forcefields:
            # center mean
            u = (g.nodes['g'].data['u_%s' % forcefield] - g.nodes['g'].data['u_%s' % forcefield].mean(dim=-1, keepdims=True)).detach().cpu().flatten()
            u_prime = g.nodes['n1'].data['u_%s_prime' % forcefield].detach().cpu().flatten()
            # rmse
            e_rmse = esp.metrics.rmse(u_qm, u) * HARTREE_TO_KCALPERMOL
            f_rmse = esp.metrics.rmse(u_qm_prime, u_prime) * (HARTREE_TO_KCALPERMOL/BOHR_TO_ANGSTROMS)
            # dataframe
            row[forcefield + "_ENERGY_RMSE"] = e_rmse.item()
            row[forcefield + "_FORCE_RMSE"] = f_rmse.item()
            # mydict
            mydict["u_%s" % forcefield].append(u)
            mydict["u_%s_prime" % forcefield].append(u_prime)
            #print(forcefield, u_qm_prime.shape, u_prime.shape)
        df.loc[len(df)] = row
    df.to_csv(f"inspect_{suffix}.csv")

    return mydict


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
    category = kwargs['category']
    path_to_data = kwargs['path_to_data']
    ckptfile = kwargs['ckptfile']
    configfile = kwargs['configfile']
    forcefields = kwargs['forcefields']

    # Load dataset
    ds_tr, ds_vl, ds_te = load_dataset(train_ratio, path_to_data, category)

    # Load model
    model = EspalomaModel.from_toml(configfile)
    net = update_net(model.net)
    model.net = net
    state_dict = torch.load(ckptfile, map_location=torch.device('cpu'))
    net.load_state_dict(state_dict)
    net.eval()

    for suffix in ['tr', 'vl', 'te']:
        if suffix == 'tr':
            ds = ds_tr
            del ds_tr
        elif suffix == 'vl':
            ds = ds_vl
            del ds_vl
        else:
            ds = ds_te
            del ds_te

        mydict = inspect_rmse(ds, forcefields, suffix)
        del ds

        # Report summary
        with open(f"summary_{suffix}.csv", "w") as wf:
            wf.write("# energy / force\n")
            for forcefield in forcefields:
                ci_e, ci_f = bootstrap_mol(mydict["u_qm"], mydict["u_%s" % forcefield], mydict["u_qm_prime"], mydict["u_%s_prime" % forcefield])
                wf.write(f"{forcefield}: {ci_e} / {ci_f}\n")


@click.command()
@click.option("--train_ratio", required=True, type=float, help="train ratio")
@click.option("--ckptfile", required=True, help="Checkpoint file (e.g. ckpt.pt)", type=str)
@click.option("--configfile", required=True, help="Config file (e.g. config.toml)", type=str)
@click.option('--path_to_data', required=True, help='Path to QC datasets')
@click.option("--category", required=True, type=click.Choice(['debug', 'small', 'peptide', 'rna']), help="category of the dataset")
@click.option("--forcefields", required=True, multiple=True, type=click.Choice(['amber14sb', 'openff-2.1.0', 'openff-2.0.0', 'gaff-2.11']), help="Baseline forcefields")
def cli(**kwargs):
    print(espfit.__version__)
    run(kwargs)


if __name__ == "__main__":
    cli()
