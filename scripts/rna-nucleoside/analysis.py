#!/usr/bin/env python
# coding: utf-8
"""
Analyze nucleosides.

Reference
---------
Reparameterization of RNA χ Torsion Parameters for the AMBER Force Field and Comparison to NMR Spectra for Cytidine and Uridine, JCTC, 2010
doi: 10.1021/ct900604a
"""
import click
import mdtraj
import logging
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
#from matplotlib.gridspec import GridSpec
import barnaba as bb
#from barnaba.nucleic import Nucleic


# =====================
# DEFINE
# =====================
logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore")

pd.options.display.max_rows = None
pd.options.display.max_columns = None
pd.options.display.precision = 1
pd.options.display.float_format = '{:.1f}'.format

params = {'legend.fontsize': 20, 
          'font.size': 20, 
          'axes.labelsize': 28,
          'axes.titlesize': 28,
          'xtick.labelsize': 20,
          'ytick.labelsize': 20,
          'savefig.dpi': 600, 
          'figure.figsize': [16, 8],
          'xtick.major.size': 10,
          'xtick.minor.size': 7,
          'ytick.major.size': 10,
          'ytick.minor.size': 7}

plt.rcParams.update(params)

UNIT_NM_TO_ANGSTROMS = 10
UNIT_PS_TO_NS = 1/1000
STRIDE = 1


def radian_to_degree(a):    
    a[np.where(a<0.0)] += 2.*np.pi
    a *= 180.0/np.pi
    # same as above
    #a = a*(180./np.pi)
    #a[np.where(a<0.0)] += 360
    return a


def load_data(kwargs):
    """ Load trajectory
    """
    pdbfile = kwargs["pdbfile"]
    trajfile = kwargs["trajfile"]
    exclude_first_nsteps = kwargs["exclude_first_nsteps"]

    init_pdb = mdtraj.load(pdbfile)
    atom_indices = init_pdb.topology.select('not (water or symbol Na or symbol Cl or symbol K)')
    init_pdb = init_pdb.atom_slice(atom_indices)    
    traj = mdtraj.load(trajfile, top=init_pdb.topology, stride=STRIDE)
    return traj[exclude_first_nsteps:]


def plot_sugar(pdbfile, traj, suffix):
    """ Generate wheel plot for pucker angles.
    """
    # Get data
    init_angles, res = bb.pucker_angles(pdbfile, topology=pdbfile)
    angles, res = bb.pucker_rao_traj(traj)
    angles = angles.reshape(angles.shape[0], angles.shape[-1])
    
    # Plot
    fig = plt.figure(figsize=(6,6))
    #ax = fig.add_subplot(2, int(len(water_models)/2), i+1, polar=True)
    ax = fig.add_subplot(1,1,1, polar=True)
    ax.scatter(angles[:,0], angles[:,1], s=10, c=np.arange(len(angles)), cmap='Blues')
    ax.scatter(init_angles[0,0,0], init_angles[0,0,1], marker="X", c="orange", edgecolors="black", s=150, linewidths=0.5)

    p3 = np.pi/5
    #ax.text(0.5*p3, 1.6, "C3'-endo", ha='center', fontsize=16, fontweight='bold')
    ax.text(0.5*p3, 1.6, "C3'-endo", ha='center', fontsize=16)
    ax.text(1.3*p3, 1.5, "C4'-exo",  ha='center', fontsize=16)
    ax.text(2.5*p3, 1.5, "O4'-endo", ha='center', fontsize=16)
    ax.text(3.7*p3, 1.5, "C1'-exo",  ha='center', fontsize=16)
    #ax.text(4.5*p3, 1.6, "C2'-endo", ha='center', fontsize=16, fontweight='bold')
    ax.text(4.5*p3, 1.6, "C2'-endo", ha='center', fontsize=16)
    ax.text(5.5*p3, 1.5, "C3'-exo",  ha='center', fontsize=16)
    ax.text(6.5*p3, 1.5, "C4'-endo", ha='center', fontsize=16)
    ax.text(7.5*p3, 1.6, "O4'-exo",  ha='center', fontsize=16)
    ax.text(8.5*p3, 1.5, "C1'-endo", ha='center', fontsize=16)
    ax.text(9.5*p3, 1.5, "C2'-exo",  ha='center', fontsize=16)

    xt = np.arange(0, 2*np.pi, p3)
    ax.set_xticks(xt)
    ax.set_yticks([])
    ax.set_ylim(0, 1.2)
    ax.tick_params(axis='both', labelsize=12)

    plt.tight_layout()
    plt.savefig(f"pucker{suffix}.png")
    plt.close()


def plot_chi(traj, suffix):
    """ Generate histogram plot for chi angle distribution.
    """
    # Get data
    angles, res = bb.backbone_angles_traj(traj)
    angles = radian_to_degree(angles).reshape(angles.shape[0], angles.shape[-1])
    angles = angles[:,-1]   # ['alpha', 'beta', 'gamma', 'delta', 'eps', 'zeta', 'chi']
    
    # Plot
    fig = plt.figure(figsize=(8,4))
    gs = fig.add_gridspec(1, 2,  width_ratios=(3, 1), wspace=0.1, bottom=0.2)

    x = np.arange(len(angles))
    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(x, angles, s=20)
    ax.set_xlabel("Time (ns)", fontsize=18)
    ax.set_ylabel(r"$\chi$ angle", fontsize=18)
    ax.set_xlim(0, len(x))
    ax.set_ylim(0, 360)
    ax.yaxis.set_ticks(np.arange(0, 361, 60))
    ax.tick_params(axis='both', labelsize=14)

    bins = np.arange(0, 361, 5)
    ax_hist = fig.add_subplot(gs[0, 1], sharey=ax)
    ax_hist.hist(angles, bins=bins, orientation = 'horizontal', density=True)
    ax_hist.set_xlabel("Fraction", fontsize=18)
    ax_hist.xaxis.set_ticks(np.arange(0.01, 0.031, 0.01))
    ax_hist.tick_params(axis='both', labelsize=14)
    
    plt.setp(ax_hist.get_yticklabels(), visible=False)
    plt.tight_layout()
    plt.savefig(f"chi{suffix}.png")
    plt.close()


def to_csv(traj, suffix):
    """ Compute anti/syn fraction, C3'-endo fraction, and J coupling.
    """
    # J coupling
    couplings, rr = bb.jcouplings_traj(traj, couplings=["H1H2", "H2H3", "H3H4"])
    couplings = couplings.reshape(couplings.shape[0], couplings.shape[-1])    
    j12 = couplings[:,0].mean()
    j12_std = couplings[:,0].std()
    j23 = couplings[:,1].mean()
    j23_std = couplings[:,1].std()
    j34 = couplings[:,-1].mean()
    j34_std = couplings[:,-1].std()
    print(f"# {suffix}")
    print("-------------")
    print(f"J12: {j12:.2f}+-{j12_std:.2f}")
    print(f"J23: {j23:.2f}+-{j23_std:.2f}")
    print(f"J34: {j34:.2f}+-{j34_std:.2f}")
    
    # RNA backbone: Consensus all-angle conformers and modular string nomenclature (an RNA Ontology Consortium contribution), RNA, 2008
    # doi: 10.1261/rna.657708

    # Anti/syn fraction
    angles_b, rr = bb.backbone_angles_traj(traj, angles=["chi", "delta"])
    angles_b = radian_to_degree(angles_b)

    syn, anti, high_anti, not_syn_anti = 0, 0, 0, 0
    for angle in angles_b:
        chi = angle[0][0]
        if chi <= 120:
            syn += 1
        elif chi >= 180 and chi < 240:
            anti += 1
        elif chi >= 240:
            high_anti += 1
        else:
            not_syn_anti += 1
    n = syn + anti + high_anti + not_syn_anti    
    syn = 100 * syn/n
    anti = 100 * anti/n
    high_anti = 100 * high_anti/n
    not_syn_anti = 100 * not_syn_anti/n
    print(f"syn: {syn:.2f}")
    print(f"anti: {anti:.2f}")
    print(f"high_anti: {high_anti:.2f}")
    print(f"not_syn_anti {not_syn_anti:.2f}")
    # C2, C3 endo population
    angles_p, rr = bb.pucker_rao_traj(traj)
    angles_p = radian_to_degree(angles_p)
    c3, c2, not_c2_c3 = 0, 0, 0
    for angle_b, angle_p in zip(angles_b, angles_p):
        delta = angle_b[0][1]
        phase = angle_p[0][0]
        if (delta >= 55 and delta < 110) or (phase >=0 and phase < 36):
            c3 += 1
        elif (delta >= 120 and delta < 175) or (phase >= 144 and phase < 180):
            c2 += 1
        else:
            not_c2_c3 += 1
    n = c3 + c2 + not_c2_c3
    c3 = 100 * c3/n
    c2 = 100 * c2/n
    not_c2_c3 = 100 * not_c2_c3/n
    print(f"c3: {c3:.2f}")
    print(f"c2: {c2:.2f}")
    print(f"not_c2_c3: {not_c2_c3:.2f}\n")

    # TWO END STATE MODEL
    #print("(TWO END STATE MODEL)")
    #if nucleoside.upper().startswith("C") or nucleoside.upper().startswith("U"):
    #    atom_indices = init_traj.topology.select("name H5 or name H6")   # H5 should be H1'
    #    atom_indices = atom_indices.reshape(1,2)
    #    d = mdtraj.compute_distances(traj, atom_indices) * UNIT_NM_TO_ANGSTROMS
    #    d_anti = 3.48
    #    d_syn = 2.12
    #    a = d_syn**6 - d_anti**6
    #    b = ((d_syn * d_anti)/d.mean())**6
    #    anti = 100 * (1/a) * (b - d_anti**6)
    #    syn = 100 - anti
    #    print("anti: {:.2f}".format(anti))
    #    print("syn: {:.2f}".format(syn))
    #c3 = 100*(j34/(j12+j34))
    #c2 = 100 - c3
    #print("C3'-endo: {:.2f}".format(c3))
    #print("C2'-endo: {:.2f}\n".format(c2))

    # Export
    result_dict = {}
    result_dict[suffix] = { "anti-all": anti + high_anti, "syn": syn, "anti": anti, "high-anti": high_anti, "non-syn-anti": not_syn_anti, \
                                "j12": j12, "j12_std": j12_std, "j12_err": j12_std/np.sqrt(len(couplings[:,0])), \
                                "j23": j23, "j23_std": j23_std, "j23_err": j23_std/np.sqrt(len(couplings[:,0])), \
                                "j34": j34, "j34_std": j34_std, "j34_err": j34_std/np.sqrt(len(couplings[:,0])), \
                                "c3": c3, "c2": c2, "non-c2-c3": not_c2_c3
                            }
    import pandas as pd
    df = pd.DataFrame.from_dict(result_dict)
    df.to_csv(f'summary{suffix}.txt', sep='\t', float_format='%.2f')


def report(kwargs):
    # Parameters
    pdbfile = kwargs["pdbfile"]
    suffix = kwargs["suffix"]
    if suffix is not None:
        suffix = "_" + suffix
    # Load trajectory
    traj = load_data(kwargs)
    # Wheel plot for pucker angles
    plot_sugar(pdbfile, traj, suffix)
    # Plot anti/chi angles population
    plot_chi(traj, suffix)
    # Export
    to_csv(traj, suffix)


@click.command()
@click.option('--suffix', default=None, help='Suffix to add to exported files', type=str)
@click.option('--pdbfile', default='solvated.pdb', help='Initial pdb file', type=str)
@click.option('--trajfile', default='traj.nc', help='Trajectory file', type=str)
@click.option('--exclude_first_nsteps', default=0, help='Exclude first n steps', type=int)
def cli(**kwargs):
    report(kwargs)


if __name__ == "__main__":
    cli()