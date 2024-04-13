#!/usr/bin/env python
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator


# Settings
params = {'legend.fontsize': 40, 
          'font.size': 40, 
          'axes.labelsize': 48,
          'axes.titlesize': 48,
          'xtick.labelsize': 40,
          'ytick.labelsize': 40,
          'savefig.dpi': 600, 
          'figure.figsize': [64, 8],
          'xtick.major.size': 10,
          'xtick.minor.size': 7,
          'ytick.major.size': 10,
          'ytick.minor.size': 7}

plt.rcParams.update(params)

import matplotlib.colors as mcolors
mycolors = mcolors.TABLEAU_COLORS


def early_stopping(df):
    """Early stopping w/ joint loss"""
    print(">Early stopping")

    vl = df['Validation_Total']
    e_vl = df['Validation_Energy']
    f_vl = df['Validation_Force']

    min_epochs = 80
    limit = 5
    count = 0
    index = []
    loss_min = vl[min_epochs:].max()
    for i in df.index[min_epochs:]:
        if count == limit:
            break
        loss_current = vl.iloc[i]
        # reset counter
        if loss_current < loss_min:
            count = 0
            index = [{"Epoch": int(df.iloc[i]["Epoch"]), "val_loss": loss_current, "e_loss": e_vl.iloc[i], "f_loss": f_vl.iloc[i]}]
            loss_min = loss_current
        # update counter
        else:
            index.append({"Epoch": int(df.iloc[i]["Epoch"]), "val_loss": loss_current, "e_loss": e_vl.iloc[i], "f_loss": f_vl.iloc[i]})
            count += 1

    print(index)
    epoch_es = index[0]['Epoch']
    print(f">Epochs_es: joint_es={epoch_es}")

    source = f"../checkpoints/ckpt{epoch_es}.pt"
    destination = f"./ckpt.pt"
    shutil.copy(source, destination)


def plot(df):
    epochs = df['Epoch']
    vl = df['Validation_Total']
    e_vl = df['Validation_Energy']
    f_vl = df['Validation_Force']

    # y-axis: https://matplotlib.org/3.1.1/gallery/ticks_and_spines/tick-locators.html
    f, ax = plt.subplots(figsize=(20, 10))
    ax.plot(epochs, e_vl, lw=4, label='Energy')
    ax.plot(epochs, f_vl, lw=4, label='Force')
    ax.plot(epochs, vl, lw=4, label='Energy+Force')
    ax.set_xlim(0, epochs.max())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.set_ylim(pow(10,0), pow(15,2))
    ax.set_yscale("log", base=10)
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))    
    ax.set_xlabel("Epochs")
    ax.set_ylabel("Loss")    
    ax.legend(loc="upper right", bbox_to_anchor=(1.00, 1.00), fontsize=28)
    plt.tight_layout()
    plt.savefig("rmse.png")
    plt.close()


def run():
    # Load and reformat
    df = pd.read_csv('rmse.csv', sep='\t')
    df = df.sort_values(by=['Epoch'])
    df = df.reset_index(drop=True)
    early_stopping(df)
    plot(df)


if __name__ == "__main__":
    run()
