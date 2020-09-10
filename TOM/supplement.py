from methanogen_implementer import efficiencies
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import NutMEG as es
import sys, os, math

sys.path.append(os.path.dirname(__file__)+'/../../NutMEG/')
from NutMEG.culture.base_organism.synthesis.cell_synthesis import cell_synthesis as synth

# change matplotlib rcparams updating plot fonts etc.
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cmr10'
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['lines.markersize'] = 8



def concs_vs_t(plot=True, Trange=range(275,375),CH4=3e-8):
    """ plot the concentration of CO2, H2, and CH4 with temperature"""

    fig = plt.figure(figsize = (7,5))
    ax= fig.add_subplot(111)

    CO2lst=[]
    H2lst=[]
    for T in Trange:
        peff = efficiencies.get_eff('averageMethanogen', Temp=T)
        CO2lst.append(peff.params['mol_CO2'])
        H2lst.append(peff.params['mol_H2'])
    if plot:
        ax.plot(Trange, CO2lst, label=r'[$\mathregular{CO}_2$]', linewidth=4)
        ax.plot(Trange, H2lst, label=r'[$\mathregular{H}_2$]', linewidth=4)
        ax.plot(Trange, [CH4]*len(Trange), label=r'[$\mathregular{CH}_4$]', linewidth=4)
        ax.set_ylabel('Concentration [M]', fontsize=14)
        ax.set_xlabel('Temperature [K]', fontsize=14)
        # plt.yscale('log')
        ax.legend()
        plt.tight_layout()
        plt.savefig('concs.pdf')
        # plt.show()
    return CO2lst, H2lst

# concs_vs_t()

def ESynth():
    """Plot energetic cost of synthesis reactions per dry g of cells"""
    AAl = []
    Pl = []
    tot=[]
    for i in range(270,400):
        s = synth.get_ESynth_density(i)

        AAl.append(s["AAsynth"])
        Pl.append(s['Psynth'])
        tot.append(s['Psynth']+s["AAsynth"])

    fig = plt.figure(figsize = (7,5))
    ax= fig.add_subplot(111)

    ax.plot(range(270,400), AAl, label='Amino acid synthesis', c='#2c7bb6', linewidth=4)
    ax.plot(range(270,400), Pl, label='Polymerisation', c='#fdae61', linewidth=4)
    ax.plot(range(270,400), tot, label='Total', c='#d7191c', linewidth=4)


    ax.set_ylabel(r'Energetic cost [J per dry g]', fontsize=14)
    ax.set_xlabel('Temperature [K]', fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=14)


    ax.set_xlim(270, 400)
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig('ESynth.pdf')
    # plt.show()

# ESynth()

def CH4get(growthrate, CH4conc=3e-8):
    """Use reverse Powell method to compute CH4 production rate from growth rate"""
    return CH4conc*growthrate*(math.exp(growthrate) - 1)

def avgCH4_GR():
    """Plot range of CH4 production rates an optimal growth rates for the
    empirical methanogens, and their linear regressions.
    """
    file = 'data/methanogens.csv'
    df = pd.read_csv(file, header=0)
    Ts, CH4s, Pressures, GRs = [],[],[],[]
    for index, row in df.iterrows():
        try:
            Ts.append(273+((float(row['Min. optimal growth temp.'])+float(row['Max. optimal growth temp.'])))/2)
            Pressures.append(float(row['Pressure'])*1000)
            GRs.append(row['Growth rate']/3600)
            CH4s.append(CH4get(row['Growth rate']/3600))
        except Exception as e:
            print(str(e))
            continue

    CH4pol = np.polyfit(Ts, np.log(CH4s), 1)
    GRspol = np.polyfit(Ts, np.log(GRs), 1)

    fig, axs = plt.subplots(nrows=2, figsize=(7,7))
    axs[0].scatter(Ts, CH4s, label='Empirical methanogens', marker='x')
    axs[0].plot([min(Ts), max(Ts)], np.exp([CH4pol[0]*(min(Ts))+CH4pol[1], CH4pol[0]*(max(Ts))+CH4pol[1]]), c='k', label='Typical optimal methanogen', linewidth=3)

    axs[1].scatter(Ts, GRs, label='Empirical methanogens', marker='x')
    axs[1].plot([min(Ts), max(Ts)], np.exp([GRspol[0]*(min(Ts))+GRspol[1], GRspol[0]*(max(Ts))+GRspol[1]]), c='k',  label='Typical optimal methanogen', linewidth=3)

    axs[0].set_yscale('log')
    axs[1].set_yscale('log')

    axs[0].set_ylabel('CH4 uptake rate $[\mathregular{M\ s}^{-1}]$', fontsize=14)
    axs[0].set_ylim(1e-20,1e-13)

    axs[1].set_ylabel('Growth rate $[\mathregular{s}^{-1}]$', fontsize=14)
    axs[1].set_ylim(1e-6,1e-3)
    axs[1].set_xlabel('Temperature [K]', fontsize=14)
    plt.tight_layout()
    plt.savefig('avgCH4GR.pdf')
    # plt.show()

# avgCH4_GR()
