import sys, os, ast, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
sys.path.append(os.path.dirname(__file__)+'/../')

from pH_T_LineGenerator import  pH_T_LineGenerator as pHT_LG
from DataFrameFetcher import DataFrameFetcher as DFF

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(18,6), sharex=True, sharey=True)


""" Begin code from Higgins et al., 2021. Plot only int pH values """

mpl.rcParams['contour.negative_linestyle'] = 'solid'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cmr10'
mpl.rcParams['font.serif'] = 'Times.tcc'
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rc('axes', unicode_minus=False)
mpl.rcParams['hatch.linewidth'] = 2.0
mpl.rcParams['font.size'] = 14


pHnames = ['7.0','8.0','9.0','10.0','11.0','12.0']
Tfloats = np.linspace(273.15, 473.15, num=21)

nomcols = plt.get_cmap('Set2', 6)
cmaplist = [nomcols(i) for i in range(nomcols.N)]
cmaplist[4] = 'chocolate'
cmaplist[3] = 'yellowgreen'


for pH, c in zip(pHnames, cmaplist):

    df=pd.read_csv('../data/speciation/Higgins2021supplement/nominalCO2/pH'+pH+'.csv', sep=',')
    df2=pd.read_csv('../data/speciation/Higgins2021supplement/highsalt/pH'+pH+'.csv', sep=',')
    df3=pd.read_csv('../data/speciation/Higgins2021supplement/lowsalt/pH'+pH+'.csv', sep=',')

    ax[0].plot(Tfloats, df['a_CO2'][:len(Tfloats)], c=c)
    ax[0].fill_between(Tfloats, df2['a_CO2'][:len(Tfloats)], df3['a_CO2'][:len(Tfloats)], facecolor=c, alpha=0.5)

    ax[0].text(270, df['a_CO2'][0], 'pH '+pH, horizontalalignment='right',
    verticalalignment='center', fontsize=12, color=c)


""" end code from Higgins et al., 2021 """

ax[0].set_title('Higgins et al. (2021) model')

# pHT_LG needs somewhere to plot, use this as a burner figure
# we will manually plot what we want in the existing figure defined above
_nofig, _noaxs = plt.subplots(nrows=1)

EMsalts = [np.float64(0.05), np.float64(0.1), np.float64(0.2)]
datadir = 'spec_T_273-473_pH_7-12'

Ts = np.linspace(273.15, 473.15, num=41)

for pH, c in zip(pHnames, cmaplist):

    SUPCRT_salt_vals = []
    for salt in EMsalts:
        # setup DataFrameFetcher at single pH and across T range.
        _DFF = DFF(model='SUPCRTnoGases', salts=[salt], fixGases=False,
          P=1, pHs=[np.float64(pH)], dr=datadir, Ts = Ts)

        LG = pHT_LG(_DFF)

        _noaxs, _df = LG.plot_spec_lines(_noaxs, 'aCO2(aq)', salt_lvl=salt, return_df=True)
        SUPCRT_salt_vals.append(_df['aCO2(aq)'].tolist())

    # now plot and fill in the space between the upper and lower salinities.
    ax[1].plot(Ts, SUPCRT_salt_vals[1], c=c)
    ax[1].fill_between(Ts, SUPCRT_salt_vals[0], SUPCRT_salt_vals[2], facecolor=c, alpha=0.5)
    ax[1].set_title('This work: SUPCRT model')
    ax[1].text(270, SUPCRT_salt_vals[1][0], 'pH '+pH, horizontalalignment='right',
    verticalalignment='center', fontsize=12, color=c)

    pitzer_salt_vals = []
    for salt in EMsalts:
        _DFF = DFF(model='pitzerPHREEQCnoGases', salts=[salt], fixGases=False,
          P=1, pHs=[np.float64(pH)], dr=datadir, Ts = Ts)

        LG = pHT_LG(_DFF)

        _noaxs, _df = LG.plot_spec_lines(_noaxs, 'aCO2(aq)', salt_lvl=salt, return_df=True)
        pitzer_salt_vals.append(_df['aCO2(aq)'].tolist())

    ax[2].plot(Ts, pitzer_salt_vals[1], c=c)
    ax[2].fill_between(Ts, pitzer_salt_vals[0], pitzer_salt_vals[2], facecolor=c, alpha=0.5)
    ax[2].set_title('This work: Pitzer model')
    ax[2].text(270, pitzer_salt_vals[1][0], 'pH '+pH, horizontalalignment='right',
    verticalalignment='center', fontsize=12, color=c)

for a in ax:
    a.set_yscale('log')
    a.set_ylabel('CO$_{2}$ Activity')
    a.yaxis.set_tick_params(labelleft=True)
    a.set_xlabel('Seawater temperature [K]')
    a.set_ylim(1e-10,1e-1)
    a.set_xlim(250,473)

fig.tight_layout()
fig.savefig('figures/S8_H21comparison_473K.pdf')
