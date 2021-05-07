import sys, os, ast, math
sys.path.append(os.path.dirname(__file__)+'../../../NutMEG')
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd

import NutMEG as nm

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
# print(mpl.rcParams.keys()) # to see what can be changed
mpl.rcParams['font.size'] = 14



pHnames = ['7.0','8.0','9.0','10.0','11.0','12.0']
pHnames5 = ['7.5','8.5','9.5','10.5','11.5', 'k']

Tfloats = np.linspace(273.15, 473.15, num=21)
  # pHlist = np.ndarray((len(pHfloats), len(Tfloats)))

fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(14,10))

nomcols = plt.get_cmap('Set2', 6)
cmaplist = [nomcols(i) for i in range(nomcols.N)]
cmaplist[4] = 'chocolate'
cmaplist[3] = 'yellowgreen'


for pH, pH5, c in zip(pHnames, pHnames5, cmaplist):
    # df=pd.read_csv('../HTHeatingdata/pH'+pH+'.csv', sep=',')
    df=pd.read_csv('../E21data/Speciation/nominalCO2/pH'+pH+'.csv', sep=',')
    df2=pd.read_csv('../E21data/Speciation/highsalt/pH'+pH+'.csv', sep=',')
    df3=pd.read_csv('../E21data/Speciation/lowsalt/pH'+pH+'.csv', sep=',')

    ax[0][0].plot(Tfloats, df['a_CO2'], c=c)
    ax[0][0].fill_between(Tfloats, df2['a_CO2'], df3['a_CO2'], facecolor=c, alpha=0.5)
    ax[1][0].plot(Tfloats, df['pH'], c=c)
    ax[1][0].fill_between(Tfloats, df3['pH'], df2['pH'], facecolor=c, alpha=0.5)

    ax[0][0].text(270, df['a_CO2'][0], 'pH '+pH, horizontalalignment='right',
    verticalalignment='center', fontsize=12, color=c)
    ax[1][0].text(270, df['pH'][0], 'pH '+pH, horizontalalignment='right',
    verticalalignment='center', fontsize=12, color=c)

    if pH5 =='k':
        break

    _df=pd.read_csv('../E21data/Speciation/nominalCO2/pH'+pH5+'.csv', sep=',')
    _df2=pd.read_csv('../E21data/Speciation/highsalt/pH'+pH5+'.csv', sep=',')
    _df3=pd.read_csv('../E21data/Speciation/lowsalt/pH'+pH5+'.csv', sep=',')

    ax[0][1].plot(Tfloats, _df['a_CO2'], c=c)
    ax[0][1].fill_between(Tfloats, _df2['a_CO2'], _df3['a_CO2'], facecolor=c, alpha=0.5)

    ax[1][1].plot(Tfloats, _df['pH'], c=c)
    ax[1][1].fill_between(Tfloats, _df2['pH'], _df3['pH'], facecolor=c, alpha=0.5)
    ax[0][1].text(270, _df['a_CO2'][0], 'pH '+pH5, horizontalalignment='right',
    verticalalignment='center', fontsize=12, color=c)
    ax[1][1].text(270, _df['pH'][0], 'pH '+pH5, horizontalalignment='right',
    verticalalignment='center', fontsize=12, color=c)

  # for the label
  # ax.plot([0,0], [0,0], c=color, linestyle='dashed', label='pH when warmed from 273.15 K')

R = nm.environment() # 1bar, 298 K

H2O = nm.reaction.reagent('H2O(l)', R, phase='l')
H = nm.reaction.reagent('H+', R, phase='aq')
OH = nm.reaction.reagent('OH-', R, phase='aq')

rxn = nm.reaction.reaction({H2O:1}, {H:1, OH:1}, R)

Ts = np.linspace(273,500, num=(500-273)*2)
pHs = []
for T in Ts:
 rxn.env.T = T
 rxn.rto_current_env()
 pHs.append(-math.log10(math.sqrt(math.exp(rxn.lnK))))

H2Oc = 'rebeccapurple'
ax[1][0].plot(Ts, pHs, c=H2Oc, linewidth=2, linestyle='dashed')
ax[1][0].text(470, 6.1, 'pH of pure water', horizontalalignment='right',
verticalalignment='center', fontsize=12, color=H2Oc)
ax[1][1].plot(Ts, pHs, c=H2Oc, linewidth=2, linestyle='dashed')
ax[1][1].text(470, 6.1, 'pH of pure water', horizontalalignment='right',
verticalalignment='center', fontsize=12, color=H2Oc)


for i, aa in enumerate(ax):
    for j, a in enumerate(aa):
        a.set_xlabel('Temperature [K]')
        a.set_xlim(240,475)
        if i ==0:
            a.set_yscale('log')
            a.set_ylabel('CO2 Activity')
            a.set_ylim(1e-10,1e-1)
        if i == 1:
            a.set_ylabel('pH')
            a.set_ylim(5.5,12.5)

plt.tight_layout()
plt.savefig('CO2act.pdf')
