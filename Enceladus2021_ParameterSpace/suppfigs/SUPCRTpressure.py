import sys, os, ast, math
sys.path.append(os.path.dirname(__file__)+'../../../NutMEG')
import matplotlib.pyplot as plt
import matplotlib as mpl

import pandas as pd
import numpy as np
import NutMEG as nm
from NutMEG.reactor.saved_systems.Enceladus import Enceladus as Enc
from NutMEG.culture.saved_organisms.TypicalOptimalMethanogen import TypicalOptimalMethanogen as TOM

mpl.rcParams['contour.negative_linestyle'] = 'solid'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cmr10'
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rc('axes', unicode_minus=False)
mpl.rcParams['hatch.linewidth'] = 2.0
# print(mpl.rcParams.keys()) # to see what can be changed
mpl.rcParams['font.size'] = 14

## set up reaction objects for the equilibrium constants

envi = nm.environment(T=273.15, P=1e5)
CO2 = nm.reaction.reagent('CO2(aq)', envi, phase='aq')
H2O = nm.reaction.reagent('H2O(l)', envi, phase='l')
HCO3 = nm.reaction.reagent('HCO3-', envi, phase='aq')
CO3 = nm.reaction.reagent('CO3--', envi, phase='aq')
H = nm.reaction.reagent('H+', envi, phase='aq')

r1 = nm.reaction.reaction(reactants={CO2:1, H2O:1}, products={HCO3:1, H:1}, env=envi)
r2 = nm.reaction.reaction(reactants={CO2:1, H2O:1}, products={CO3:1, H:2}, env=envi)

def get_rate_consts(T, P):
    r1.env.T = T
    r1.env.P = P
    r2.env.T = T
    r2.env.P = P

    r1.rto_current_env()
    r2.rto_current_env()
    return math.exp(r1.lnK), math.exp(r2.lnK)

fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(14,15), constrained_layout=True)
for a in ax.flatten():
    a.set_xlim(240, 475)
    a.set_xlabel('Seawater temperature [K]')


pHnames = ['7.0','8.0','9.0','10.0','11.0','12.0']
Tfloats = np.linspace(273.15, 473.15, num=21)

nomcols = plt.get_cmap('Set2', 6)
cmaplist = [nomcols(i) for i in range(nomcols.N)]
cmaplist[4] = 'chocolate'
cmaplist[3] = 'yellowgreen'

K1s_1, K2s_1 = [],[]
K1s_100, K2s_100 = [],[]
for T in Tfloats:
    _K1_1, _K2_1 =  get_rate_consts(T, 1e5)
    _K1_100, _K2_100 =  get_rate_consts(T, 1e7)
    K1s_1.append(_K1_1)
    K2s_1.append(_K2_1)
    K1s_100.append(_K1_100)
    K2s_100.append(_K2_100)

for pH, c in zip(pHnames, cmaplist):
    # Hconc = 10**(-float(pH))

    _df=pd.read_csv('../E21data/Speciation/nominalCO2/pH'+pH+'.csv', sep=',')

    _CO2_1, _HCO3_1, _CO3_1 = np.zeros(len(Tfloats)), np.zeros(len(Tfloats)), np.zeros(len(Tfloats))
    _CO2_100, _HCO3_100, _CO3_100 = np.zeros(len(Tfloats)), np.zeros(len(Tfloats)), np.zeros(len(Tfloats))


    ii = 0
    for T, K1_1, K2_1, K1_100, K2_100 in zip(Tfloats, K1s_1, K2s_1, K1s_100, K2s_100):
        Hconc = 10**(-_df['pH'][ii])
        thisCO2_1 = (1+(K1_1/Hconc)+(K2_1/(Hconc*Hconc)))**(-1)
        thisCO2_100 = (1+(K1_100/Hconc)+(K2_100/(Hconc**2)))**(-1)

        _CO2_1[ii] = thisCO2_1
        _HCO3_1[ii] = (K1_1*thisCO2_1/Hconc)
        _CO3_1[ii] = (K2_1*thisCO2_1/(Hconc**2))

        _CO2_100[ii] = (thisCO2_100)
        _HCO3_100[ii] = (K1_100*thisCO2_100/Hconc)
        _CO3_100[ii] = (K2_100*thisCO2_100/(Hconc**2))

        ii += 1

    CO2err = 100*np.abs(_CO2_100-_CO2_1)/_CO2_100
    ax[0][0].plot(Tfloats, _CO2_1, c=c)
    ax[0][0].plot(Tfloats, _CO2_100, c=c, ls='--')
    ax[0][1].plot(Tfloats, CO2err, c=c)
    ax[0][0].text(270, _CO2_1[0], 'pH '+pH, horizontalalignment='right',
    verticalalignment='center', fontsize=12, color=c)
    ax[0][1].text(270, CO2err[0], 'pH '+pH, horizontalalignment='right',
    verticalalignment='center', fontsize=12, color=c)

    HCO3err = 100*np.abs(_HCO3_100-_HCO3_1)/_HCO3_100
    ax[1][0].plot(Tfloats, _HCO3_1, c=c)
    ax[1][0].plot(Tfloats, _HCO3_100, c=c, ls='--')
    ax[1][1].plot(Tfloats, HCO3err, c=c)
    if pH != '9.0' and pH != '8.0':
        ax[1][0].text(270, _HCO3_1[0], 'pH '+pH, horizontalalignment='right',
          verticalalignment='center', fontsize=12, color=c)
        if pH == '7.0':
            ax[1][1].text(270, 2, 'pH '+pH, horizontalalignment='right',
              verticalalignment='center', fontsize=12, color=c)
        else:
            ax[1][1].text(270, HCO3err[0], 'pH '+pH, horizontalalignment='right',
              verticalalignment='center', fontsize=12, color=c)
    elif pH == '9.0':
        ax[1][0].text(270, 1, 'pH '+pH, horizontalalignment='right',
          verticalalignment='center', fontsize=12, color=c)
        ax[1][1].text(270, -0.1, 'pH '+pH, horizontalalignment='right',
          verticalalignment='center', fontsize=12, color=c)
    elif pH == '8.0':
        ax[1][0].text(270, 0.91, 'pH '+pH, horizontalalignment='right',
          verticalalignment='center', fontsize=12, color=c)
        ax[1][1].text(270, 0.3, 'pH '+pH, horizontalalignment='right',
          verticalalignment='center', fontsize=12, color=c)

    # ax[1][1].text(270, HCO3err[0], 'pH '+pH, horizontalalignment='right',
      # verticalalignment='center', fontsize=12, color=c)

    CO3err = 100*np.abs(_CO3_100-_CO3_1)/_CO3_100
    ax[2][0].plot(Tfloats, _CO3_1, c=c)
    ax[2][0].plot(Tfloats, _CO3_100, c=c, ls='--')
    ax[2][1].plot(Tfloats, CO3err, c=c)
    ax[2][0].text(270, _CO3_1[0], 'pH '+pH, horizontalalignment='right',
    verticalalignment='center', fontsize=12, color=c)
    ax[2][1].text(270, CO3err[0], 'pH '+pH, horizontalalignment='right',
    verticalalignment='center', fontsize=12, color=c)

#axis housekeeping
ax[0][0].set_ylabel('[CO$_{2}$] / [DIC]')
ax[0][1].set_ylabel('[CO$_{2}$] / [DIC] % error (1-100 bar)')
ax[1][0].set_ylabel('[HCO$_{3}^{-}$] / [DIC]')
ax[1][1].set_ylabel('[HCO$_{3}^{-}$] / [DIC] % error (1-100 bar)')
ax[2][0].set_ylabel('[CO$_{3}^{2-}$] / [DIC]')
ax[2][1].set_ylabel('[CO$_{3}^{2-}$] / [DIC] % error (1-100 bar)')

for i, aa in enumerate(ax):
    for j, a in enumerate(aa):
        if j == 0:
            a.set_yscale('log')

# plt.suptitle('Using [H] from speciation (changing with T)')
# plt.suptitle('Using constant [H] with T')

plt.savefig('SUPCRTconcs.pdf')
