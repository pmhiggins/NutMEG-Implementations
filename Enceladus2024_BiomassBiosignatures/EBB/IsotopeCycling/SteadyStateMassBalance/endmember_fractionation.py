import sys, os; sys.path.append(os.path.dirname(__file__)+'/../../../../../NutMEG')
sys.path.append(os.path.dirname(__file__)+'/../../')

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from copy import deepcopy
from EncBmBs_utils import IsotopeConversions as ICtools
from itertools import chain
from SSMBtest import SimpleMethaneBox


plt.rcParams['errorbar.capsize'] = 8
plt.rcParams['lines.linewidth'] = 3
ls = ['dashed', '-', 'dotted']

# yield of CH4 (frac. of CO2 outbound flux) to be used for the steady
# state modules. Arbitrarily low value used here to avoid using up CO2 stock.
yield_CH4 = 0.01


dC_CO2_space = 60
dC_CO2_sf = dC_CO2_space + 1 + 2 # max enrichment at seafloor from ice sheet + ocean trans.

# maximum and minimum epsilons for CH4 formation on Enceladus
# these are a broad summary of results computed in the IsotopeCycling
# directory and/or described in the manuscript.
epsilon_habitat = [22.68, 4.4]
epsilon_thermo = [24, -5]
epsilon_FTT = [49.5, 15]


# net isotope difference owing to carbonate speciation is similar across
# Enceladus' parameter space as shown in CarbonateIsotopeSpeciation.py.
# Below are characteristic changes at temperatures relevant to Enceladus' interior.
DdC_CO2_DIC = 10. # DIC is +10 pm vs CO2 at T=0.
DdC_CO2_60_0 = 4.
DdC_CO2_100_0 = 7.



def DdC_CH4_x1_x4(dC_CH4_sf):
    # max and min enrichment to CH4 on ocean transit.
    return (dC_CH4_sf - 4 - 1, dC_CH4_sf + 1)


fig, ax = plt.subplots(figsize=(6,4), nrows=1)
axs = [ax]
for ax in axs:
    ax.set_xlim(0,16)
axs[0].set_ylim(0,100)

ax = axs[0]



def get_CH4_minmax(dC_CO2_outs, e_lst):

    CH4_sf_all = []
    CH4_space_all = []
    for dC_CO2_out in dC_CO2_outs:
        for i in range(2):
            SMB = SimpleMethaneBox(ICtools.dCtoR(dC_CO2_out), yield_CH4, ICtools.etoa(e_lst[i]/1000))
            CH4_sf_all.append(ICtools.RtodC(SMB.CH4_SB.df['R']['CH4_out']))
    for CH4 in CH4_sf_all:
        for sp_DdC in DdC_CH4_x1_x4(CH4):
            CH4_space_all.append(sp_DdC)

    return min(CH4_space_all), max(CH4_space_all)


def get_m_eb(ls):
    """ get mean and errorbars from a list for plotting """
    m = (min(ls)+max(ls))/2
    eb = abs(m-min(ls))
    return m, eb

def addlines(ax, x1,x2, CO2lst, CH4lst, c):
    _m, _eb = get_m_eb(CO2lst)
    ax.errorbar([x1], _m, yerr=_eb, c='slategray', markeredgewidth=2, capsize=4, alpha=0.8)
    _m, _eb = get_m_eb(CH4lst)
    _ = ax.errorbar([x2], _m, yerr=_eb, c=c, markeredgewidth=2, alpha=0.9)


#### First endmember: FTT/Sabatier

# relevant CO2: >100C
# space is the min, in case no changes on transit
dC_FTT_CO2 = [dC_CO2_space + DdC_CO2_100_0, dC_CO2_sf + DdC_CO2_DIC]
dC_FTT_CH4 = get_CH4_minmax(dC_FTT_CO2, epsilon_FTT)

addlines(axs[0], 1.5, 2.5, dC_FTT_CO2, dC_FTT_CH4, 'k')

#### Second endmember: Thermogenesis of abiotic OM

dC_abiothermo_CO2 = dC_FTT_CO2
dC_abiothermo_CH4 = get_CH4_minmax(dC_abiothermo_CO2, epsilon_thermo)

addlines(axs[0], 5.5, 6.5, dC_abiothermo_CO2, dC_abiothermo_CH4, 'tab:red')


#### Third endmember: Biotic Methane

dC_bio_CO2 = [dC_CO2_space + DdC_CO2_60_0, dC_CO2_sf + DdC_CO2_60_0]
dC_bio_CH4 = get_CH4_minmax(dC_bio_CO2, epsilon_habitat)

addlines(axs[0], 9.5, 10.5, dC_bio_CO2, dC_bio_CH4, 'tab:blue')

#### Fourth endmember: Thermogenesis of biotic OM

dC_biothermo_OM = [dC_bio_CH4[0]+5, dC_bio_CH4[1] -1]
dC_biothermo_CH4 = get_CH4_minmax(dC_biothermo_OM, epsilon_thermo)

addlines(axs[0], 13.5, 14.5, dC_biothermo_OM, dC_biothermo_CH4 , 'tab:purple')





ax.axhline(dC_CO2_space, c='k', ls='dotted')

for _ax in axs:
    _ax.axvline(4, c='k')
    _ax.axvline(8, c='k')
    _ax.axvline(12, c='k')
    _ax.tick_params(axis='x', bottom=False, labelbottom=False)


ax.set_ylabel(r'$\delta^{13}$C $[\perthousand]$ when $\delta^{13}$C$_{x_{4}, \mathregular{CO2}}$ is 60 $\perthousand$')
# axs[1].set_ylabel(r'$\Delta\delta^{13}$C$_{x4, CO2/CH4}$ = $\delta^{13}$C$_{x4, CO2}$ - $\delta^{13}$C$_{x4, CH4}$ $[\perthousand]$')
plt.savefig(str(dC_CO2_space)+'a_.png', dpi=800)
plt.savefig(str(dC_CO2_space)+'a_.pdf')

plt.close()
