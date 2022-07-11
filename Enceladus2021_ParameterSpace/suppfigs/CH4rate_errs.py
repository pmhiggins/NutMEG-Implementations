import sys, os, ast, math
sys.path.append(os.path.dirname(__file__)+'../../../NutMEG')

import NutMEG
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from NutMEG.reactor.saved_systems.Enceladus import Enceladus as Enc
from NutMEG.reactor.saved_systems.VenusDrop import VenusDrop
from NutMEG.culture.saved_organisms.TypicalOptimalMethanogen import TypicalOptimalMethanogen as TOM

# change matplotlib rcparams updating plot fonts etc.
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cmr10'
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['lines.markersize'] = 8


def CH4get(growthrate, Pres=182000, CH4conc=3e-8):
    """Use reverse Powell method to compute CH4 production rate from growth rate"""
    # CH4conc = 10**((Pres*2.5e-5)-12.5) #Toggle me to include a range of CH4
    # print(CH4conc)
    return CH4conc*growthrate*(math.exp(growthrate) - 1)

### For loading data from teh empirical methanogens

file = '../../../NutMEG-Implementations/TOM/data/methanogens.csv'
df = pd.read_csv(file, header=0)
empTs, CH4s, CO2s, H2s, Pressures, GRs = [],[],[],[],[],[]
for index, row in df.iterrows():
    try:
        T = 273+((float(row['Min. optimal growth temp.'])+float(row['Max. optimal growth temp.'])))/2
        empTs.append(T)
        P = float(row['Pressure'])*1000
        Pressures.append(P)
        GRs.append(row['Growth rate']/3600)
        CH4s.append(CH4get(row['Growth rate']/3600, Pres=P))
        CO2s.append(VenusDrop.getgasconc('CO2(aq)', 0.2*P, T, P_bar=P, S=0))
        H2s.append(VenusDrop.getgasconc('H2(aq)', 0.8*P, T, P_bar=P, S=0))
    except Exception as e:
        print(str(e))
        continue

# CH4pol = np.polyfit(Ts, np.log(CH4s), 1)
# GRspol = np.polyfit(Ts, np.log(GRs), 1)

def concupdate(R):
    """Change the concentrations of H2 and CO2 in reactor R to replicate TOM conditions at the reactor's temperature and pressure."""
    R.composition['H2(aq)'].activity = VenusDrop.getgasconc('H2(aq)', 0.8*R.env.P,
      R.env.T, P_bar=R.env.P)
    R.composition['CO2(aq)'].activity = VenusDrop.getgasconc('CO2(aq)', 0.2*R.env.P,
      R.env.T, P_bar=R.env.P)
    return R



Ts=range(273,373)
ks, kTs, met = [],[], []
kTs30, kTs10 = [], []
CO2_30, CO2_10, CO2_18, CO2_500 = [],[],[],[]
H2_30, H2_10, H2_18, H2_500 = [],[],[],[]

for T in Ts:
    E = Enc('Enceladus', T=T, nominals=True)
    E.env.P=182000
    E = concupdate(E)
    TOMobj = TOM(E, paramchange={'Basal':1e-15, 'Tdef':'None'}, fromdata=False)

    ks.append( TOMobj.respiration.net_pathway.rate_constant_RTP)
    met.append(TOMobj.max_metabolic_rate)
    kTs.append(TOMobj.respiration.net_pathway.rate_constant_env)
    for r, mr in TOMobj.respiration.net_pathway.reactants.items():
        if r.name=='CO2(aq)':
            CO2_18.append(r.activity)
        elif r.name=='H2(aq)':
            H2_18.append(r.activity)

    E.env.P = 300000
    E = concupdate(E)
    TOM300 = TOM(E, paramchange={'Basal':1e-15, 'Tdef':'None', 'TOMpressure':300000}, fromdata=True)
    kTs30.append(TOM300.respiration.net_pathway.rate_constant_env)
    for r, mr in TOM300.respiration.net_pathway.reactants.items():
        if r.name=='CO2(aq)':
            CO2_30.append(r.activity)
        elif r.name=='H2(aq)':
            H2_30.append(r.activity)

    E.env.P = 100000
    E = concupdate(E)
    TOM100 = TOM(E, paramchange={'Basal':1e-15, 'Tdef':'None', 'TOMpressure':100000}, fromdata=True)
    kTs10.append(TOM100.respiration.net_pathway.rate_constant_env)
    for r, mr in TOM100.respiration.net_pathway.reactants.items():
        if r.name=='CO2(aq)':
            CO2_10.append(r.activity)
        elif r.name=='H2(aq)':
            H2_10.append(r.activity)

    E.env.P = 5000000
    E = concupdate(E)
    TOM5000 = TOM(E, paramchange={'Basal':1e-15, 'Tdef':'None', 'TOMpressure':5000000}, fromdata=True)
    for r, mr in TOM100.respiration.net_pathway.reactants.items():
        if r.name=='CO2(aq)':
            CO2_500.append(r.activity)
        elif r.name=='H2(aq)':
            H2_500.append(r.activity)



"""  For plotting variation of k with pressure (substrate conc)
empK = []
for i in range(len(empTs)):
    empK.append(CH4s[i] / (CO2s[i]*(H2s[i]**4)))

# fig=plt.figure()
# ax = fig.add_subplot(111)#, projection='3d')

plt.scatter(empTs, empK, c=np.array(Pressures)*0.001)
plt.plot(Ts, kTs, c=mpl.cm.get_cmap()(6.5/12), lw=3, label='TOM match: 182 kPa')
plt.plot(Ts, kTs30, c=mpl.cm.get_cmap()(0.98), lw=3, label='300 kPa')
plt.plot(Ts, kTs10, c=mpl.cm.get_cmap()(1.5/12), lw=3, label='100 kPa')

plt.yscale('log')
plt.xlabel('Temperature (K)', fontsize=14)
plt.ylabel('Rate Const (/s)', fontsize=14)
# plt.ylim(-17,3)
cb = plt.colorbar()
cb.ax.set_ylabel('Pressure [kPa]',fontsize=14)
plt.legend(loc='best')
plt.tight_layout()
plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.8)
plt.savefig('kT_oldCH4_newTOMPress.pdf')
# plt.show()
"""

""" GROWTH RATES """


""" For plotting metabolic and growth rate changing with pressure (substrate conc) """

# first we need the polyfit
# GRpoly = np.polyfit(empTs, np.log(GRs), 1)
# GRlin = np.exp((GRpoly[0]*Ts) + GRpoly[1])

CH4poly = np.polyfit(empTs, np.log(CH4s), 1)
CH4lin = np.exp((CH4poly[0]*Ts) + CH4poly[1])
CH4_p1 = 10**(np.log10(CH4lin)+1)
CH4_m1 = 10**(np.log10(CH4lin)-1)


fig, axs = plt.subplots(nrows=1, figsize=(7,5))
cm = axs.scatter(empTs, CH4s, c=np.array([Pressures])*0.001)
axs.plot(Ts, CH4lin, lw=3, c='k', label='nominal k')
axs.plot(Ts, CH4_p1, lw=3, c='k', linestyle='dashed', label='bounds of k')
axs.plot(Ts, CH4_m1, lw=3, c='k', linestyle='dashed')
axs.legend(fontsize=14)


# cm = axs[1].scatter(empTs, GRs, c=np.array([Pressures])*0.001)
# axs[1].plot(Ts, GRlin, lw=3, c='k')

axs.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.8)
axs.set_yscale('log')

axs.set_ylabel('$\mathregular{CH}_4$ production rate $[\mathregular{M\ (cell s)}^{-1}]$', fontsize=14)
# axs[1].set_ylabel('Growth rate $[\mathregular{s}^{-1}]$', fontsize=14)
axs.set_xlabel('Temperature [K]', fontsize=14)
# plt.ylim(-17,3)


# fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.85, 0.097, 0.03, 0.979-0.097])
cb = fig.colorbar(cm, label='Pressure [Pa]', orientation='vertical')
cb.ax.tick_params(labelsize=14)
cb.ax.set_ylabel('Pressure [kPa]',fontsize=14)

plt.tight_layout()
# plt.legend(loc='best')
plt.savefig('CH4rate.pdf')
# plt.show()






""" Concentrations at different pressures """
"""
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)
ax.plot(Ts, CO2_18, label='180 kPa', c='k', linewidth=4)
ax.plot(Ts, H2_18, c='r', linewidth=4)
ax.plot(Ts, CO2_30, label='300 kPa', c='k', ls='--', linewidth=4)
ax.plot(Ts, H2_30, c='r', ls='--', linewidth=4)
ax.plot(Ts, CO2_10, label='100 kPa', c='k', ls='dotted', linewidth=4)
ax.plot(Ts, H2_10, c='r', ls='dotted', linewidth=4)
ax.plot(Ts, CO2_500, label='500 kPa', c='k', ls='dashdot', linewidth=4)
ax.plot(Ts, H2_500, c='r', ls='dashdot', linewidth=4)
ax.set_ylabel('Concentration [M]', fontsize=14)
ax.set_xlabel('Seawater temperature [K]', fontsize=14)
ax.set_yscale('log')
ax.legend(loc='best')
plt.tight_layout()
plt.savefig('concs_withPbar.pdf') # change VenusDrop to change Pbar
"""
