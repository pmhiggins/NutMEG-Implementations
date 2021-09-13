import sys, os, math
sys.path.append(os.path.dirname(__file__)+'../../NutMEG')
import ThesisSetup

import NutMEG as nm

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

R = nm.environment() # 1bar, 298 K

H2O = nm.reaction.reagent('H2O(l)', R, phase='l')
H = nm.reaction.reagent('H+', R, phase='aq')
OH = nm.reaction.reagent('OH-', R, phase='aq')

rxn = nm.reaction.reaction({H2O:1}, {H:1, OH:1}, R)

Ts = np.linspace(273,500, num=(500-273)*2)
pHs_01b = []
pHs_1b = []
pHs_100b = []
pHs_1000b = []

for T in Ts:
    rxn.env.T = T
    for P, pH in zip([1e5,1e7,1e8], [pHs_1b, pHs_100b, pHs_1000b]):
        rxn.env.P = P # NutMEG stores Pressures in Pa, 1 bar = 10^5 Pa
        rxn.rto_current_env() # updates std Gibbs and eq. constant
        pH.append(-math.log10(math.sqrt(math.exp(rxn.lnK))))

fig, ax = plt.subplots(figsize=(6,4), constrained_layout=True)
ax.plot(Ts, pHs_1b, c='rebeccapurple', label='Pressure = 1 bar')
ax.plot(Ts, pHs_100b, c='tab:orange', ls='dotted', label='Pressure = 100 bar')
ax.plot(Ts, pHs_1000b, c='tab:red', ls='dashed', label='Pressure = 1000 bar')

ax.set_ylabel('Neutral pH')
ax.set_xlabel('Temperature [K]')
ax.set_xlim(273,500)
ax.set_ylim(5.5,7.5)
ax.legend(loc='upper right', facecolor=(0.98,0.98, 0.98,1), framealpha=1.)
ax.grid(b=True, which='major', axis='both', color='#666666', linestyle='-', alpha=0.8)
plt.savefig('figs/pH_with_T.pdf', transparent=True)
