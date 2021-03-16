import sys, os, ast, math
sys.path.append(os.path.dirname(__file__)+'../../../NutMEG')
import matplotlib.pyplot as plt
import matplotlib as mpl

import NutMEG
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


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(5,8))

Ts = range(273,473) #range(273,373)
GM10, GATP10 = [],[]
GM80, GATP80 = [],[]

for T in range(273,473):
    E = Enc('Enceladus', T=T, depth=0)
    TOMobj = TOM(E, paramchange={'Basal':1e-15, 'Tdef':'None'})

    GM10.append(TOMobj.respiration.net_pathway.std_molar_gibbs)
    GATP10.append(TOMobj.respiration.ATP_production.std_molar_gibbs)

    E = Enc('Enceladus', T=T, depth=10)
    TOMobj = TOM(E, paramchange={'Basal':1e-15, 'Tdef':'None'})

    GM80.append(TOMobj.respiration.net_pathway.std_molar_gibbs)
    GATP80.append(TOMobj.respiration.ATP_production.std_molar_gibbs)

axs[0].set_title('Methanogenesis')
axs[0].plot(Ts, GM10, color='b', label='10 bar')
axs[0].plot(Ts, GM80, color = 'b', linestyle='dashed', label='80 bar')

axs[1].set_title('ATP Production')
axs[1].plot(Ts, GATP10, color='b', label='10 bar')
axs[1].plot(Ts, GATP80, color = 'b', linestyle='dashed', label='80 bar')

for ax in axs:
    ax.set_xlabel('Temperature [K]')
    ax.set_xlim(273,473)
    ax.set_ylabel('Standard free energy [J/mol]')
    ax.legend()
plt.tight_layout()
plt.savefig('stds_Meth_ATP.pdf')
plt.show()
