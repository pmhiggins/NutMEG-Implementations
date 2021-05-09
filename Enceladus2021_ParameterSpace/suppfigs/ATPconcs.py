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
mpl.rcParams['lines.linewidth'] = 3

mpl.rc('axes', unicode_minus=False)
mpl.rcParams['hatch.linewidth'] = 2.0
# print(mpl.rcParams.keys()) # to see what can be changed
mpl.rcParams['font.size'] = 14


fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5,5))

Ts = range(273,473)
GATP, worseGATP, minGATP, betterGATP, maxGATP = [],[],[], [],[]

for T in Ts:
    E = Enc('Enceladus', T=T, depth=0)
    TOMobj = TOM(E, paramchange={'Basal':1e-15, 'Tdef':'None'})

    GATP.append(TOMobj.respiration.G_P/1000)

    # TOMobj.respiration.build_ATP_reaction([0.0002, 0.004, 0.0056, 7.])#[0.00001,0.0004,0.05,7.])

    worseGATP.append(TOMobj.respiration.G_P*0.5/1000)
    minGATP.append(TOMobj.respiration.G_P*0.25/1000)


    # TOMobj.respiration.build_ATP_reaction([0.002,0.022,0.0007,7.])#[0.001,0.04,0.0005,7.])
    betterGATP.append(TOMobj.respiration.G_P*1.5/1000)
    maxGATP.append(TOMobj.respiration.G_P*2.0/1000)



axs.set_title('ATP Production')

axs.plot(Ts, minGATP, color='m', linestyle='dotted', label=r'$n_\mathregular{ATP} = 0.25$')
axs.plot(Ts, worseGATP, color='m', linestyle='dashed', label=r'$n_\mathregular{ATP} = 0.5$')

axs.plot(Ts, GATP, color='b', label=r'$n_\mathregular{ATP} = 1.0$')
axs.plot(Ts, betterGATP, color='c', linestyle='dashed', label=r'$n_\mathregular{ATP} = 1.5$')
axs.plot(Ts, maxGATP, color='c', linestyle='dotted', label=r'$n_\mathregular{ATP} = 2.0$')
# axs.plot(Ts, worseGATP, color='b', linestyle='dashed', label='0.0002, 0.004, 0.0056')
# axs.plot(Ts, betterGATP, color='b', linestyle='dotted', label='0.002,0.022,0.0007')

# axs[1].plot(Ts, GATP80, color = 'b', linestyle='dashed', label='80 bar')

for ax in [axs]:
    ax.set_xlabel('Temperature [K]')
    ax.set_xlim(273,473)
    ax.set_ylabel('Gibbs free energy [kJ/mol]')
    ax.legend(bbox_to_anchor=(0., 1.12, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
plt.tight_layout()
plt.savefig('ATPconcs_minmax.pdf')
plt.show()
