import sys, os, math
sys.path.append(os.path.dirname(__file__)+'../../NutMEG')
import ThesisSetup

import NutMEG as nm
from NutMEG.reactor.saved_systems.Enceladus import Enceladus
from NutMEG.culture.saved_organisms.TypicalOptimalMethanogen import TypicalOptimalMethanogen as TOM

import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np

from copy import deepcopy
from itertools import chain

def SAget(v):
    return (4*math.pi*((v*3)**2))**(1/3)

def default_org():
    rtr = Enceladus('EncT', nominals=True, T=298)
    # make an organism, for this example, use the preset methanogen.
    org = TOM(rtr)
    org.maintenance.Tdef='None'
    org.maintenance.pHdef='FluxPerm'
    return org

def get_total_power(org, T, pHext, orgparams):
    org.memb_pot = orgparams.pop('memb_pot', -1e-5)
    org.PermH = orgparams.pop('PermH', 1e-10)
    org.PermOH = orgparams.pop('PermOH', 1e-10)
    org.base_volume = orgparams.pop('base_volume', 3.44e-18)
    org.surfacearea = SAget(org.base_volume)
    org.pH_interior = orgparams.pop('pH_interior', 7.)
    org.maintenance.Tdef = orgparams.pop('Tdef', 'TijhuisAnaerobe')
    org.locale.env.T = T
    org.locale.update_pH(pHext, _from='pH')
    pHad = nm.pHadaptations(org)
    org.maintenance.get_P_pH()
    org.maintenance.get_P_T()
    pow_pH = org.maintenance.net_dict['pH']
    pow_T = org.maintenance.net_dict['T']
    # print(org.PermH, pow_pH)
    pow_tot = pow_pH + pow_T
    return pow_pH, pow_T, pow_tot


def get_power_grid(orgparams, Tspace=None, pHspace=None):
    if Tspace==None:
        Tspace = np.linspace(273,400, num=100)
    if pHspace==None:
        pHspace = np.linspace(0,14, num=100)

    org = default_org()

    powers = np.ndarray((len(Tspace), len(pHspace)))

    yindex = 0
    for T in Tspace:
        xindex = 0
        for pH in pHspace:
            powers[yindex][xindex] = get_total_power(org, T, pH, deepcopy(orgparams))[2]
            xindex += 1
        yindex += 1
    pHmesh, Tmesh = np.meshgrid(pHspace, Tspace)
    return pHmesh, Tmesh, powers

fig, axs2 = plt.subplots(nrows=3, ncols=3, figsize=(8,8))
axs = axs2.flatten()

T6 = get_power_grid({'Tdef':'TijhuisAnaerobe', 'PermH':1e-6, 'PermOH':1e-6})
T9 = get_power_grid({'Tdef':'TijhuisAnaerobe', 'PermH':1e-9, 'PermOH':1e-9})
T12 = get_power_grid({'Tdef':'TijhuisAnaerobe', 'PermH':1e-12, 'PermOH':1e-12})
L10_6 = get_power_grid({'Tdef':'Lever10pc', 'PermH':1e-6, 'PermOH':1e-6})
L10_9 = get_power_grid({'Tdef':'Lever10pc', 'PermH':1e-9, 'PermOH':1e-9})
L10_12 = get_power_grid({'Tdef':'Lever10pc', 'PermH':1e-12, 'PermOH':1e-12})
LAA_6 = get_power_grid({'Tdef':'Lever1/250', 'PermH':1e-6, 'PermOH':1e-6})
LAA_9 = get_power_grid({'Tdef':'Lever1/250', 'PermH':1e-9, 'PermOH':1e-9})
LAA_12 = get_power_grid({'Tdef':'Lever1/250', 'PermH':1e-12, 'PermOH':1e-12})

for xyz, ax in zip([T6,T9,T12,LAA_6,LAA_9,LAA_12,L10_6,L10_9,L10_12], axs):
    contf = ax.contourf(xyz[0], xyz[1], np.log10(xyz[2]), cmap='coolwarm', levels=np.arange(-21, -8, 1), vmin=-21, vmax=-7)
# contf = axs[0].contourf(LAA_12[0], LAA_12[1], np.log10(LAA_12[2]), cmap='coolwarm')

axs[0].annotate('Perm = 10$^{-6}$', xy=(0.5, 1.1), xytext=(0.5, 1.2), xycoords='axes fraction',
        fontsize=12, ha='center', va='center', rotation='horizontal',
        arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=2))

axs[1].annotate('Perm = 10$^{-9}$', xy=(0.5, 1.1), xytext=(0.5, 1.2), xycoords='axes fraction',
        fontsize=12, ha='center', va='center', rotation='horizontal',
        arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=2))

axs[2].annotate('Perm = 10$^{-12}$', xy=(0.5, 1.1), xytext=(0.5, 1.2), xycoords='axes fraction',
        fontsize=12, ha='center', va='center', rotation='horizontal',
        arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=2))

axs[2].annotate('Tijhuis (1993)', xy=(1.1, 0.5), xytext=(1.2, 0.5), xycoords='axes fraction',
        fontsize=12, ha='left', va='center', rotation='vertical',
        arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=2))

axs[5].annotate('Lever et al. (2015) 1 AA', xy=(1.1, 0.5), xytext=(1.2, 0.5), xycoords='axes fraction',
        fontsize=12, ha='left', va='center', rotation='vertical',
        arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=2))

axs[8].annotate('Lever et al. (2015) 10%', xy=(1.1, 0.5), xytext=(1.2, 0.5), xycoords='axes fraction',
        fontsize=12, ha='left', va='center', rotation='vertical',
        arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=2))

for ax in axs2[-1]:
    ax.set_xlabel('pH')
for ax in [axs[0], axs[3], axs[6]]:
    ax.set_ylabel('Temperature [K]')

for a in chain(axs2[0], axs2[1]):
    a.get_xaxis().set_ticks([])

for i in axs2:
    for a in i[1:]:
        a.get_yaxis().set_ticks([])

cbaxes = fig.add_axes([0.075, 0.065, 0.92-0.075, 0.02])

fig.colorbar(contf, cax=cbaxes, label='Total power requirement [W cell$^{-1}$]', orientation='horizontal', pad=0.5, aspect=7, extend='both')
fig.subplots_adjust(left=0.075, bottom=0.15, right=0.92, top=0.933, wspace=0.05, hspace=0.05)
plt.savefig('pHTgrid.pdf')
