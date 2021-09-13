import sys, os, math
sys.path.append(os.path.dirname(__file__)+'../../NutMEG')
import ThesisSetup

import NutMEG as nm
from NutMEG.reactor.saved_systems.Enceladus import Enceladus
from NutMEG.culture.saved_organisms.TypicalOptimalMethanogen import TypicalOptimalMethanogen as TOM

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
from uncertainties import ufloat as uf

Ts = range(270, 400)

# make a reactor object
rtr = Enceladus('EncT')
# make an organism, for this example, use the preset methanogen.
org = TOM(rtr)
# org.dry_mass=29e-18 # set volume to 1 cubic micron.

TE = nm.apps_theory_estimates(org, rtr)
Tijhuis, TijhuisAerobe, TijhuisAnaerobe = [],[], []
Lever10pc, Lever2pc, Lever1AA = [],[],[] # in-built energy
Lever10pc_ces, Lever2pc_ces, Lever1AA_ces = [],[],[] # constant energy

for i, T in enumerate(Ts):
    TE.loc.change_T(T)
    TE.org.get_ESynth(AA=True) # update the synthesis energy
    td = TE.temperature_defenses(T)
    Tijhuis.append(uf(td['Tijhuis'], 0.29*td['Tijhuis']))
    TijhuisAerobe.append(uf(td['TijhuisAerobe'], 0.41*td['TijhuisAerobe']))
    TijhuisAnaerobe.append(uf(td['TijhuisAnaerobe'], 0.32*td['TijhuisAnaerobe']))
    Lever1AA.append(td['Lever1/250'])
    Lever10pc.append(td['Lever10pc'])
    Lever2pc.append(td['Lever2pc'])
    print(T, TE.org.E_synth/(1000*org.dry_mass), 1.7e-11 *(org.dry_mass/29e-18)/(1000*org.dry_mass))
    TE.org.E_synth = 1.7e-11 *(org.dry_mass/29e-18) # convert so dry mass matches this organism
    td = TE.temperature_defenses(T)
    Lever1AA_ces.append(td['Lever1/250'])
    Lever10pc_ces.append(td['Lever10pc'])
    Lever2pc_ces.append(td['Lever2pc'])


fig, ax = plt.subplots(figsize=(6,7), constrained_layout=True)

for Ti, c, l in zip([TijhuisAerobe, TijhuisAnaerobe], ['tab:olive','tab:blue'], ['aerobes', 'anaerobes']):
    Ti = np.array(Ti)
    ax.fill_between(Ts,
                    unp.nominal_values(Ti)-unp.std_devs(Ti),
                    unp.nominal_values(Ti)+unp.std_devs(Ti),
                    facecolor=c,
                    alpha=0.6,
                    linewidth=0,
                    label='Tijhuis et al. (1993) estimate for '+l)
ax.plot(Ts, Lever10pc, c='tab:green', label='Lever et al (2015), 10% racemization; our $E_\mathregular{syn}$')
ax.plot(Ts, Lever2pc, c='tab:green', linestyle='dashed', label='Lever et al (2015), 2% racemization; our $E_\mathregular{syn}$')
ax.plot(Ts, Lever1AA, c='tab:green', linestyle='dotted', label='Lever et al (2015), single amino acid racemization; our $E_\mathregular{syn}$')

ax.plot(Ts, Lever10pc_ces, c='tab:pink', label='Lever et al (2015), 10% racemization; their $E_\mathregular{syn}$')
ax.plot(Ts, Lever2pc_ces, c='tab:pink', linestyle='dashed', label='Lever et al (2015), 2% racemization; their $E_\mathregular{syn}$')
ax.plot(Ts, Lever1AA_ces, c='tab:pink', linestyle='dotted', label='Lever et al (2015), single amino acid racemization; their $E_\mathregular{syn}$')

ax.plot([278,278],
  [unp.nominal_values(TijhuisAnaerobe)[8]-unp.std_devs(TijhuisAnaerobe)[8],
    unp.nominal_values(TijhuisAerobe)[8]+unp.std_devs(TijhuisAerobe)[8]],
  c='k', solid_capstyle='round', label='limits of Tijhuis et al. (1993) data set')

ax.plot([330,330],
  [unp.nominal_values(TijhuisAnaerobe)[60]-unp.std_devs(TijhuisAnaerobe)[60],
    unp.nominal_values(TijhuisAerobe)[60]+unp.std_devs(TijhuisAerobe)[60]],
  c='k', solid_capstyle='round')

ax.set_title('Temperature-specific maintenance powers')
ax.set_ylabel('Maintenance Power [W cell$^{-1}$]')
ax.set_xlabel('Temperature [K]')
ax.set_yscale('log')
ax.set_xlim(Ts[0],Ts[-1])
ax.legend(bbox_to_anchor=(0., 1.08, 1., .102), loc=3,
       ncol=1, mode="expand", borderaxespad=0.)
# plt.show()
plt.savefig('Tmaintenance.pdf')
