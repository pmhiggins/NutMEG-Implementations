import sys, os, ast
sys.path.append(os.path.dirname(__file__)+'../../../NutMEG')
sys.path.append(os.path.dirname(__file__)+'../../../NutMEG-Implementations/TOM')
sys.path.append(os.path.dirname(__file__)+'../')


import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import EnceladusPlotStyles as EPS


##### Needs titles, labels, ticks on colorbar, legend? sorted

mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13
mpl.rcParams['font.size'] = 13

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10,10))

axs[0][0], contf = EPS.MethanogenesisEnergyContourPlot(axs[0][0],
  pHrange=np.linspace(7,12, num=11), CO2origin='HTHeatingSalts', quotienttype='salty_nominal')

axs[1][0], contf2 = EPS.MethanogenesisEnergyContourPlot(axs[1][0],
  pHrange=np.linspace(7,12, num=11), CO2origin='HTHeatingSalts', quotienttype='salty_low')

axs[1][1], contf3 = EPS.MethanogenesisEnergyContourPlot(axs[1][1],
  pHrange=np.linspace(7,12, num=11), CO2origin='HTHeatingSalts', quotienttype='salty_high')

fig.delaxes(axs[0][1])

# fig.subplots_adjust(bottom=0.08, left=0.08, right=0.92, top=0.92, wspace=0.05, hspace=0.08)

plt.subplots_adjust(wspace=0.2, hspace=0.33, left=0.08, right=0.98, bottom=0.05, top=0.95)

axs[0][0].set_title('Nominal salt level')
axs[1][0].set_title('Low salt level')
axs[1][1].set_title('High salt level')


cbaxes = fig.add_axes([0.58, 0.6, 0.3975, 0.05])
fig.colorbar(contf, cax=cbaxes, ticks=[-160,-120,-80,-40,0,40,80], label='Free Energy of methanogenesis [kJ/mol]', orientation='horizontal', extend='both')


axs[0][0].legend(bbox_to_anchor=(1.2, 0.8, 1.0, .102), loc=3,
       ncol=1, mode="expand", borderaxespad=0.)

for ax in axs:
    for a in ax:
        a.set_xlim(7,12)
        a.set_ylim(273.15,473)
        a.set_ylabel('Temperature [K]')
        a.set_xlabel('Bulk ocean pH')

plt.savefig('freeenergysalts.pdf')

# EPS.make_MGEContourPlot(CO2origin='HTHeatingSalts', save='energyplottest3.pdf', show=False, pHax=False, pHbars=False, quotienttype='salty_nominal')
#
# EPS.make_MGEContourPlot(CO2origin='HTHeatingSalts', save='energyplottest_highsalt.pdf', show=False, pHax=False, pHbars=False, quotienttype='salty_high')
#
# EPS.make_MGEContourPlot(CO2origin='HTHeatingSalts', save='energyplottest_lowsalt.pdf', show=False, pHax=False, pHbars=False, quotienttype='salty_low')
