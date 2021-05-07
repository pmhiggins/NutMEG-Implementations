import sys, os, ast
sys.path.append(os.path.dirname(__file__)+'../../NutMEG')
sys.path.append(os.path.dirname(__file__)+'../../NutMEG-Implementations/TOM')

import EnceladusPlotStyles as EPS

"""
This file was used to generate various plots for the manuscript.
To replicate, uncomment as needed.
"""

"""Methanogenesis energy plots (Figure 1)"""

## this generates figure 1 in the manuscript
EPS.make_MGEContourPlot(CO2origin='HTHeatingSalts', save='figs/energyplot_endmembersalt.pdf', show=False, pHax=True, pHbars=True, quotienttype='salty_endmember')

## this generates an energy plot with the salts fixed at nominal values
# EPS.make_MGEContourPlot(CO2origin='HTHeatingSalts', save=None, show=True, pHax=False, pHbars=False, quotienttype='salty_nominal')

## this generates an energy plot with the salts fixed at highest values
# EPS.make_MGEContourPlot(CO2origin='HTHeatingSalts', save='energyplottest_highsalt.pdf', show=False, pHax=False, pHbars=False, quotienttype='salty_high')

## this generates an energy plot with the salts fixed at lowest values
# EPS.make_MGEContourPlot(CO2origin='HTHeatingSalts', save='energyplottest_lowsalt.pdf', show=False, pHax=False, pHbars=False, quotienttype='salty_low')




"""Power Supply Nominal Plot (Figure xx)"""

# EPS.make_PSContourPlot(CO2origin='HTHeatingSalts', save='figs/PowerSupply_nominal.pdf', show=False, mesh=True)


"""Plus / minus power supply plots (Figure xx)"""

## endmember power supplies (manuscript figure xx)
# EPS.PSunc_plot(CO2origin='HTHeatingSalts', mesh=True, save='PowerSupply_comp_k.pdf')

## endmember plots with different nATP values
# EPS.PSunc_plot(CO2origin='HTHeatingSalts', mesh=True, save='suppfigs/PowerSupply_comp_k05.pdf', nATP=0.5)
# EPS.PSunc_plot(CO2origin='HTHeatingSalts', mesh=True, save='suppfigs/PowerSupply_comp_k15.pdf', nATP=1.5)
# EPS.PSunc_plot(CO2origin='HTHeatingSalts', mesh=True, save='suppfigs/PowerSupply_comp_k025.pdf', nATP=0.25)
# EPS.PSunc_plot(CO2origin='HTHeatingSalts', mesh=True, save='suppfigs/PowerSupply_comp_k200.pdf', nATP=2.0)





""" Habitability grids (Figure xx) """

# EPS.PShabitabilityPlot(nATP = 1.0)

## supplemental habitability plots
# EPS.PShabitabilityPlot(nATP = 2.0)
# EPS.PShabitabilityPlot(nATP = 0.25)
