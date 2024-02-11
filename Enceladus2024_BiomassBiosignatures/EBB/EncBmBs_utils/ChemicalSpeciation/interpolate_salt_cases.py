import numpy as np
import math
import matplotlib.pyplot as plt
import reaktoro as rkt
import pandas as pd
import sys, os


Carbs = np.logspace(math.log10(0.01),math.log10(0.1), num=21) #DIC
Cl = np.logspace(math.log10(0.05),math.log10(0.2), num=21)
Cl2 = np.logspace(math.log10(0.05),math.log10(0.4), num=31)


pHs = np.linspace(7,12,num=51)
dom_ion = ['HCO3-', 'HCO3-','HCO3-','HCO3-','HCO3-','HCO3-','HCO3-','HCO3-','CO3-2','CO3-2','CO3-2']
def dom_ion(pH_val):
    if pH_val <11:
        return 'HCO3-'
    else:
        return 'CO3-2'


DIC_tots = np.zeros((len(Carbs), len(pHs)))
CO2_tots = np.zeros((len(Carbs),len(pHs)))


for i in range(len(Carbs)):
    inc = 0.
    _this_DIC_tots = []
    for pH_i, pH in enumerate(pHs):
        CO2 = 0.
        for _ in range(10):

            db = rkt.PhreeqcDatabase('phreeqc.dat')
            solution = rkt.AqueousPhase(rkt.speciate("H O Na Cl C"), rkt.exclude("organic"))
            solution.setActivityModel(rkt.ActivityModelPitzerHMW())


            system = rkt.ChemicalSystem(db, solution)

            specs = rkt.EquilibriumSpecs(system)
            specs.temperature()
            specs.pressure()
            specs.pH()
            specs.charge()
            specs.openTo('Na+')

            solver = rkt.EquilibriumSolver(specs)

            state = rkt.ChemicalState(system)
            state.temperature(273.15, "kelvin")
            state.pressure(1, "bar")

            state.set('H2O', 1., "kg")
            # state.set('Na+', mNa, "mol")
            state.set('Cl-', Cl[i], "mol")
            # state.set('K+', mK, "mol")
            state.set(dom_ion(pH_i), Carbs[i]+inc, "mol")
            state.set('H+', 10**-pH, "mol")

            conditionsT = rkt.EquilibriumConditions(specs)
            conditionsT.temperature(273.15, "kelvin")
            conditionsT.pressure(1, "bar")
            conditionsT.pH(pH)
            conditionsT.charge(0.)

            solver.solve(state, conditionsT)

            p = rkt.ChemicalProps(state)
            # print(p)

            CO2 = float(p.speciesConcentration('CO2'))
            HCO3 = float(p.speciesConcentration('HCO3-'))
            CO3 = float(p.speciesConcentration('CO3-2'))
            NaHCO3 = float(p.speciesConcentration('NaHCO3'))
            NaCO3 = float(p.speciesConcentration('NaCO3-'))

            inc += Carbs[i] - HCO3 - CO3 - NaHCO3 - NaCO3
            # inc += inc_by
        DIC_tots[i][pH_i] = Carbs[i]+inc
        CO2_tots[i][pH_i] = CO2


np.save('DIC_Cl_correlation/DIC_tots.npy', DIC_tots)
np.save('DIC_Cl_correlation/CO2_tots.npy', CO2_tots)


DIC_tots = np.load('DIC_Cl_correlation/DIC_tots.npy')
CO2_tots = np.load('DIC_Cl_correlation/CO2_tots.npy')



# THIS IS FOR EXTENDING THE SALT CASES UP TO 0.4 M Cl

DIC_tots_HS = np.zeros((len(Cl2), len(pHs)))
CO2_tots_HS = np.zeros((len(Cl2), len(pHs)))

for i in range(len(pHs)):
    polyDIC = np.polyfit(np.log10(Cl), np.log10(DIC_tots.T[i]), deg=1)
    polyCO2 = np.polyfit(np.log10(Cl), np.log10(CO2_tots.T[i]), deg=1)
    for j in range(len(Cl2)):
        if j<21:
            # there is a ~ 1% error in the poly fit for DIC, so keep the
            # original values for the fost 21 salinities. This will also keep it
            # aligned with any previous results.
            DIC_tots_HS[j][i] = DIC_tots[j][i]
            CO2_tots_HS[j][i] = CO2_tots[j][i]
        else:
            DIC_tots_HS[j][i] = 10**(polyDIC[1] + (polyDIC[0]*math.log10(Cl2[j])))
            CO2_tots_HS[j][i] = 10**(polyCO2[1] + (polyCO2[0]*math.log10(Cl2[j])))

np.save('DIC_tots_HS.npy', DIC_tots_HS)
np.save('CO2_tots_HS.npy', CO2_tots_HS)



for DIC, pH in zip(DIC_tots.T, pHs):
    poly = np.polyfit(np.log10(Cl), np.log10(DIC), deg=1)
    plt.plot(Cl, DIC, label='pH '+str(pH))
    plt.plot(Cl2, 10**(poly[1] + (poly[0]*np.log10(Cl2))), ls='dashed')
    # plt.plot(Cl2, DIC_tots_HS, ls='dotted')

    plt.xlabel('[Cl]')
    plt.ylabel('[DIC]')
    # plt.yscale('log')
    # plt.xscale('log')

    plt.legend()
plt.show()
