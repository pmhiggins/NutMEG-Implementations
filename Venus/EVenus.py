"""
Code for calaculations performed in Cockell et al. (2020).
Estmating dissolved concentrations and energetic availability
of known metabolisms in Venusian cloud droplets. Edit the
docstings from line ~430 and downwards to recreate
calculations and plots from the study.

Tested and working on NutMEG master branch and all active
branches as of April 2020.

Conducted at the University of Edinburgh, 2020.

authors : Pete Higgins & Andrew Johnstone
contact : p.m.higgins@ed.ac.uk

"""

import sys, os, math, sqlite3, csv
# the below is needed if you have NutMEG in an adjacent directory
sys.path.append(os.path.dirname(__file__)+'../../')
from copy import deepcopy
import numpy as np
import PIL # only required if you want to save images as .tif
import matplotlib.pyplot as plt

import NutMEG as es
from NutMEG.reactor.saved_systems.VenusDrop import VenusDrop
import NutMEG.reaction as reaction


def setup_sulfatereduction(R, k_RTP=1):
    """ Set up a sulate reduction reaction to pass to the sulfate reducer
    as its metabolism. Use the reactor R."""

    # set up reagents. No need for composition as that is saved in the reactor.
    H2aq = reaction.reagent('H2(aq)', R.env, phase='aq')
    SO4 = reaction.reagent('SO4--', R.env, phase='aq', charge=-2)
    H = reaction.reagent('H+', R.env, phase='aq', charge=1)
    HS = reaction.reagent('HS-', R.env, phase='aq', charge=-1)
    H2O = reaction.reagent('H2O(l)', R.env, phase='l')

    # put reagents into a reaction. Dictionary values are the molar ratios.
    thermalSR = reaction.reaction({H2aq:4, SO4:1, H:1}, {HS:1, H2O:4},
      R.env)

    # not important for this code, but a rate constant is required
    # to set up an organism in NutMEG.
    thermalSR.rate_constant_RTP = k_RTP

    return thermalSR

def setup_methanogenesis(R, k_RTP=1):
    """ Set up a methanogenesis reaction to pass to the methanogen
    as its metabolism. Use the reactor R."""

    # set up reagents. No need for composition as that is saved in the reactor.
    H2aq = reaction.reagent('H2(aq)', R.env, phase='aq')
    H2O = reaction.reagent('H2O(l)', R.env, phase='l')
    CO2 = reaction.reagent('CO2(aq)', R.env, phase='aq')
    CH4 = reaction.reagent('CH4(g)', R.env, phase='aq')

    # put reagents into a reaction. Dictionary values are the molar ratios.
    thermalSR = reaction.reaction({CO2:1, H2aq:4}, {CH4:1, H2O:2},
      R.env)

    # not important for this code, but a rate constant is required
    # to set up an organism in NutMEG.
    thermalSR.rate_constant_RTP = k_RTP


    return thermalSR

def setup_h2oxidation(R, k_RTP=1):
    """ Set up a H2 oxidation reaction to pass to the H2 oxidiser
    as its metabolism. Use the reactor R."""

    # set up reagents. No need for composition as that is saved in the reactor.
    H2aq = reaction.reagent('H2(aq)', R.env, phase='aq')
    O2 = reaction.reagent('O2(aq)', R.env, phase='aq')
    H2O = reaction.reagent('H2O(l)', R.env, phase='l')

    # put reagents into a reaction. Dictionary values are the molar ratios.
    thermalSR = reaction.reaction({H2aq:2, O2:1,}, {H2O:2}, R.env)

    # not important for this code, but a rate constant is required
    # to set up an organism in NutMEG.
    thermalSR.rate_constant_RTP = k_RTP

    return thermalSR

def update_all_comp_ppm(R):
    """Set all the relevant reagents to be their `standard` ppm
    values at the current temperature etc.

    R : VenusDrop
        VenusDrop-like reactor to update
    """
    R.update_reagent('H2(aq)', 25)
    R.update_reagent('CO2(aq)', 940000)
    R.update_reagent('CH4(g)', 5)
    R.update_reagent('O2(aq)', 43.6)


def Venus_atm_energy(org_name, rxn_func, variable,temps=(278, 314, 350)):
    """ Calculate dissolved concentration of constituents and molar
    free energy of metabolism for organism org_name, at temperatures
    temps. Return three lists: the atmospheric ppms used, the three concentrations used at the temperatures

    org_name : str
        Name of the organism
    rxn_func : func
        function which initialises metabolism
    variable : str
        Name of reagent to vary
    temps : tuple
        Three temperatures to use (representative of heights
        in the atmosphere)
    """
    # setup specific reactor modelled on venusian cloud drops
    V=VenusDrop(H2ppm=25, workoutID=False)
    update_all_comp_ppm(V)

    # initialise the organism .Passing V means we
    # put it in our VenusDrop
    SR = es.base_organism(org_name, V,
      rxn_func(V), workoutID=False)

    # temperature lower bound
    mol_L = [] # conc variable
    Q_L = [] # reaction quotient
    K_L = [] # equilibrium constant
    stdG_L = [] # standard gibbs per mol
    G_L = [] #Available energy per mol

    # temperature middle value
    mol_M=[]
    Q_M =[]
    K_M = []
    stdG_M = []
    G_M =[]

    # temperature higher bound
    mol_H=[]
    Q_H =[]
    K_H = []
    stdG_H = []
    G_H =[]

    # iterate over ppm values dependent on varaible selected
    if variable=='CO2(aq)':
        ppmlst=np.linspace(500000, 1050000, num=101)
    elif variable=='H2(aq)':
        ppmlst=np.linspace(1, 50, num=101)
    elif variable=='O2(aq)':
        ppmlst=np.linspace(1, 100, num=101)
    elif variable=='CH4(g)':
        ppmlst=np.linspace(0.1, 30, num=301)

    for ppm in ppmlst:


        V.env.T = temps[0]
        # find conc at this ppm and update CO2 concentration
        # also update other constituent concs, as beta changes with T
        update_all_comp_ppm(V)
        V.update_reagent(variable, ppm)
        V.unify_reaction(SR.respiration.net_pathway, overwrite=False)

        # update the quotient, std gibbs, and molar gibbs in
        # the current environment.
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        # extract the parameters we're interested in from the organism
        mol_L.append(V.composition[variable].conc)
        Q_L.append(SR.respiration.net_pathway.quotient)
        K_L.append(math.exp(SR.respiration.net_pathway.lnK))
        stdG_L.append(SR.respiration.net_pathway.std_molar_gibbs)
        G_L.append(SR.respiration.net_pathway.molar_gibbs)


        ##### now repeat for the other two temperatures
        V.env.T = temps[1]
        update_all_comp_ppm(V)
        V.update_reagent(variable, ppm)
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        mol_M.append(V.composition[variable].conc)
        Q_M.append(SR.respiration.net_pathway.quotient)
        K_M.append(math.exp(SR.respiration.net_pathway.lnK))
        stdG_M.append(SR.respiration.net_pathway.std_molar_gibbs)
        G_M.append(SR.respiration.net_pathway.molar_gibbs)


        V.env.T = temps[2]
        update_all_comp_ppm(V)
        V.update_reagent(variable, ppm)
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        mol_H.append(V.composition[variable].conc)
        Q_H.append(SR.respiration.net_pathway.quotient)
        K_H.append(math.exp(SR.respiration.net_pathway.lnK))
        stdG_H.append(SR.respiration.net_pathway.std_molar_gibbs)
        G_H.append(SR.respiration.net_pathway.molar_gibbs)


    # return ppm variation, concs, and molar gibbs.
    # edit to return K, Q, or standard gibbs if desired.
    return ppmlst, [mol_L, mol_M, mol_H], [G_L, G_M, G_H]

def Venus_overall_methanogenesis():
    """ maximum and minimum molar free energies from methanogenesis """
    # setup specific reactor modelled on venusian cloud drops
    V=VenusDrop(H2ppm=25, workoutID=False)
    update_all_comp_ppm(V)

    SR = es.base_organism('methanogen', V,
      setup_methanogenesis(V), workoutID=False)

    # set up lists for plotting
    minG=[]
    midG=[]
    maxG=[]

    Tlst=np.linspace(278,350, num=100)

    for T in Tlst:
         ### 25C
        V.env.T = T

        #for minG, maximise reactants, minimise products
        V.update_reagent('H2(aq)', 35)
        V.update_reagent('CO2(aq)', 1000000)
        V.update_reagent('CH4(g)', 0.1)
        V.unify_reaction(SR.respiration.net_pathway, overwrite=False)

        # this method updates the quotient, std gibbs, and molar gibbs in
        # the current environment.
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        minG.append(SR.respiration.net_pathway.molar_gibbs)

        #for midG, use recorded values
        V.update_reagent('H2(aq)', 25)
        V.update_reagent('CO2(aq)', 940000)
        V.update_reagent('CH4(g)', 5)
        V.unify_reaction(SR.respiration.net_pathway, overwrite=False)

        # this method updates the quotient, std gibbs, and molar gibbs in
        # the current environment.
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        midG.append(SR.respiration.net_pathway.molar_gibbs)

        #for maxG, maximise products, minimise reactants
        V.update_reagent('H2(aq)', 15)
        V.update_reagent('CO2(aq)', 750000)
        V.update_reagent('CH4(g)', 10)
        V.unify_reaction(SR.respiration.net_pathway, overwrite=False)

        # this method updates the quotient, std gibbs, and molar gibbs in
        # the current environment.
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        maxG.append(SR.respiration.net_pathway.molar_gibbs)

    return Tlst, minG, midG, maxG

def Venus_overall_sulfatereduction():
    """ maximum and minimum molar free energies from sulfate reduction """
    # setup specific reactor modelled on venusian cloud drops
    V=VenusDrop(H2ppm=25, workoutID=False)
    update_all_comp_ppm(V)

    SR = es.base_organism('SulfateReducer', V,
      setup_sulfatereduction(V), workoutID=False)

    # set up lists for plotting
    minG=[]
    midG=[]
    maxG=[]

    Tlst=np.linspace(278,350, num=100)

    for T in Tlst:
         ### 25C
        V.env.T = T
        update_all_comp_ppm(V)

        #for minG, maximise reactants, minimise products
        V.update_reagent('H2(aq)', 35)
        V.unify_reaction(SR.respiration.net_pathway, overwrite=False)

        # this method updates the quotient, std gibbs, and molar gibbs in
        # the current environment.
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        minG.append(SR.respiration.net_pathway.molar_gibbs)

        #for midG, use recorded values
        V.update_reagent('H2(aq)', 25)
        V.unify_reaction(SR.respiration.net_pathway, overwrite=False)

        # this method updates the quotient, std gibbs, and molar gibbs in
        # the current environment.
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        midG.append(SR.respiration.net_pathway.molar_gibbs)

        #for maxG, maximise products, minimise reactants
        V.update_reagent('H2(aq)', 15)
        V.unify_reaction(SR.respiration.net_pathway, overwrite=False)

        # this method updates the quotient, std gibbs, and molar gibbs in
        # the current environment.
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        maxG.append(SR.respiration.net_pathway.molar_gibbs)

    return Tlst, minG, midG, maxG

def Venus_overall_h2oxidation():
    """ maximum and minimum molar free energies from H2 oxidation """
    # setup specific reactor modelled on venusian cloud drops
    V=VenusDrop(H2ppm=25, workoutID=False)
    update_all_comp_ppm(V)

    SR = es.base_organism('H2 oxidiser', V,
      setup_h2oxidation(V), workoutID=False)

    # set up lists for plotting
    minminG=[]
    minG=[]
    midG=[]
    maxG=[]

    Tlst=np.linspace(278,350, num=100)

    for T in Tlst:
         ### 25C
        V.env.T = T
        update_all_comp_ppm(V)

        #for minG, maximise reactants, minimise products
        V.update_reagent('H2(aq)', 15)
        V.update_reagent('O2(aq)', 3)
        V.unify_reaction(SR.respiration.net_pathway, overwrite=False)

        # this method updates the quotient, std gibbs, and molar gibbs in
        # the current environment.
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        minminG.append(SR.respiration.net_pathway.molar_gibbs)

        #for minG, maximise reactants, minimise products
        V.update_reagent('H2(aq)', 35)
        V.update_reagent('O2(aq)', 68.8)
        V.unify_reaction(SR.respiration.net_pathway, overwrite=False)

        # this method updates the quotient, std gibbs, and molar gibbs in
        # the current environment.
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        minG.append(SR.respiration.net_pathway.molar_gibbs)

        #for midG, use recorded values
        V.update_reagent('H2(aq)', 25)
        V.update_reagent('O2(aq)', 43.6)
        V.unify_reaction(SR.respiration.net_pathway, overwrite=False)

        # this method updates the quotient, std gibbs, and molar gibbs in
        # the current environment.
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        midG.append(SR.respiration.net_pathway.molar_gibbs)

        #for maxG, maximise products, minimise reactants
        V.update_reagent('H2(aq)', 15)
        V.update_reagent('O2(aq)', 18.4)
        V.unify_reaction(SR.respiration.net_pathway, overwrite=False)

        # this method updates the quotient, std gibbs, and molar gibbs in
        # the current environment.
        SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

        maxG.append(SR.respiration.net_pathway.molar_gibbs)

    return Tlst, minG, midG, maxG, minminG

def plot_concentrations(ppmlst, concs_with_T, target_ppm, temps, reagent_name):

    fig, ax = plt.subplots()
    ax.plot(ppmlst, concs_with_T[0], 'g', label=str(temps[0])+' K')
    ax.plot(ppmlst, concs_with_T[1], 'b', label=str(temps[1])+' K')
    ax.plot(ppmlst, concs_with_T[2], 'r', label=str(temps[2])+' K')

    ax.axvline(target_ppm[0], ls='--', c='k')
    ax.axvline(target_ppm[1], c='k')
    ax.axvline(target_ppm[2], ls='--', c='k')

    plt.xlabel("Atmospheric Concentration of "+reagent_name+' [ppm]', fontsize=14)
    plt.ylabel("Concentration of Dissolved "+reagent_name+' [mol/L]', fontsize=14)
    ax.legend()
    plt.show()

def plot_energies(ppmlst, gibbs_with_T, target_ppm, temps, reaction_name, reagent_name):

    fig, ax = plt.subplots(figsize=(7,5), dpi=150)
    ax.plot(ppmlst, gibbs_with_T[0], c='#2c7bb6', label=str(temps[0])+' K', linewidth=2)
    ax.plot(ppmlst, gibbs_with_T[1], c='#fdae61', label=str(temps[1])+' K', linewidth=2)
    ax.plot(ppmlst, gibbs_with_T[2], c='#d7191c', label=str(temps[2])+' K', linewidth=2)

    ax.axvline(target_ppm[0], ls='--', c='k', lw=2)
    ax.axvline(target_ppm[1], c='k', lw=2)
    ax.axvline(target_ppm[2], ls='--', c='k', lw=2)

    plt.xlabel("Atmospheric concentration of "+reagent_name+" [ppm]", fontsize=14)
    plt.ylabel("Free energy [J/mol]", fontsize=14)
    ax.legend()
    plt.tight_layout()
    plt.savefig(reaction_name+'_'+reagent_name+'.tif')



##### tuples of target ppms to use
H2targets=(15,25,35)
O2targets=(18.4,43.6,68.8)
CO2targets=(750000, 940000,1000000)
CH4targets=(0.1, 5, 10)

##############   Methanogenesis vary CO2   ##############
#"""
# CO2 conc range, and corresponding gibbs range. Keeps H2 at 25 ppm and CH4 at 5 ppm
CO2ppmlst, CO2concs, CO2gibbs = Venus_atm_energy('Methanogen', setup_methanogenesis, 'CO2(aq)', temps=(278, 314, 350))

plot_concentrations(CO2ppmlst, CO2concs, CO2targets, (278, 314, 350), r'$\mathregular{CO}_2$')
plot_energies(CO2ppmlst, CO2gibbs, CO2targets, (278, 314, 350), 'methanogenesis', r'$\mathregular{CO}_2$')
#"""

##############   Methanogenesis vary H2   ##############
"""
# H2 conc range, and corresponding gibbs range. Keeps CO2 at 94% and CH4 at 5 ppm
H2ppmlst, H2concs, H2gibbs = Venus_atm_energy('Methanogen', setup_methanogenesis, 'H2(aq)', temps=(278, 314, 350))

plot_concentrations(H2ppmlst, H2concs, H2targets, (278, 314, 350), r'$\mathregular{H}_2$')
plot_energies(H2ppmlst, H2gibbs, H2targets, (278, 314, 350), 'methanogenesis', r'$\mathregular{H}_2$')
"""

##############   Methanogenesis vary CH4   ##############
"""
CH4ppmlst, CH4concs, CH4gibbs = Venus_atm_energy('Methanogen', setup_methanogenesis, 'CH4(g)', temps=(278, 314, 350))
plot_concentrations(CH4ppmlst,CH4concs, CH4targets, (278, 314, 350), r'$\mathregular{CH}_4$')
plot_energies(CH4ppmlst, CH4gibbs, CH4targets, (278, 314, 350), 'methanogenesis', r'$\mathregular{CH}_4$')
"""

##############   Sulfate Reduction   ##############
"""
SRH2ppmlst, SRH2concs, SRH2gibbs = Venus_atm_energy('Sulfate Reducer', setup_sulfatereduction, 'H2(aq)', temps=(278, 314, 350))

plot_concentrations(SRH2ppmlst, SRH2concs, H2targets, (278, 314, 350), r'\mathregular{H}_2$')
plot_energies(SRH2ppmlst, SRH2gibbs, H2targets, (278, 314, 350), 'sulfate reduction', r'$\mathregular{H}_2$')
"""

##############   H2 oxidation vary O2   ##############
"""
# O2 conc range, and corresponding gibbs range. Keeps H2 at 25 ppm
O2ppmlst, O2concs, O2gibbs = Venus_atm_energy('H2 oxidiser', setup_h2oxidation, 'O2(aq)', temps=(278, 314, 350))

plot_concentrations(O2ppmlst, O2concs, O2targets, (278, 314, 350), r'$\mathregular{O}_2$')
plot_energies(O2ppmlst, O2gibbs, O2targets, (278, 314, 350), 'H2 oxidation', r'$\mathregular{O}_2$')
# """

##############    H2 oxidation vary H2   ##############
"""
# H2 conc range, and corresponding gibbs range. Keeps O2 at 43.6 ppm
H2ppmlst, H2concs, H2gibbs = Venus_atm_energy('H2 oxidiser', setup_h2oxidation, 'H2(aq)', temps=(278, 314, 350))

plot_concentrations(H2ppmlst, H2concs, H2targets,(278, 314, 350), r'$\mathregular{H}_2$')
plot_energies(H2ppmlst, H2gibbs, H2targets, (278, 314, 350), r'$\mathregular{H}_2$ oxidation', r'$\mathregular{H}_2$')
"""

def overall_plot():
    Tlst, MminG, MmidG, MmaxG = Venus_overall_methanogenesis()
    Tlst, SminG, SmidG, SmaxG = Venus_overall_sulfatereduction()
    Tlst, HminG, HmidG, HmaxG, HminminG = Venus_overall_h2oxidation()

    fig, ax = plt.subplots(figsize=(7,5), dpi=150)



    ax.plot(Tlst, MmidG, c='#d7191c', label='Methanogenesis', lw=2)
    ax.fill_between(Tlst, MminG, MmaxG, color='#d7191c', alpha=0.5)

    ax.plot(Tlst, SmidG, c='#fdae61', label='Sulfate Reduction', lw=2)
    ax.fill_between(Tlst, SminG, SmaxG, color='#fdae61', alpha=0.5)

    ax.plot(Tlst, HmidG, c='#2c7bb6', label=r'$\mathregular{H}_2$ oxidation', lw=2)
    ax.fill_between(Tlst, HminG, HmaxG, color='#2c7bb6', alpha=0.5)

    ax.plot(Tlst, HminminG, c='tab:olive', label=r'$\mathregular{H}_2$ oxidation ($\mathregular{O}_2$ 3ppm)', lw=2)

    ax.axhline(0, c='k', lw=0.5)



    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Free energy of reaction [J/mol]')
    # ax.set_title("Free energy of metabolic reactions within composition uncertainties \n in the Venusian atmosphere's temperate region")
    ax.legend()
    plt.tight_layout()
    plt.savefig('summary.tif')

overall_plot()
