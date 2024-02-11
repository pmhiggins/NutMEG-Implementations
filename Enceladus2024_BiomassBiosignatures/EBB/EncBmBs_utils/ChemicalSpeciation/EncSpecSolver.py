import EncSpecSetup as ESSetup
import numpy as np
import matplotlib.pyplot as plt
import reaktoro as rkt
import pandas as pd
import sys, os

# rkt.Warnings.disable(525)

_models = {
  'pitzerPHREEQC': {
    'db' : rkt.PhreeqcDatabase('phreeqc.dat'),
    'ActivityMethod' : rkt.ActivityModelPitzerHMW(),
    'species' : {
      'input_aq' :rkt.speciate("H O Na Cl C"),
      'input_gases' : '',
      'full_rkt' : ["CO3-2", "H+", "H2O", "CO2", "(CO2)2", "HCO3-","CH4","Cl-", "H2", "Na+", "NaCO3-", "NaHCO3", "OH-", "NaOH","CH4(g)", "CO2(g)", "H2(g)","H2O(g)", "Halite", "O2", "O2(g)"], #O2 and O2(g)
      'full_nm' : ["CO3-2", "H+", "H2O", "CO2(aq)", "(CO2)2", "HCO3-","CH4(aq)","Cl-", "H2", "Na+", "NaCO3-", "NaHCO3", "OH-", "NaOH","CH4(g)", "CO2(g)", "H2(g)","H2O(g)", "NaCl", "O2", "O2(g)"],
      'H2O': 'H2O',
      'H+':'H+',
      'OH-':'OH-',
      'Na+':'Na+',
      'Cl-':'Cl-',
      'HCO3-':'HCO3-',
      'CO2':'CO2',
      'CO3-2':'CO3-2',
      'CH4':'CH4',
      'H2':'H2',
      'NaOH' : 'NaOH',
      'CO2(g)':'CO2(g)',
      'H2(g)':'H2(g)',
      'CH4(g)':'CH4(g)'
      }
    },
  'pitzerPHREEQCnoGases': {
    'db' : rkt.PhreeqcDatabase('phreeqc.dat'),
    'ActivityMethod' : rkt.ActivityModelPitzerHMW(),
    'species' : {
      'input_aq' :"H2O H+ OH- Na+ Cl- HCO3- CO2 CO3-2 NaCO3- NaHCO3 NaOH",
      'input_gases' : '',
      'full_rkt' : ["CO3-2", "H+", "H2O", "CO2", "HCO3-","Cl-", "Na+", "NaCO3-", "NaHCO3", "OH-", "NaOH", "Halite"],
      'full_nm' : ["CO3-2", "H+", "H2O", "CO2(aq)", "HCO3-","Cl-", "Na+", "NaCO3-", "NaHCO3", "OH-", "NaOH", "NaCl"],
      'H2O': 'H2O',
      'H+':'H+',
      'OH-':'OH-',
      'Na+':'Na+',
      'Cl-':'Cl-',
      'HCO3-':'HCO3-',
      'CO2':'CO2',
      'CO3-2':'CO3-2',
      'H2':'H2',
      'NaOH' : 'NaOH'
      }
    },
  'SUPCRT' : {
    'db' : rkt.SupcrtDatabase('supcrt07'),
    'ActivityMethod' : rkt.ActivityModelHKF(),
    'species' : {
      'input_aq' :"H2O(aq) H+ OH- Na+ Cl- HCO3- CO2(aq) CO3-2 Methane(aq) H2(aq) NaCl(aq) NaOH(aq)",
      'input_gases': '',#'CH4(g) H2(g) CO2(g) O2(g) H2O(g)',
      'full_rkt': ["H2O(aq)", "H+", "OH-", "Na+", "Cl-", 'HCO3-','CO2(aq)', 'CO3-2', 'Methane(aq)', 'H2(aq)', 'NaCl(aq)', 'NaOH(aq)', 'CH4(g)', 'H2(g)', 'CO2(g)', 'O2(g)', 'H2O(g)'],
      'full_nm': ["H2O", "H+", "OH-", "Na+", "Cl-", 'HCO3-','CO2(aq)', 'CO3-2', 'CH4(aq)', 'H2(aq)', 'NaCl(aq)', 'NaOH(aq)', 'CH4(g)', 'H2(g)', 'CO2(g)',  'O2(g)', 'H2O(g)'],
      'H2O': 'H2O(aq)',
      'H+':'H+',
      'OH-':'OH-',
      'Na+':'Na+',
      'Cl-':'Cl-',
      'HCO3-':'HCO3-',
      'CO2':'CO2(aq)',
      'CO3-2':'CO3-2',
      'CH4':'Methane(aq)',
      'H2':'H2(aq)',
      'NaCl' : 'NaCl(aq)',
      'NaOH' : 'NaOH(aq)',
      'CO2(g)':'CO2(g)',
      'H2(g)':'H2(g)',
      'CH4(g)':'CH4(g)'
      }
    },
  'SUPCRTnoGases' : {
    'db' : rkt.SupcrtDatabase('supcrt07'),
    'ActivityMethod' : rkt.ActivityModelHKF(),
    'species' : {
      'input_aq' :"H2O(aq) H+ OH- Na+ Cl- HCO3- CO2(aq) CO3-2 NaCl(aq) NaOH(aq)",
      'input_gases': '',
      'full_rkt': ["H2O(aq)", "H+", "OH-", "Na+", "Cl-", 'HCO3-','CO2(aq)', 'CO3-2', 'NaCl(aq)', 'NaOH(aq)'],
      'full_nm': ["H2O", "H+", "OH-", "Na+", "Cl-", 'HCO3-','CO2(aq)', 'CO3-2', 'NaCl(aq)', 'NaOH(aq)'],
      'H2O': 'H2O(aq)',
      'H+':'H+',
      'OH-':'OH-',
      'Na+':'Na+',
      'Cl-':'Cl-',
      'HCO3-':'HCO3-',
      'CO2':'CO2(aq)',
      'CO3-2':'CO3-2',
      'NaCl' : 'NaCl(aq)',
      'NaOH' : 'NaOH(aq)',
      }
    }
}


# savedirectory = 'rkt_Mar2023'
savedirectory = 'wide_pH'



def getAqueousPhase(T, pH_bo, salt, P=1e5, model='pitzerPHREEQC', fixGases=True, stdout=False,
T_bo=273.15, P_bo=1e5):
    """
    For set T pH and salt level perform a quick Enceladus speciation
    mimicing the one in Higgins et al 2021
    """

    # bulk ocean first estimates from H+C 2021
    mCO2_bo, aH2O_bo, _pH_bo, mDIC_bo, mCl_bo, mNa_bo = ESSetup.get_setup_spec_params(T_bo, pH_bo, salt)

    # mCO3, pH_bo, mDIC_bo, mCl_bo, mNa, mK = 0.002, 12., 0.03, 0.55, 0.1, 0.05

    # nominal plume molar ratios for the first (top) speciation
    # not important when the H2 and CH4 dissolved gases are not simulated.
    mCH4_bo = (0.2 / 0.55) * mCO2_bo
    mH2_bo =  (0.9 / 0.55) * mCO2_bo

    db = _models[model]['db']
    solution = None
    if model == 'pitzerPHREEQC':
        solution = rkt.AqueousPhase(_models[model]['species']['input_aq'], rkt.exclude("organic"))
    else:
        solution = rkt.AqueousPhase(_models[model]['species']['input_aq'])
    solution.setActivityModel(_models[model]['ActivityMethod'])

    minerals = rkt.MineralPhases()

    # where properties are missing in the database e.g. for saturation
    # approximate them as He. This is mostly needed for organics and SUPCRT
    critprops = rkt.CriticalProps()
    critprops.setMissingAs('He')

    system = None

    if  model != 'pitzerPHREEQCnoGases' and model != 'SUPCRTnoGases':
        gases = rkt.GaseousPhase(_models[model]['species']['input_gases'])
        gases.setActivityModel(rkt.ActivityModelPengRobinson())

        # the SUPCRT mdoels do not like minerals
        if model != 'SUPCRT':
            system = rkt.ChemicalSystem(db, solution, gases, minerals)
        else:
            system = rkt.ChemicalSystem(db, solution, gases)
    else:
        if model != 'SUPCRTnoGases':
            system = rkt.ChemicalSystem(db, solution, minerals)
        else:
            system = rkt.ChemicalSystem(db, solution)


    specs = rkt.EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.pH()
    specs.charge()
    specs.openTo(_models[model]['species']['Na+'])

    solver = rkt.EquilibriumSolver(specs)

    # setup the initial 'guess' at the top of the ocean based on Higgins et al 2021
    state = rkt.ChemicalState(system)
    state.temperature(T_bo, "kelvin")
    state.pressure(P_bo, "Pa")
    state.set(_models[model]['species']["H2O"], 1., "kg")
    # state.set(_models[model]['species']["Na+"], mNa_bo, "mol")
    state.set(_models[model]['species']["Cl-"], mCl_bo, "mol")
    # state.set(_models[model]['species']["CO2"], mCO2_bo, "mol")
    state.set(_models[model]['species']["HCO3-"], mDIC_bo, "mol")

    # state.set(_models[model]['species']["CO3-2"], mDIC_bo-mCO2_bo, "mol")
    state.set(_models[model]['species']["H+"], 10**-pH_bo, "mol")
    if model != 'pitzerPHREEQCnoGases' and model != 'SUPCRTnoGases':
        # only SUPCRT and pitzer with PHREEQC can include these species
        #if statement may no longer be needed if we are not using stock pitzer
        state.set(_models[model]['species']["CH4"], mCH4_bo, "mol")
        state.set(_models[model]['species']["H2"], mH2_bo, "mol")

        # have NO gaseous concs (all aqueous)
        state.set(_models[model]['species']["CH4(g)"], 1e-16, "mol")
        state.set(_models[model]['species']["H2(g)"], 1e-16, "mol")
        state.set(_models[model]['species']["CO2(g)"], 1e-16, "mol")



    if fixGases and model != 'pitzerPHREEQCnoGases' and model != 'SUPCRTnoGases':
        # prevent the gases from reacting
        specs.activity(_models[model]['species']['CH4(g)'])
        specs.activity(_models[model]['species']['H2(g)'])
        specs.activity(_models[model]['species']['CO2(g)'])


    aprops = rkt.AqueousProps(state)

    if stdout:
        print('############# Initial Conditions, 273 K, 1 Bar #############')
        print(aprops)


    # now set the conditions for the speciation
    conditions_bo = rkt.EquilibriumConditions(specs)
    conditions_bo.temperature(T_bo, "kelvin")
    conditions_bo.pressure(P_bo, "Pa")
    conditions_bo.pH(pH_bo)
    conditions_bo.charge(0.)


    restrictions = rkt.EquilibriumRestrictions(system)

    if fixGases and model != 'pitzer' and model != 'pitzerPHREEQCnoGases' and model != 'SUPCRTnoGases':
        # add a restrictions on the gases again
        restrictions.cannotReact(_models[model]['species']['CH4'])
        restrictions.cannotReact(_models[model]['species']['H2'])

        restrictions.cannotReact(_models[model]['species']['CH4(g)'])
        restrictions.cannotReact(_models[model]['species']['H2(g)'])
        restrictions.cannotReact(_models[model]['species']['CO2(g)'])

    # uncomment this to close the system to Na+
    # restrictionsT.cannotReact(_models[model]['species']['Na+'])

    # do the bulk ocean specation
    solver.solve(state, conditions_bo, restrictions)
    aprops = rkt.AqueousProps(state)

    if stdout:
        print('############# First Speciation, 273 K, 1 Bar #############')
        print(aprops)

    # now we have [Na+] and charge balance, open to OH and respeciate at
    # your updated conditions
    specs_T = rkt.EquilibriumSpecs(system)
    specs_T.temperature()
    specs_T.pressure()
    specs_T.charge()
    specs_T.openTo('OH-')

    solver_T = rkt.EquilibriumSolver(specs_T)

    # continue trying to keep these constant
    if fixGases and model != 'pitzerPHREEQCnoGases' and model != 'SUPCRTnoGases':
        # continue trying to keep these constant
        specs_T.activity(_models[model]['species']['CH4(g)'])
        specs_T.activity(_models[model]['species']['H2(g)'])
        specs_T.activity(_models[model]['species']['CO2(g)'])


    # specs.closedTo('Na+')
    # specs.openTo('OH-')

    conditions_T = rkt.EquilibriumConditions(specs_T)
    conditions_T.temperature(T, "kelvin")
    conditions_T.pressure(P, "Pa")
    # conditionsT.pH(pH)
    conditions_T.charge(0.)

    restrictions_T = rkt.EquilibriumRestrictions(system)

    if fixGases and model != 'pitzerPHREEQCnoGases' and model != 'SUPCRTnoGases':
        restrictions_T.cannotReact(_models[model]['species']['CH4'])
        restrictions_T.cannotReact(_models[model]['species']['H2'])

        restrictions_T.cannotReact(_models[model]['species']['CH4(g)'])
        restrictions_T.cannotReact(_models[model]['species']['H2(g)'])
        restrictions_T.cannotReact(_models[model]['species']['CO2(g)'])


    solver_T.solve(state, conditions_T, restrictions_T)

    # now get the aqueous properties (molal etc.)
    aprops = rkt.AqueousProps(state)
    props = state.props()

    if stdout:
        print("############# SPECIATED STATE AQUEOUS PROPERTIES ############# ")
        print(aprops)
        # print("############# ALL PROPERTIES ############# ")
        # print(state.props())

    return aprops, props




def get_DIC(T, pH_bo, salt, model='pitzerPHREEQC'):
    """
    return a dict of DIC molal properties at given T, 273 pH and
    salt concentration. DIC here is the sumo of dissolved carbonate species
    CO2, HCO3-, CO3-2.

    """
    ap, p = getAqueousPhase(T, pH_bo, salt, P=1e5, model=model)
    _CO2 = float(ap.speciesMolality(_models[model]['species']['CO2']))
    _HCO3 = float(ap.speciesMolality(_models[model]['species']['HCO3-']))
    _CO3 = float(ap.speciesMolality(_models[model]['species']['CO3-2']))

    return {'DIC': _CO2 + _HCO3 + _CO3,
      'CO2' : _CO2,
      'HCO3' : _HCO3,
      'CO3' : _CO3}





def get_spec_params(T, pH_bo, salt, P=1e5, model='pitzerPHREEQC', fixGases=False, stdout=False):
    """
    return a dict of all relevant speciated properties at
    given T, bulk ocean pH (ie at 273K) and salt concentration.

    If the speciation fails, return a dictionary of the speciation parameters
    populated with np.nan

    """

    _d = {}

    try:
        # perform speciation, get AqueousPhase properties
        ap, p = getAqueousPhase(T, pH_bo, salt, P=P, model=model,
          fixGases=fixGases, stdout=stdout)

        print(T, pH_bo, salt, ' successful' )

        # populate a dictionary with the speciation information
        for species, db_species in zip(
          _models[model]['species']['full_nm'],
          _models[model]['species']['full_rkt']):
            _d['m'+species] = float(p.speciesConcentration(db_species))
            _d['g'+species] = float(p.speciesActivityCoefficient(db_species))
            _d['a'+species] = float(p.speciesActivity(db_species))

        _d['mDIC'] = _d['mCO2(aq)'] + _d['mHCO3-'] + _d['mCO3-2']

        _d['I'] = float(ap.ionicStrength())
        _d['pE'] = float(ap.pE())
        _d['Eh'] = float(ap.Eh())
        _d['pH'] = float(ap.pH())

        _d['logSI_CO2g'] = float(ap.saturationIndexLg('CO2(g)'))
        _d['logSI_H2Og'] = float(ap.saturationIndexLg('H2O(g)'))
        _d['logSI_O2g'] = float(ap.saturationIndexLg('O2(g)'))
        _d['logSI_CH4g'] = float(ap.saturationIndexLg('CH4(g)'))
        _d['logSI_H2g'] = float(ap.saturationIndexLg('H2(g)'))
        _d['logSI_NaCl'] = float(ap.saturationIndexLg('Halite'))

    except:
        # for fringe cases where the speciation fails, return a dict of np.nan
        # as not to disrupt a series of speciations.

        print(model, P, T, pH_bo, salt, ' failed, skipping' )

        for species, db_species in zip(_models[model]['species']['full_nm'], _models[model]['species']['full_rkt']):
            _d['m'+species] = np.nan
            _d['g'+species] = np.nan
            _d['a'+species] = np.nan

        _d['mDIC'] = _d['mCO2(aq)'] + _d['mHCO3-'] + _d['mCO3-2']

        _d['I'] = np.nan
        _d['pE'] = np.nan
        _d['Eh'] = np.nan
        _d['pH'] = np.nan

        _d['logSI_CO2g'] = np.nan
        _d['logSI_H2Og'] = np.nan
        _d['logSI_NaCl'] = np.nan
        _d['logSI_O2g'] = np.nan
        _d['logSI_CH4g'] = np.nan
        _d['logSI_H2g'] = np.nan


    return _d
