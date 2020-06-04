"""
Here we're going to simulate some competition between two organisms - a
methanogen and a sulfate reducer - and see what changing some of the model
parameters does to their growth curves.

There's a bug when you run the lifespan part on the first run, this is something
to do with stopping simulations and I will get to it (one day...) It won't affect
results much, so just run the code again after getting the error and it wil be fine.
"""

# set std_dbpath to modelparams_draft

import sys, os, ast
sys.path.append(os.path.dirname(__file__)+'../../')

import NutMEG as nm
import NutMEG.reaction as reaction
import NutMEG.plotter as nutplt
import matplotlib.pyplot as plt
import pandas as pd
import sqlite3
import numpy as np

# choose a path and filename to a database to save model output. This will save
# having to rerun potentially many simulations again just for replotting.
# By default, NutMEG creates a database called NutMEG_db in the directory of
# your code, but let's make a custom one for this.
dbpath='modelparams'

# This method sets up the new database. Note we have to pass the dbpath whenever
# we use it, otherwise NutMEG will use the default.
nm.db_helper.create_major_tables(replace=False, dbpath=dbpath)


def setup_methanogenesis(R, k_RTP=0.0001):
    """Create a simple methanogenesis reaction as a standard thermodynamic
    reaction inside reactor R.
    The reagents have no special properties, we can add those later
    """

    # introduce some reagents
    CO2 = reaction.reagent('CO2(aq)', R.env, phase='aq')
    H2aq = reaction.reagent('H2(aq)', R.env, phase='aq')
    CH4aq = reaction.reagent('CH4(g)', R.env, phase='g')
    H2O = reaction.reagent('H2O(l)', R.env, phase='l')

    # the overall reaction is CO2 + 4H2 -> CH4 + 2H2O
    thermalMG = reaction.reaction({CO2:1, H2aq:4}, {CH4aq:1, H2O:2},
      R.env)

    # we need to manually set the RTP rate constant for use in biochemistry
    thermalMG.rate_constant_RTP = k_RTP

    return thermalMG

def setup_sulfatereduction(R, k_RTP=0.0001):
    """Create a simple sulfate reduction reaction as a standard thermodynamic
    reaction inside reactor R.
    The reagents have no special properties, we can add those later
    """

    H2aq = reaction.reagent('H2(aq)', R.env, phase='aq')
    SO4 = reaction.reagent('SO4--', R.env, phase='aq', charge=-2)
    H = reaction.reagent('H+', R.env, phase='aq', charge=1)
    HS = reaction.reagent('HS-', R.env, phase='aq', charge=-1)
    H2O = reaction.reagent('H2O(l)', R.env, phase='l')

    # the overall reaction is 4H2 + SO4(2-) + H(+) -> HS(-) + 4H2O
    thermalSR = reaction.reaction({H2aq:4, SO4:1, H:1}, {HS:1, H2O:4},
      R.env)

    # we need to manually set the RTP rate constant for use in biochemistry
    thermalSR.rate_constant_RTP = k_RTP

    return thermalSR


def initial_conditions(R, comp={}):
    """Set up the reactor R to use by populating it with reagents. To change
    the concentrations used pass then in the dict comp in the  format
    {Name : conc}
    """
    # metabolic concentrations
    mol_CO2 = comp.pop('CO2(aq)', 0.001)
    mol_CH4 = comp.pop('CH4(g)', 1e-7)
    mol_H2 =  comp.pop('H2(aq)', 0.001)
    mol_SO4 = comp.pop('SO4--', 0.001)
    mol_HS = comp.pop('HS-', 1e-7)

    # life also needs a source of N and P
    mol_NH3 = comp.pop('NH3(aq)', 0.001)
    mol_H2PO4 = comp.pop('H2PO4-', 0.001)

    # concentration of H is tied to pH.
    mol_H=10**(-R.pH)

    # now set up as reagents
    CO2 = reaction.reagent('CO2(aq)', R.env, phase='aq', conc=mol_CO2,
      activity=mol_CO2)
    H2aq = reaction.reagent('H2(aq)', R.env, phase='aq', conc=mol_H2,
      activity=mol_H2)
    CH4aq = reaction.reagent('CH4(g)', R.env, phase='g', conc=mol_CH4,
      activity=mol_CH4)
    H2O = reaction.reagent('H2O(l)', R.env, phase='l', conc=55.5,
      activity=1, phase_ss=True)
    H = reaction.reagent('H+', R.env, charge=1, conc=mol_H,
      phase='aq', activity=mol_H)
    SO4 = reaction.reagent('SO4--', R.env, phase='aq', charge=-2,
      conc=mol_SO4, activity= mol_SO4)
    HS = reaction.reagent('HS-', R.env, phase='aq', charge=-1, conc=mol_HS,
      activity=mol_HS)

    # add these to the composition of R. This will be shared between
    # the organisms
    R.composition = {CO2.name:CO2, H2aq.name:H2aq,
      CH4aq.name:CH4aq, H2O.name:H2O, H.name:H, SO4.name:SO4, HS.name:HS}

    # add in the extra nutrients
    R.composition['NH3(aq)'] = reaction.reagent('NH3(aq)', R.env,
      phase='aq', conc=mol_NH3, activity=mol_NH3)
    R.composition['H2PO4-'] = reaction.reagent('H2PO4-', R.env, phase='aq',
      conc=mol_H2PO4, activity=mol_H2PO4)


def simulate_competition(comp={}, k=[0.001,5000], inflow={'H+':0},
  methanogen_changes={}, SR_changes={}, reactor_changes={}):
  # nATP=[1,1],
  # life=[float('inf'),float('inf')], maint=[0,0], inflow={'H+':0},
  # Pconc=1e-6, Puptakes=(1e-10,1e-10)):

    # create a reactor object. reactor names relate them to where they get
    # saved in the database, so let's just call this reactor1.
    # passing workoutID=False means it doesn't save parameters to the database
    # yet. We'll make some changes and do so later.
    R = nm.reactor('reactor1', workoutID=False, pH=7.0, dbpath=dbpath, **reactor_changes)

    initial_conditions(R, comp=comp) # sets up composition.

    # # implement a source or sink of reagents.
    # for ke,v in inflow.items():
    #     R.composition_inputs[ke] = v

    # create a horde or methanogens, and give it any unique parameters
    # You could also create a base_organism and put that in a colony for a
    # slightly more accurate population, but expect that to take much longer
    # and probably be buggy as it's not been thoroughly tested.
    H = nm.horde('Methanogen', R, setup_methanogenesis(R, k_RTP=k[0]), 10,
      Tdef='None', dbpath=dbpath, **methanogen_changes)
    # Tdef shows which adapations against temperature to use, here we'll ignore
    # it and put in our own maitnenance powers.

    # create a horde of sulfate reducers, and give it any unique parameters
    H2 = nm.horde('SulfateReducer', R, setup_sulfatereduction(R, k_RTP=k[1]), 10,
      Tdef='None', dbpath=dbpath, **SR_changes)


    # setting up the organisms adds their metabolisms to R so it's ready to save!
    R.dbh.workoutID()

    # put both hordes in a culture object, which keeps organisms together.
    Cu = nm.culture(hordes=[H2, H])

    # put the culture and the reactor together. Now we have an ecosystem ready
    # to go!
    ES = nm.ecosystem(R, Cu, dbpath=dbpath)

    # this runs a growth prediction with the two hordes together in the reactor
    # with initial conditions as defined and saves the output in the database.
    # Various criteria cause the simulation
    # to stop checkout the ecosystem class to see what you can change. We'll
    # add a maximum time (equivalent to ~6yrs).
    ES.predict_growth(tmax=2e8)

    return ES.dbh.SimID


cSeries=['b', 'r']

def orgcurves(orgfig, SimIDs, param='no_alive', ax=None, ls='-'):
    """plot the ecosystem simulation parameter param against time for the simulation(s)
    in SimIDs. You can add any parameter from the Full_Results_Sim table.
    When you run a simulation, you'll see some options output every 100 steps.
    """

    TSeries = [] # to hold the time data
    paramSeries =[] # to hold the parameter
    num_orgs = 2 # number of organism types (for the colors)

    for S in SimIDs:
        OrgIDs, LocID = nm.ecosystem_dbhelper.db_helper.findOrgIDsLocID(S, dbpath=dbpath)
        T = nm.ecosystem_dbhelper.db_helper.extract_param_db_Sim(
          S, 'Time', dbpath=dbpath)
        OrgIDs = ast.literal_eval(OrgIDs)
        num_orgs = max(len(OrgIDs), num_orgs)

        for O in OrgIDs:
            paramSeries.append(
              nm.ecosystem_dbhelper.db_helper.extract_param_db_Sim(
                S, param+'_'+O, dbpath=dbpath))


        if T ==[0]:
            # this is a simulation which had no growth or an error, skip
            continue
        else:
            # we need as many T series as we have organisms.
            for O in range(len(OrgIDs)):
                TSeries.append(T)

    # this plots on the axis and passes the axes back. You could just use
    # matplotlib if you prefer
    return orgfig.linearplot(xvals=TSeries, yvals=paramSeries, ax=ax,
      colors=cSeries, show=False, ls=[ls])


def compcurves(orgfig, SimIDs, ax=None, ls='-'):
    """Plot some of the composition of the environment with time for the
    simualtion SimID.
    """

    TSeries = [] # times
    CO2, SO4, H2, H, CH4 = [],[],[],[], [] # list for some reagents
    LSeries = ['CO2', 'SO4--', 'H2', 'H+', 'CH4']# list of labels for them
    colSeries = ['b', 'r', 'g', 'c', 'k'] # colors for each line
    for S in SimIDs:
        T = nm.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S, 'Time', dbpath=dbpath)

        comp = tuple(nm.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S, 'Composition', dbpath=dbpath))
        # this is a list of the composition dictionary at each time step!
        # turn these into lists for each reagent for easy plotting
        for c in comp:
            CO2.append(ast.literal_eval(c[0])['CO2(aq)'])
            SO4.append(ast.literal_eval(c[0])['SO4--'])
            H2.append(ast.literal_eval(c[0])['H2(aq)'])
            H.append(ast.literal_eval(c[0])['H+'])
            CH4.append(ast.literal_eval(c[0])['CH4(g)'])

        if T ==[0]:
            # ignore a simulation which didn't work or had no growth
            continue
        else:
            # copy the time list as many times as we need...
            [TSeries.append(T) for i in range(len(LSeries))]


    return orgfig.linearplot(xvals=TSeries, yvals=[CO2, SO4, H2, H, CH4],
      ax=ax, labels=LSeries, colors=colSeries, show=False, ls=[ls])



fig, axs = plt.subplots(3, 2, figsize=(6,8))
orgfig = nutplt.growthparams(fig=fig)

Sim = simulate_competition(methanogen_changes={'Basal':1e-15}, SR_changes={'Basal':1e-15})
orgcurves(orgfig, [Sim], ax=axs[0][0])
compcurves(orgfig, [Sim], ax=axs[0][1])

axs[0][1].set_yscale('log')
axs[0][0].set_yscale('log')
axs[0][0].set_title('Control')
axs[0][1].set_title('Control Composition')
# Now let's change some model parameters to show how it can change the growth
# curves.

# let's increase the maintenance power of the methanogen, making survival harder
for MP in [2e-15, 5e-15, 1e-14, 5e-14]:
    Sim = simulate_competition(methanogen_changes={'Basal':MP}, SR_changes={'Basal':1e-15})
    orgcurves(orgfig, [Sim], ax=axs[1][0], ls='--')
axs[1][0].set_title('Maintenance Power')
axs[1][0].set_yscale('log')



# change the yield of ATP per mol metabolism - for the sulfate reducer
for n in np.linspace(0.5,1.5, num=5):
    Sim = simulate_competition(methanogen_changes={'Basal':1e-15}, SR_changes={'Basal':1e-15, 'n_ATP':n})
    orgcurves(orgfig, [Sim], ax=axs[1][1], ls='--')
axs[1][1].set_title('ATP Yield')


# let's add an inflow of H+
for h in np.logspace(-16,-12, num=5):
    Sim = simulate_competition(methanogen_changes={'Basal':1e-15}, SR_changes={'Basal':1e-15}, reactor_changes={'composition_inputs':{'H+':h}})
    orgcurves(orgfig, [Sim], ax=axs[2][0], ls='--')
axs[2][0].set_title('H+ inflow')



# give both organisms a life span which they can't exceed.
for l in [1e3,1e4,5e4]:
    Sim = simulate_competition(methanogen_changes={'Basal':1e-15, 'base_life_span':l}, SR_changes={'Basal':1e-15, 'base_life_span':l})
    orgcurves(orgfig, [Sim], ax=axs[2][1], ls='--')
axs[2][1].set_title('Vary lifespans')
axs[2][1].set_xlim(0,1e6)


# tidy up the axes a little
for ax in axs:
    for a in ax:
        a.set_yscale('log')
        a.set_xlabel('Time [s]')
plt.tight_layout()

plt.show()


#this prints out the tables in the database for the orgnaisms and reactor
# which might help you understand how the data is stored. When extracting
# data (e.g. the param to send to orgcurves), used one of the column names,
# or alternatively look in the base_organism_dbhelper or reactor_dbhelper classes
# for dbdict to see what you could use. Hopefully you'll recognise some from the
# code above!
nm.db_helper.print_table('Methanogen', dbpath=dbpath)
nm.db_helper.print_table('SulfateReducer', dbpath=dbpath)
nm.db_helper.print_table('reactor1', dbpath=dbpath)
nm.db_helper.print_table('Summary', dbpath=dbpath)
nm.db_helper.print_table('FullResults_Sim_'+str(Sim), dbpath=dbpath)


# why not experiment with different combinations?
