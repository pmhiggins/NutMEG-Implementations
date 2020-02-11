import sys, os, math, sqlite3, csv
from copy import deepcopy
import numpy as np
sys.path.append(os.path.dirname(__file__)+'../..')

import NutMEG as es
from NutMEG.reactor.saved_systems.VenusDrop import VenusDrop
import NutMEG.reaction as reaction

import matplotlib.pyplot as plt

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
      # to set up an organism.
    thermalSR.rate_constant_RTP = k_RTP

    return thermalSR


# setup specific reactor modelled on venusian cloud drops
V=VenusDrop(H2ppm=50, workoutID=False)

# initialise a basic sulfate reducer. It's an organism with sulfate
# redution as its net energy-yielding reaction. Passing V means we
# put it in our VenusDrop
SR = es.base_organism('SulfateReducer', V,
  setup_sulfatereduction(V), workoutID=False)

# set up lists for plotting

# 25C
molH2_25 = [] # conc H2
Q_25 = [] # reaction quotient
K_25 = [] # equilibrium constant
stdG_25 = [] # standard gibbs per mol
G_25 = [] #Available energy per mol

# do the same for 75C
molH2_75=[]
Q_75 =[]
K_75 = []
stdG_75 = []
G_75 =[]

# iterate over H2 ppm values
H2ppmlst=np.linspace(0.01, 200, num=101)

for ppm in H2ppmlst:
    ### 25C
    V.env.T = 273+25

    # find conc at this ppm and update
    V.update_conditions(ppm, 1e-10)

    # this method updates the quotient, std gibbs, and molar gibbs in
    # the current environment.
    SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

    # extract the parameters we're interested in from the organism
    molH2_25.append(V.composition['H2(aq)'].conc)
    Q_25.append(SR.respiration.net_pathway.quotient)
    K_25.append(math.exp(SR.respiration.net_pathway.lnK))
    stdG_25.append(SR.respiration.net_pathway.std_molar_gibbs)
    G_25.append(SR.respiration.net_pathway.molar_gibbs)

    ### 75C
    V.env.T = 273+75
    # repeat updates for new temperature
    V.update_conditions(ppm, 1e-10)
    SR.respiration.net_pathway.update_molar_gibbs_from_quotient()

    molH2_75.append(V.composition['H2(aq)'].conc)
    Q_75.append(SR.respiration.net_pathway.quotient)
    K_75.append(math.exp(SR.respiration.net_pathway.lnK))
    stdG_75.append(SR.respiration.net_pathway.std_molar_gibbs)
    G_75.append(SR.respiration.net_pathway.molar_gibbs)



# plot the H2 concentration with ppm
plt.plot(H2ppmlst, molH2_25)
plt.plot(H2ppmlst, molH2_75)
plt.show()

# plot the free energy available with ppm
plt.plot(H2ppmlst, G_25)
plt.plot(H2ppmlst, G_75)
plt.show()
