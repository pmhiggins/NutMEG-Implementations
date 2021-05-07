import sys, os, ast, math, statistics
sys.path.append(os.path.dirname(__file__)+'../../NutMEG')
sys.path.append(os.path.dirname(__file__)+'../../NutMEG-Implementations/TOM')
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp

# imports from NutMEG
import NutMEG as es
from NutMEG.reactor.saved_systems.Enceladus import Enceladus as Enc
from NutMEG.culture.saved_organisms.TypicalOptimalMethanogen import TypicalOptimalMethanogen as TOM

# other classes relevant to this project
import theory_emp_match as tem
import QuotientUncertainties as Qunc
from EnergyCalculations import getPS
from EnergyCalculations import getE_TOM



es.db_helper.create_major_tables(replace=False)


gEE_default_output = ['Gibbs_Methanogenesis', 'PowerSupply']
# extra options are:
#   'Quotient', 'Composition', 'ATPGibbs', 'max_k'

default_sigmaorder = ['CO2s', 'H2s', 'CH4s', 'nATPs',  'k_corr', 'T', 'pH']
default_nomdict = {'T':300,'pH':8.5, 'CH4':0., 'CO2':0., 'H2':0., 'nATP':1., 'k_corr':0.0}


def boolify_output(requests):
    """
    From a list of requests, build a dictionary for the outputs to present.
    """
    output_options = ['Gibbs_Methanogenesis', 'PowerSupply', 'Quotient',
      'Composition', 'ATPGibbs', 'max_k', 'CO2Tiger']
    outputdict={}
    for r in output_options:
        if r in requests:
            outputdict[r] = True
        else:
            outputdict[r] = False
    return outputdict

def getPSfromSigmas(sigmas, sigmaorder=default_sigmaorder, dummyvals=False,
  salttype='nom',):
    """
    From the sigmas of the shape returned in getsimgas(), get a list of power
    supplies for those sigma values.
    """
    CO2s, H2s, CH4s, nATPs, k_corr, Ts, pHs = sigmas

    E = Enc(name='En', pH=pHs, CO2origin='HTHeatingSalts', T=Ts, workoutID=True)

    E = Qunc.get_salty_Enc(E, salttype)

    E.composition['CO2(aq)'].activity = E.composition['CO2(aq)'].activity.n + (CO2s * E.composition['CO2(aq)'].activity.s)
    E.composition['H2(aq)'].activity = E.composition['H2(aq)'].activity.n + (H2s * E.composition['H2(aq)'].activity.s)
    E.composition['Methane(aq)'].activity = E.composition['Methane(aq)'].activity.n + (CH4s * E.composition['Methane(aq)'].activity.s)

    # H2O errors not supported yet
    E.composition['H2O(l)'].activity = E.composition['H2O(l)'].activity.n

    return getPS(E, nATPs, k_corr)


def getsigmas(num, fixed={'T':None,'pH':None, 'CH4':None, 'CO2':None, 'H2':None, 'nATP':None, 'k_corr':None}, Tlims=(273,473), pHlims=(8,12), Encel=None, dummyvals=False):
    """
    get a uniform sample of the default_nomdict, unless some parameters are
    fixed (passed as the keys of fixed). For activities, these are only relative:
    e.g. between -1 and +1 of the uncertainty. For nATP these are the actual
    values of nATP. For kRTP they are the change in  log-space.

    Encel and dummyvals are from legacy and testing so are obsolete. But these
    are left in for legacy code which used them.
    """
        # for now Encel is a placeholder - I realised it doesn't work.
    Tmid = statistics.mean(Tlims)
    pHmid = statistics.mean(pHlims)
    pHerr = 0.5*(pHlims[1]-pHlims[0])
    Terr = 0.5*(Tlims[1]-Tlims[0])
    CH4mid, CO2mid, H2mid = 0,0,0
    CH4err, CO2err, H2err = 1.,1.,1.
    nATPmid, nATPerr = 1., 0.5
    k_corrmid, k_correrr = 0,1.

    nomdict = {'T':statistics.mean(Tlims),'pH':statistics.mean(pHlims), 'CH4':0., 'CO2':0., 'H2':0., 'nATP':1.125, 'k_corr':0.}
    errdict = {'T':Terr,'pH':pHerr, 'CH4':1., 'CO2':1., 'H2':1., 'nATP':0.875, 'k_corr':1.}
    for key, val in fixed.items():
        if isinstance(val, float) or isinstance(val, int):
            nomdict[key] = val
            errdict[key] = 0

    ndim = len(default_sigmaorder)

    nominal = np.array([nomdict['CO2'], nomdict['H2'],
      nomdict['CH4'],
      nomdict['nATP'],
      nomdict['k_corr'],
      nomdict['T'],
      nomdict['pH']])
    lowbounds = np.array([nomdict['CO2'] - errdict['CO2'],
      nomdict['H2'] - errdict['H2'],
      nomdict['CH4'] - errdict['CH4'],
      nomdict['nATP'] - errdict['nATP'],
      nomdict['k_corr'] - errdict['k_corr'],
      nomdict['T'] - errdict['T'],
      nomdict['pH'] - errdict['pH']])
    upbounds = np.array([nomdict['CO2'] + errdict['CO2'],
      nomdict['H2'] + errdict['H2'],
      nomdict['CH4'] + errdict['CH4'],
      nomdict['nATP'] + errdict['nATP'],
      nomdict['k_corr'] + errdict['k_corr'],
      nomdict['T'] + errdict['T'],
      nomdict['pH'] + errdict['pH']])

    # nb no gaussians, all flat
    return [nominal + (upbounds - lowbounds)*(-0.5+np.random.rand(ndim)) for i in range(num)]



def getEncEnergetics(ocean_pH=8., T=273.15, nATP=1.0, k_corr=0.,
  CO2origin='pH', # change to HTHeating later
  output=gEE_default_output,
  quotienttype='salty_endmember'):
    """
    Get the requested parameters in output at the selected parameters of
    ocean_pH, T, nATP, and k_corr.
    """

    outputbools = boolify_output(output)
    outputdict = {}

    E = Enc(name='En', pH=ocean_pH, CO2origin=CO2origin, T=T, workoutID=False)

    if outputbools['Composition']:
        _concs = {}
        for species, rx in E.composition.items():
            _concs[species] = np.array([[rx.activity.n,
              rx.activity.n + rx.activity.s,
              rx.activity.n - rx.activity.s]])
        outputdict['Composition'] = _concs

    if outputbools['CO2Tiger']:
        CO2T = E.get_tigerstripe_CO2(logform=True)
        CO2Tiger = 10**CO2T.n
        CO2Tiger_up = 10**(CO2T.n+CO2T.s)
        CO2Tiger_down = 10**(CO2T.n-CO2T.s)
        outputdict['CO2Tiger'] = np.array([[CO2Tiger, CO2Tiger_do, CO2Tiger_up]])

    if outputbools['Gibbs_Methanogenesis'] or outputbools['PowerSupply'] or outputbools['ATPGibbs'] or outputbools['max_k'] or outputbools['Quotient']:

        # need to do the energetic calculations. This is outsourced to the
        # QuotientUncertainties class.

        Qdict = {} # dictionary outputs from Qunc

        if quotienttype == 'allunc':
            Qdict = Qunc.Q_allunc(E, nATP, k_corr)
        elif quotienttype == 'salty_endmember':
            Qdict = Qunc.Q_salty_endmember(E, nATP, k_corr)
        elif quotienttype == 'salty_nominal':
            Qdict = Qunc.Q_salty(E, nATP, k_corr, 'nom')
        elif quotienttype == 'salty_high':
            Qdict = Qunc.Q_salty(E, nATP, k_corr, 'high')
        elif quotienttype == 'salty_low':
            Qdict = Qunc.Q_salty(E, nATP, k_corr, 'low')
        else:
            raise ValueError('unrecognised quotienttype, how do you want to consider salt effects?')

        if outputbools['Gibbs_Methanogenesis']:
            outputdict['Gibbs_Methanogenesis'] = Qdict['Gibbs_Methanogenesis']
        if outputbools['Quotient']:
            outputdict['Quotient'] = Qdict['Quotient']
        if outputbools['PowerSupply']:
            outputdict['PowerSupply'] = Qdict['PowerSupply']
        if outputbools['ATPGibbs']:
            outputdict['ATPGibbs'] = Qdict['ATPGibbs']
        if outputbools['max_k']:
            outputdict['max_k'] = Qdict['max_k']

    return outputdict


def getMesh(Trange, pHrange, params=['Gibbs_Methanogenesis'],
  nATP=1.0, k_corr=0.0, CO2origin='pH', quotienttype='salty_endmember'):
    """
    Get a dictionary of meshgrids, including the temperatures, ocean pH values,
    and entries in params.
    """
    Meshes = {p:[] for p in params}
    for p in params:
        op = np.ndarray((len(Trange), len(pHrange)))
        op_do = np.ndarray((len(Trange), len(pHrange)))
        op_up = np.ndarray((len(Trange), len(pHrange)))

        yindex = 0
        for T in Trange:
            xindex = 0
            for pH in pHrange:
                _theseparams = getEncEnergetics(ocean_pH=pH, T=T, nATP=nATP,
                  k_corr=k_corr,
                  CO2origin=CO2origin, # change to Glein later
                  output=params,
                  quotienttype=quotienttype)
                op[yindex][xindex] = _theseparams[p][0][0]
                op_do[yindex][xindex] = _theseparams[p][0][1]
                op_up[yindex][xindex] = _theseparams[p][0][2]

                xindex += 1
            yindex += 1
        Meshes[p] = np.array([op, op_do, op_up])

    Meshes['pH'], Meshes['T'] = np.meshgrid(pHrange, Trange)

    return Meshes


######### for maintenance estimates ###########

def get_maintenances(T,pH, E=None, TOM=None, nATP=1.0):
    """
    Get the four maintenance updates at T and pH.

    This may take time to run the first time, it will be quicker if you have
    the database set up in the TOM directory. You can also contact me for a
    copy of the full database.
    """
    PT=1
    if T <375:
        if nATP == 1.0:
            PT = tem.MaintenanceRange_nATPs(Trange=[int(T)], mCH4=3e-8, Tlst=False, fraction=False, Perform=False, dbpath='../TOM/allMtestc')[0][0]
        elif nATP == 0.5 or nATP==0.25:
            PT = tem.MaintenanceRange_nATPs(Trange=[int(T)], mCH4=3e-8, Tlst=False, fraction=False, Perform=False, dbpath='../TOM/allMtestc')[1][0]
            if PT <1e-20:
                #correct up for the v low estimates at near 0C
                PT=1e-12
        elif nATP > 1.5:
            PT = tem.MaintenanceRange_nATPs(Trange=[int(T)], mCH4=3e-8, Tlst=False, fraction=False, Perform=False, dbpath='../../NutMEG-Implementations/TOM/allMtestc')[2][0]

    # make a TOM
    if E==None or TOM==None:
        E = Enc('EncT')
        TOM = getE_TOM(E)
    # make an Enceladus
    TE = es.applications.theory_estimates(TOM, E)
    TE.loc.change_T(T)
    TE.org.get_ESynth(AA=True)
    td = TE.temperature_defenses(T)
    Ti = td['TijhuisAnaerobe']
        # L10.append(td['Lever10pc'])
    L2 = td['Lever2pc']
    return PT, Ti, L2

def maintenancemesh(Trange, pHrange, nATP=1.0):
    """Get a meshgrid of maintenance powers."""

    E = Enc('EncT')
    TOM = getE_TOM(E)

    HC = np.ndarray((len(Trange), len(pHrange)))
    Ti = np.ndarray((len(Trange), len(pHrange)))
    L = np.ndarray((len(Trange), len(pHrange)))
    B = np.ndarray((len(Trange), len(pHrange)))

    yindex = 0
    for T in Trange:
        xindex = 0
        for pH in pHrange:
            MPs = get_maintenances(T, pH, E, TOM, nATP=nATP)
            HC[yindex][xindex] = MPs[0]
            Ti[yindex][xindex] = MPs[1]
            L[yindex][xindex] = MPs[2]
            B[yindex][xindex] = 1e-18

            xindex += 1
        yindex += 1
    return HC, Ti, L, B
