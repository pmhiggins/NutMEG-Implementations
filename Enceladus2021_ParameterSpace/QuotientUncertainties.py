from copy import deepcopy
from EnergyCalculations import MGparams
from NutMEG.reactor.saved_systems.Enceladus import Enceladus as Enc
import numpy as np
from uncertainties import ufloat as uf

def Q_allunc(E, nATP, k_corr):
    """
    Get growth parameters for a methanogenesis quotient in which all the
    reagents have symettric uncertainties (ie they are all ufloats)
    """

    QE = deepcopy(E)#(name='En', pH=E.ocean_pH, CO2origin='pH', T=T, workoutID=False)

    ######### Nominal parameters ##########

    QE.composition['CO2(aq)'].activity = E.composition['CO2(aq)'].activity.n
    QE.composition['H2(aq)'].activity = E.composition['H2(aq)'].activity.n
    QE.composition['Methane(aq)'].activity = E.composition['Methane(aq)'].activity.n
    QE.composition['H2O(l)'].activity = E.composition['H2O(l)'].activity.n
    #quotient, gibbs, ATPgibbs, PowerSupply, max_k
    Q, G, GP, PS, mk = MGparams(QE, k_corr=k_corr, nATP=nATP)


    ######### min products, max reactants  #########
    QE.composition['CO2(aq)'].activity = (
      E.composition['CO2(aq)'].activity.n + E.composition['CO2(aq)'].activity.s)
    QE.composition['H2(aq)'].activity = (
      E.composition['H2(aq)'].activity.n + E.composition['H2(aq)'].activity.s)
    QE.composition['Methane(aq)'].activity = (
      E.composition['Methane(aq)'].activity.n - E.composition['Methane(aq)'].activity.s)
    QE.composition['H2O(l)'].activity = (
      E.composition['H2O(l)'].activity.n - E.composition['H2O(l)'].activity.s)

    Q_up, G_up, GP_up, PS_up, mk_up = MGparams(QE, k_corr=k_corr, nATP=nATP)

    ######### max products, min reactants  #########
    QE.composition['CO2(aq)'].activity = (
      E.composition['CO2(aq)'].activity.n - E.composition['CO2(aq)'].activity.s)
    QE.composition['H2(aq)'].activity = (
      E.composition['H2(aq)'].activity.n - E.composition['H2(aq)'].activity.s)
    QE.composition['Methane(aq)'].activity = (
      E.composition['Methane(aq)'].activity.n + E.composition['Methane(aq)'].activity.s)
    QE.composition['H2O(l)'].activity = (
      E.composition['H2O(l)'].activity.n + E.composition['H2O(l)'].activity.s)

    Q_do, G_do, GP_do, PS_do, mk_do = MGparams(QE, k_corr=k_corr, nATP=nATP)

    return {'Gibbs_Methanogenesis' : np.array([[G, G_do, G_up]]),
      'Quotient' : np.array([[Q, Q_do, Q_up]]),
      'PowerSupply' : np.array([[PS, PS_do, PS_up]]),
      'ATPGibbs': np.array([[GP, GP_do, GP_up]]),
      'max_k' : np.array([[mk, mk_do, mk_up]])}



def Q_salty_endmember(E, nATP, k_corr):
    """
    Get growth parameters for a methanogenesis quotient considering the widest
    endmemebers possible. So, variation across the salt content of the ocean
    affecting CO2. The 'up' bounds are for a salty ocean with maximised reactants
    and minimised products.
    The 'low' bounds are for a less salty ocean with maximised products and
    minimized reactants.
    """

    QE = deepcopy(E)#nc(name='En', pH=ocean_pH, CO2origin='pH', T=T, workoutID=False)

    CO2salt, H2Osalt = Enc.get_CO2_from_HTHeating(QE.env.T, QE.ocean_pH, salts=True)[:2]
    CO2salt273, H2Osalt273 = Enc.get_CO2_from_HTHeating(273.15, QE.ocean_pH, salts=True)[:2]

    ######### Nominal parameters ##########

    QE.composition['CO2(aq)'].activity = CO2salt[0]
    QE.composition['H2(aq)'].activity = E.calc_mol_H2(CO2salt273[0]).n
    QE.composition['Methane(aq)'].activity = E.calc_mol_CH4(CO2salt273[0]).n
    QE.composition['H2O(l)'].activity = H2Osalt[0]
    # quotient, gibbs, ATPgibbs, PowerSupply, max_k
    Q, G, GP, PS, mk = MGparams(QE, k_corr=k_corr, nATP=nATP)

    ######### min products, max reactants  #########
    QE.composition['CO2(aq)'].activity = (CO2salt[1])
    H2 = E.calc_mol_H2(CO2salt273[1])
    QE.composition['H2(aq)'].activity = (H2.n + H2.s)
    CH4 = E.calc_mol_CH4(CO2salt273[1])
    QE.composition['Methane(aq)'].activity = (CH4.n - CH4.s)
    QE.composition['H2O(l)'].activity = (H2Osalt[1])

    Q_up, G_up, GP_up, PS_up, mk_up = MGparams(QE, k_corr=k_corr, nATP=nATP)

    ######### max products, min reactants  #########
    QE.composition['CO2(aq)'].activity = (CO2salt[2])
    H2 = E.calc_mol_H2(CO2salt273[2])
    QE.composition['H2(aq)'].activity = (H2.n - H2.s)
    CH4 = E.calc_mol_CH4(CO2salt273[2])
    QE.composition['Methane(aq)'].activity = (CH4.n + CH4.s)
    QE.composition['H2O(l)'].activity = (H2Osalt[2])

    Q_do, G_do, GP_do, PS_do, mk_do = MGparams(QE, k_corr=k_corr, nATP=nATP)

    return {'Gibbs_Methanogenesis' : np.array([[G, G_do, G_up]]),
      'Quotient' : np.array([[Q, Q_do, Q_up]]),
      'PowerSupply' : np.array([[PS, PS_do, PS_up]]),
      'ATPGibbs': np.array([[GP, GP_do, GP_up]]),
      'max_k' : np.array([[mk, mk_do, mk_up]])}


def get_salty_Enc(E, salttype, CO2unc=0.0, H2Ounc=0.0):
    """Get a salty, nominal or less-salty enceladus with the activities as ufloats.
    """
    QE = deepcopy(E)#nc(name='En', pH=ocean_pH, CO2origin='pH', T=T, workoutID=False)

    CO2salt, H2Osalt = Enc.get_CO2_from_HTHeating(QE.env.T, QE.ocean_pH, salts=True)[:2]
    CO2salt273, H2Osalt273 = Enc.get_CO2_from_HTHeating(273.15, QE.ocean_pH, salts=True)[:2]
    aCO2, aH20, aH2, aCH4 = 0,0,0,0
    if salttype == 'nom':
        aCO2 = CO2salt[0]
        aH2O = H2Osalt[0]
        aH2 = QE.calc_mol_H2(CO2salt273[0])
        aCH4 = QE.calc_mol_CH4(CO2salt273[0])
    elif salttype =='high':
        aCO2 = CO2salt[1]
        aH2O = H2Osalt[1]
        aH2 = QE.calc_mol_H2(CO2salt273[1])
        aCH4 = QE.calc_mol_CH4(CO2salt273[1])
    elif salttype == 'low':
        aCO2 = CO2salt[2]
        aH2O = H2Osalt[2]
        aH2 = QE.calc_mol_H2(CO2salt273[2])
        aCH4 = QE.calc_mol_CH4(CO2salt273[2])
    else:
        raise ValueError('Unknown salt type sent to Q_salty')

    QE.composition['CO2(aq)'].activity = uf(aCO2, aCO2*CO2unc)
    QE.composition['H2(aq)'].activity = aH2
    QE.composition['Methane(aq)'].activity = aCH4
    QE.composition['H2O(l)'].activity = uf(aH2O, aH2O*H2Ounc)

    return QE

def Q_salty(E, nATP, k_corr, salttype, CO2unc=0.0, H2Ounc=0.0):
    """
    Get growth parameters for a methanogenesis quotient considering an endmember
    of salt content. given by salttype.

    salttype = {'nom', 'high', 'low'}
    """

    _QE = get_salty_Enc(E, salttype, CO2unc=CO2unc, H2Ounc=H2Ounc)
    QE = deepcopy(_QE)
    aCO2 = _QE.composition['CO2(aq)'].activity
    aH2O = _QE.composition['H2O(l)'].activity
    aH2 = _QE.composition['H2(aq)'].activity
    aCH4 = _QE.composition['Methane(aq)'].activity

    ######### Nominal parameters ##########

    QE.composition['CO2(aq)'].activity = aCO2.n
    QE.composition['H2(aq)'].activity = aH2.n
    QE.composition['Methane(aq)'].activity = aCH4.n
    QE.composition['H2O(l)'].activity = aH2O.n
    # quotient, gibbs, ATPgibbs, PowerSupply, max_k
    Q, G, GP, PS, mk = MGparams(QE, k_corr=k_corr, nATP=nATP)

    ######### min products, max reactants  #########
    QE.composition['CO2(aq)'].activity = aCO2.n+aCO2.s
    QE.composition['H2(aq)'].activity = (aH2.n + aH2.s)
    QE.composition['Methane(aq)'].activity = (aCH4.n - aCH4.s)
    QE.composition['H2O(l)'].activity = aH2O.n-aH2O.s

    Q_up, G_up, GP_up, PS_up, mk_up = MGparams(QE, k_corr=k_corr, nATP=nATP)

    ######### max products, min reactants  #########
    QE.composition['CO2(aq)'].activity = aCO2.n-aCO2.s
    QE.composition['H2(aq)'].activity = (aH2.n - aH2.s)
    QE.composition['Methane(aq)'].activity = (aCH4.n + aCH4.s)
    QE.composition['H2O(l)'].activity = aH2O.n+aH2O.s

    Q_do, G_do, GP_do, PS_do, mk_do = MGparams(QE, k_corr=k_corr, nATP=nATP)

    return {'Gibbs_Methanogenesis' : np.array([[G, G_do, G_up]]),
      'Quotient' : np.array([[Q, Q_do, Q_up]]),
      'PowerSupply' : np.array([[PS, PS_do, PS_up]]),
      'ATPGibbs': np.array([[GP, GP_do, GP_up]]),
      'max_k' : np.array([[mk, mk_do, mk_up]])}
