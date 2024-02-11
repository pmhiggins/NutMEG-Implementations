import math, os, sys
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt


def ens_point5(x):
    """ Round a value to the nearest 0.5 """
    return round(x*2)/2.


def get_saltname_as_float(salt):
    """
    Convert the string identifiers from Higgins et al 2021 into floats of
    salt concentration
    """
    _sn = {'lowsalt':0.05, 'nominalsalt':0.1, 'highsalt':0.2}

    return _sn[salt]


def get_setup_spec_params(Temp=273.15, oceanpH=8.0, salt=0.1):
    """
    Get the input parameters to begin the speciation, based on an interpolation
    across the speciation from Higgins et al. 2021.

    If a string is passed for the salt parameter, return the exact value from
    Higgins et al. 2021.

    returns: molality of CO2, an activity for water, pH, molality of DIC,
      molality for Cl., molality for Na.

    Na and Cl molalities are constant here, but the speciation solver is
    open to Na+ so that should not be an issue
    """

    if type(salt) == type(''):
        # string has been passed, so use the older nominal/low/high terminology
        salt = get_saltname_as_float(salt)

    # note this allows [Cl] up to 0.4 M, the original
    # distribution is up to 0.2, the first 21 values of this one.
    Cl = np.logspace(math.log10(0.05),math.log10(0.4), num=31)
    # Carbs = np.logspace(math.log10(0.01),math.log10(0.1), num=21)
    # pHs = np.linspace(7,12,num=11)
    pHs = np.linspace(7,12,num=51)
    dom_ion = ['HCO3-', 'HCO3-','HCO3-','HCO3-','HCO3-','HCO3-','HCO3-','HCO3-','CO3-2','CO3-2','CO3-2']

    # find the closest element in the list of Cls from H+2021
    Cl = np.asarray(Cl)
    Cl_id = (np.abs(Cl - salt)).argmin()
    _salt = Cl[Cl_id]

    # same again for the pH
    # find the closest element in the list of Cls from H+2021
    pHs = np.asarray(pHs)
    pH_id = (np.abs(pHs - oceanpH)).argmin()
    _pH = pHs[pH_id]

    # although this is reading in 'high salt' fits, the first 21 entries
    # corr. to [Cl] 0.05 to 0.2 are IDENTICAL to the original distribution (H+ 2021).
    # only entries 22-31 are extrapolated.
    _DIC = np.load('DIC_Cl_correlation/DIC_tots_HS.npy').T
    _mCO2 = np.load('DIC_Cl_correlation/CO2_tots_HS.npy').T

    # for [Cl] that are not on the grid, use a polyfit
    polyDIC = np.polyfit(np.log10(Cl), np.log10(_DIC[pH_id]), deg=1)
    mDIC = 10**(polyDIC[1] + (polyDIC[0]*np.log10(salt)))

    polyCO2 = np.polyfit(np.log10(Cl), np.log10(_mCO2[pH_id]), deg=1)
    mCO2 = 10**(polyCO2[1] + (polyCO2[0]*np.log10(salt)))


    return mCO2, 1., _pH, mDIC, _salt, _salt
