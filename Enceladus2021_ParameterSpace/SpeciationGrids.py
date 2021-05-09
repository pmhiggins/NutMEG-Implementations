import numpy as np
import pandas as pd
from scipy import interpolate

import matplotlib.pyplot as plt

"""
This is for extracting information from the carbonate speciation data files.
"""

def extract_save():
    pHnames = ['7.0','7.5','8.0','8.5','9.0','9.5','10.0','10.5','11.0','11.5','12.0']

    pHfloats = np.linspace(7.,12., num=11)
    Tfloats = np.linspace(273.15, 473.15, num=21)

    aH2O = np.ndarray((len(pHfloats), len(Tfloats)))
    aCO2 = np.ndarray((len(pHfloats), len(Tfloats)))
    pHHT = np.ndarray((len(pHfloats), len(Tfloats)))

    aH2O_hs = np.ndarray((len(pHfloats), len(Tfloats)))
    aCO2_hs = np.ndarray((len(pHfloats), len(Tfloats)))
    pHHT_hs = np.ndarray((len(pHfloats), len(Tfloats)))
    aCO2_hs_factor = np.ndarray((len(pHfloats), len(Tfloats)))

    aH2O_ls = np.ndarray((len(pHfloats), len(Tfloats)))
    aCO2_ls = np.ndarray((len(pHfloats), len(Tfloats)))
    pHHT_ls = np.ndarray((len(pHfloats), len(Tfloats)))
    aCO2_ls_factor = np.ndarray((len(pHfloats), len(Tfloats)))


    fn_preamble = 'E21data/Speciation/'
    for i, pH in enumerate(pHnames):
        df=pd.read_csv(fn_preamble+'nominalCO2/pH'+pH+'.csv', sep=',')
        aH2O[i] = df['a_H2O']
        aCO2[i] = df['a_CO2']
        pHHT[i] = df['pH']
        df2=pd.read_csv(fn_preamble+'highsalt/pH'+pH+'.csv', sep=',')
        aH2O_hs[i] = df2['a_H2O']
        aCO2_hs[i] = df2['a_CO2']
        pHHT_hs[i] = df2['pH']
        aCO2_hs_factor[i] = df2['a_CO2factor']
        df3=pd.read_csv(fn_preamble+'lowsalt/pH'+pH+'.csv', sep=',')
        aH2O_ls[i] = df3['a_H2O']
        aCO2_ls[i] = df3['a_CO2']
        pHHT_ls[i] = df3['pH']
        aCO2_ls_factor[i] = df3['a_CO2factor']

    Tgrid, pHgrid = np.meshgrid(Tfloats, pHfloats)

    np.save(fn_preamble+'Tgrid.npy', Tgrid)
    np.save(fn_preamble+'pHgrid.npy', pHgrid)

    #nominals
    np.save(fn_preamble+'nominalCO2/np/aH2O.npy', aH2O)
    np.save(fn_preamble+'nominalCO2/np/aCO2.npy', aCO2)
    np.save(fn_preamble+'nominalCO2/np/pHHT.npy', pHHT)

    # high salt
    np.save(fn_preamble+'highsalt/np/aH2O.npy', aH2O_hs)
    np.save(fn_preamble+'highsalt/np/aCO2.npy', aCO2_hs)
    np.save(fn_preamble+'highsalt/np/pHHT.npy', pHHT_hs)
    np.save(fn_preamble+'highsalt/np/aCO2factor.npy', aCO2_hs_factor)

    # low salt
    np.save(fn_preamble+'lowsalt/np/aH2O.npy', aH2O_ls)
    np.save(fn_preamble+'lowsalt/np/aCO2.npy', aCO2_ls)
    np.save(fn_preamble+'lowsalt/np/pHHT.npy', pHHT_ls)
    np.save(fn_preamble+'lowsalt/np/aCO2factor.npy', aCO2_ls_factor)

# extract_save()


def interpo_CO2_H2O(Temp=273.15, oceanpH=8.0, salt='nominalCO2'):
    """
    Interpolate over the data files to estimate water and CO2 activity at a
    given Temp and ocean_pH
    """
    fn_preamble = 'E21data/Speciation/'

    aH2O = np.load(fn_preamble+salt+'/np/aH2O.npy')
    aCO2 = np.load(fn_preamble+salt+'/np/aCO2.npy')
    # pHHT = np.load('HTHeatingdata/nominalCO2/pHHT.npy')

    pHfloats = np.linspace(7.,12., num=11)
    Tfloats = np.linspace(273.15, 473.15, num=21)

    fCO2 = interpolate.interp2d(Tfloats,pHfloats,aCO2,kind='cubic')

    fH2O = interpolate.interp2d(Tfloats,pHfloats,aH2O,kind='cubic')

    return fCO2(Temp, oceanpH)[0], fH2O(Temp, oceanpH)[0]

# print(interpo_CO2_H2O(Temp=273, oceanpH=8))
# print(interpo_CO2_H2O(Temp=273, oceanpH=8, salt='highsalt'))
# print(interpo_CO2_H2O(Temp=273, oceanpH=8, salt='lowsalt'))
# print(interpo_CO2_H2O(Temp=273, oceanpH=8.25, salt='lowsalt'))
# print(interpo_CO2_H2O(Temp=273, oceanpH=8.5, salt='lowsalt'))
