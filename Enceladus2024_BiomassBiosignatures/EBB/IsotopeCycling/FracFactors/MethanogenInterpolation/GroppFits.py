import sys, os; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__)+'/../../../../../../NutMEG')


import pandas as pd
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import matplotlib.tri as tri

import NutMEG
sys.path.append(os.path.dirname(__file__)+'/../../../')
from EncBmBs_utils import EnceladusFetcher

def get_stdG(T=333):
    """
    Calculate the standard gibbs free energy of methanogenesis at given temperature
    """
    en = NutMEG.environment()
    en.T = T
    CO2aq = NutMEG.reaction.reagent('CO2(aq)', en, phase='aq')
    H2aq = NutMEG.reaction.reagent('H2(aq)', en, phase='aq')
    CH4aq = NutMEG.reaction.reagent('Methane(aq)', en, phase='aq')
    H2O = NutMEG.reaction.reagent('H2O(aq)', en, phase='l')

    thermalMG = NutMEG.reaction.reaction(
      {CO2aq:1, H2aq:4}, {CH4aq:1, H2O:2}, en)

    thermalMG.update_std_molar_gibbs_G()
    stdG = thermalMG.std_molar_gibbs
    return stdG


def conv_DeltaG_to_H2(Aff, stdG, CO2, T=333, CH4=10e-6):
    """
    Convert molar delta G to [H2], given stdG and [CO2]
    assume [CH4] is 10e-6 unless told otherwise (fixed to Gropp et al. value)
    """
    DeltaG = -Aff
    R=8.31
    lnH2 = 0.25*(((stdG-DeltaG)/(R*T))+((np.log(CH4) -np.log(CO2))))
    return np.exp(lnH2)


def getEncConcs(pH, saltlvl, H2corr, T=333.15, CH4corr=0, dr='spec_T_273-473_pH_8-9'):
    """
    Returns the activity of H2 and CO2 in an Enceladus habitat using the pH,
    saltlvl, H2corr, CH4corr and T indicators, saved in directory dr.
    """

    Enc = EnceladusFetcher().fetch(
      T,
      pH,
      np.float64(saltlvl),
      P = 1, model='pitzerPHREEQCnoGases',
      spec_dr =  dr,
      fixGases=False,
      gases= {'H2':H2corr, 'CH4':CH4corr},
      GasScale = 1.0)

    return [np.log10(Enc.composition['H2(aq)'].activity), np.log10(Enc.composition['CO2(aq)'].activity)]

def __main__(overlay=False):
    """
    Read in Gropp 2022 data, collect epsilon as function of [CO2] and [H2],
    then interpolate.

    Overplots Enceladus ranges of [CO2] and [H2]

    Save the figure.
    """
    df = pd.read_csv(os.path.dirname(__file__)+'/../../../data/MethanogenFractionation/Gropp2023_epsilons.csv')

    CO2s = ['0.1mM', '1mM', '3mM', '10mM', '100mM']
    CO2floats = [-4,-3,np.log10(0.003),-2,-1] # in log10 space M
    NaNs = [-2,0, -3,-2,-2] # number of NaNs on the end of each data fit

    #  X1, X2, Y1, Y2 (for dG vs E, [H2] vs, E)
    Gtot, Htot, EHtot, EGtot, Ctot = [],[],[],[],[]

    stdG = get_stdG() # standard free energy of methanogenesis at 60 C

    for CO2, CO2f, _NaN in zip(CO2s, CO2floats, NaNs):
        # Gropp+2022 has a line, and we have two columns, for each [CO2]

        G = df[CO2+' DeltaG']
        Y = df[CO2+' e']
        H = np.log10(conv_DeltaG_to_H2(1000*G, stdG, 10**(CO2f)))


        spl=None
        # k param gives order to fit the spline. 2 works best across the board
        if _NaN == 0:
            Hspl = interp.InterpolatedUnivariateSpline(H,Y, k=2)
            Gspl = interp.InterpolatedUnivariateSpline(G,Y, k=2)

        else:
            Hspl = interp.InterpolatedUnivariateSpline(H[:_NaN],Y[:_NaN], k=2)
            Gspl = interp.InterpolatedUnivariateSpline(G[:_NaN],Y[:_NaN], k=2)


        fitH = np.linspace(-9, -1, num=1000)
        fitG = np.linspace(0, 200, num=1000)

        fitGY = Gspl(fitG)
        fitHY = Hspl(fitH)


        # fits are no good less than or more than the end of the data range
        # they diverge significantly
        # so cap off tops and tails to the be value at the end of the range.
        for i,fG in enumerate(fitG):
            if fG < min(G):
                fitGY[i] = Y[0]
            elif fG > max(G):
                fitGY[i] = Y.tolist()[_NaN-1]
        for i,fH in enumerate(fitH):
            if fH < min(H):
                fitHY[i] = Y[0]
            elif fH > max(H):
                fitHY[i] = Y.tolist()[_NaN-1]

        Htot.extend([a for a in fitH])
        Gtot.extend([a for a in fitG])

        EHtot.extend([a for a in fitHY])
        EGtot.extend([a for a in fitGY])

        Ctot.extend([CO2f for a in fitG])



    Atot = [1 + (e/1000) for e in EHtot]


    # figure showing interpolated contour,
    # and overlaying Enceladus ranges of [CO2] and [H2]

    figC, axC = plt.subplots(nrows=1, ncols=1)
    # range of alpha to consider on colormap
    levels = np.linspace(1.0, 1.1, 101)

    f = interp.CloughTocher2DInterpolator(list(zip(Htot, Ctot)), EHtot)
    # f = interp.RectBivariateSpline(Htot, Ctot, EHtot)

    # np.save(os.path.dirname(__file__)+'/../../../data/MethanogenFractionation/GroppInterp.npy', f, allow_pickle=True)
    # f = np.load(os.path.dirname(__file__)+'/../../../data/MethanogenFractionation/GroppInterp.npy', allow_pickle=True)[()]

    x2 = np.linspace(-9, -1, num=1000)
    y2 = np.linspace(-4, -1, num=1000)
    X2, Y2 = np.meshgrid(x2, y2)
    Afit = -f(X2, Y2)


    cbm = axC.pcolormesh(X2, Y2, Afit, cmap='viridis_r', vmin=-100, vmax=0)


    figC.colorbar(cbm, label=r'Biotic methanogenesis enrichment factor $\epsilon^{bio, KIE}_{CO2/CH4} [\perthousand]$')
    axC.tricontour(Htot, Ctot, EHtot, levels=[1.061], colors=['k'])

    axC.set_ylabel(r'$\log_{10}$[CO$_{2}$] [M]')
    axC.set_xlabel(r'$\log_{10}$[H$_{2}$] [M]')
    axC.set_title('Methanogen CO$_{2}$-CH$_{4}$ C$^{13}$ enrichment at 60$\degree$C, [CH$_4$] = 10$^{-5}$ M'+'\n'+r'$\quad$')

    if overlay:
        # overlay enceladus habitats
        pH8_lowsalt_lowH2 = getEncConcs(8, 0.05, -1)
        pH8_lowsalt_highH2 = getEncConcs(8, 0.05, 1)
        pH8_highsalt_lowH2 = getEncConcs(8, 0.2, -1)
        pH8_highsalt_highH2 = getEncConcs(8, 0.2, 1)

        pH8box = [pH8_lowsalt_lowH2,
          pH8_lowsalt_highH2,
          pH8_highsalt_highH2,
          pH8_highsalt_lowH2,
          pH8_lowsalt_lowH2]


        pH9_lowsalt_lowH2 = getEncConcs(9, 0.05, -1)
        pH9_lowsalt_highH2 = getEncConcs(9, 0.05, 1)
        pH9_highsalt_lowH2 = getEncConcs(9, 0.2, -1)
        pH9_highsalt_highH2 = getEncConcs(9, 0.2, 1)

        pH9box = [pH9_lowsalt_lowH2,
          pH9_lowsalt_highH2,
          pH9_highsalt_highH2,
          pH9_highsalt_lowH2,
          pH9_lowsalt_lowH2]

        # transform the lists of coordinates into X and Y lists
        pH9xs, pH9ys = zip(*pH9box)
        pH8xs, pH8ys = zip(*pH8box)

        axC.plot(pH8xs, pH8ys, c='pink')
        axC.plot(pH9xs, pH9ys, c='tab:orange')

        axC.text(-3,-2.25, 'pH: 8', c='pink', ha='left', va='bottom')
        axC.text(-4.25,-2.75, 'pH: 9', c='tab:orange', ha='left', va='bottom')



    figC.savefig('GroppContour_EncOverlay.png', dpi=600)
    # plt.show()

__main__(overlay=True)
