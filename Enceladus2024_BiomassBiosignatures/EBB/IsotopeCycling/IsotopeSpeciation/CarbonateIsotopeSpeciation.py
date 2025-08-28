import sys, os, math

sys.path.append(os.path.dirname(__file__)+'/../../../../../NutMEG')
import NutMEG as nm
from NutMEG.reactor.saved_systems.Enceladus import Enceladus
from NutMEG.reaction.reagent import reagent as rgt

sys.path.append(os.path.dirname(__file__)+'/../../')
from EncBmBs_utils import IsotopeConversions as ICtools
from EncBmBs_utils import EnceladusFetcher

sys.path.append(os.path.dirname(__file__)+'../')
from FracFactors import CarbonateFracFactorDeines1974

from copy import deepcopy
import sympy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# import warnings
# warnings.filterwarnings('ignore')

class CarbonateIsotopeSpeciation:

    def __init__(self, rtr_T, rtr_T0, R_CO2_T0, CFracFactor='Deines'):
        self.rtr_T = rtr_T
        self.rtr_T0 = rtr_T0

        if CFracFactor=='Deines':
            self.FF = CarbonateFracFactorDeines1974()
        else:
            self.FF = CFracFactor

        self.R_CO2_T0 = R_CO2_T0
        self.set_Rs_from_R_CO2() # set R vals for DIC, CO3 and HCO3 at T0
        self.set_R_CO2_T() # set R vals for DIC, CO3 and HCO3 at T

    @staticmethod
    def rtr_skip(Ts, concs_T0, concs_T):
        """
        If you prefer to just pass known compositions at known temperatures
        instead of NutMEG reactor objects, this function will generate the
        reactors needed at the temperatures of interest.
        """
        rtr_T0 = Enceladus('rtr_T0', T=Ts[0], nominals=True)
        rtr_T0.composition['CO2(aq)'].molal = concs_T0['CO2(aq)']
        rtr_T0.add_reagent(rgt('HCO3-', rtr_T0.env, molal = concs_T0['HCO3-']))
        rtr_T0.add_reagent(rgt('CO3-2', rtr_T0.env, molal = concs_T0['CO3-2']))

        rtr_T = Enceladus('rtr_T', T=Ts[1], nominals=True)
        rtr_T.composition['CO2(aq)'].molal = concs_T['CO2(aq)']
        rtr_T.add_reagent(rgt('HCO3-', rtr_T.env, molal = concs_T['HCO3-']))
        rtr_T.add_reagent(rgt('CO3-2', rtr_T.env, molal = concs_T['CO3-2']))

        return rtr_T0, rtr_T


    def set_Rs_from_R_CO2(self):
        """
        Assuming equilibrium, use the R value of CO2 to compute the R value
        of the other dissolved carbonate species
        """

        a_CO3, a_HCO3 = self.get_fracfactors(self.rtr_T0)
        mCO2, mHCO3, mCO3, DIC = self.get_C8_molalities(self.rtr_T0)

        self.R_CO3_T0 = a_CO3 * self.R_CO2_T0
        self.R_HCO3_T0 = a_HCO3 * self.R_CO2_T0

        ans = mCO2*ICtools.RtoF(self.R_CO2_T0)
        ans += (mHCO3*ICtools.RtoF(self.R_HCO3_T0))
        ans += (mCO3* ICtools.RtoF(self.R_CO3_T0))
        ans /= DIC

        self.R_DIC = ICtools.FtoR(ans)


    def set_R_CO2_T(self):
        """
        Solve the carbonate isotope speciation at the elevated temperature
        """

        a_CO3, a_HCO3 = self.get_fracfactors(self.rtr_T)
        mCO2, mHCO3, mCO3, DIC = self.get_C8_molalities(self.rtr_T)

        x = sp.Symbol('x') # let x === R_CO2, we need to solve for this

        f = ((mCO2*x/(1+x))+(mHCO3*a_HCO3*x/(1+(a_HCO3*x)))+(mCO3*a_CO3*x/(1+(a_CO3*x)))) - (DIC*ICtools.RtoF(self.R_DIC))


        s = sp.solve(f) # get the solutions for R_CO2
        s2 = np.array(s, dtype=complex) # convert to a numpy array as complex no.s
        realposvals = s2[s2.real > 0.] # choose solutions where the real part is +ve

        if len(realposvals) >1:
            # if more than one root has positive real part,
            # choos the one with the smallest imaginary part
            realposvals = realposvals[np.where(np.abs(realposvals.imag) == np.amin(np.abs(realposvals.imag)))[0]]

        self.R_CO2_T = realposvals[0].real
        self.R_CO3_T = self.R_CO2_T * a_CO3
        self.R_HCO3_T = self.R_CO2_T * a_HCO3



    @staticmethod
    def get_C8_molalities(_rtr):
        """ Return molalities of the three carbonate species and total DIC"""

        mCO2 = _rtr.composition['CO2(aq)'].molal
        mHCO3 = _rtr.composition['HCO3-'].molal
        mCO3 = _rtr.composition['CO3-2'].molal
        DIC = mCO2 + mHCO3 + mCO3

        return mCO2, mHCO3, mCO3, DIC

    def get_fracfactors(self, _rtr):
        """
        Return fractionation factors between CO2 and CO3, and CO2 and HCO3
        at the reactor temperature.
        """

        a_CO3 = self.FF.CO3_wrt_CO2(_rtr.env.T)
        a_HCO3 = self.FF.HCO3_wrt_CO2(_rtr.env.T)

        return a_CO3, a_HCO3



def CIStest():
    """
    test implementation, compute and print out the
    carbonate speicaiton delta C values at 0C and 100C for a pH8 Enceladus
    """
    testE = EnceladusFetcher().fetch(373.15, 8.0, np.float64(0.1),
      model='pitzerPHREEQCnoGases', fixGases=False, spec_dr='spec_T_273-473_pH_8-9')
    testE0 = EnceladusFetcher().fetch(273.15, 8.0, np.float64(0.1),
      model='pitzerPHREEQCnoGases', fixGases=False, spec_dr='spec_T_273-473_pH_8-9')

    print(testE.env)
    print(testE0.env)

    testFF = CarbonateFracFactorDeines1974()

    R = ICtools.dCtoR(72)

    CIS = CarbonateIsotopeSpeciation(testE, testE0, R, testFF)
    print('T = '+str(CIS.rtr_T0.env.T)+' K')

    print('CO2', ICtools.RtodC(CIS.R_CO2_T0))
    print('CO3', ICtools.RtodC(CIS.R_CO3_T0))
    print('HCO3', ICtools.RtodC(CIS.R_HCO3_T0))
    print('DIC', ICtools.RtodC(CIS.R_DIC))

    print('T = '+str(CIS.rtr_T.env.T)+' K')
    print('CO2', ICtools.RtodC(CIS.R_CO2_T))
    print('CO3', ICtools.RtodC(CIS.R_CO3_T))
    print('HCO3', ICtools.RtodC(CIS.R_HCO3_T))
    print('DIC', ICtools.RtodC(CIS.R_DIC))


# CIStest()


def CIS_plot(stdout=False):
    Clconcs = np.logspace(math.log10(0.05),math.log10(0.2), num=21) #21
    Clconcs = np.round(Clconcs, decimals=5)
    _top_R_CO2s = np.linspace(1/(84.+13), 1/(84.-13), num=11) #11

    pHvals = [8.0,9.0,10.0,11.0] #.5 pHs are optional too
    Tvals = np.linspace(273.15, 393.15, num=13) #13
    spec_dr = 'spec_T_273-473_pH_7-12'



    DICfig, DICaxs = plt.subplots(nrows=2, ncols=2, figsize=(8,8))
    DICaxs = DICaxs.flatten()
    cmap = plt.get_cmap("viridis", len(Clconcs)+2)
    for pH, ax in zip(pHvals, DICaxs):
        for i, Cl in enumerate(Clconcs):

            this_df = pd.read_csv('../../data/speciation/spec_T_273-473_pH_7-12/Clconc_'+str(Cl)+'/spec_1bar_pitzerPHREEQCnoGases.csv')
            _df = this_df[this_df['pH_bo'] == pH]

            testE0 = EnceladusFetcher().fetch(Tvals[0], pH, Cl,
              model='pitzerPHREEQCnoGases', fixGases=False, spec_dr='spec_T_273-473_pH_7-12')
            testE = deepcopy(testE0)


            DICs = []
            CO2s = []
            HCO3s = []
            CO3s = []
            FF = CarbonateFracFactorDeines1974()
            R0 = ICtools.dCtoR(60)

            for T in Tvals:
                if stdout:
                    print('fetching pH: ',pH,'; Cl: ',Cl,' temperature: ',T)
                __df = _df[_df['T'] == T]
                testE.composition['CO2(aq)'].molal = __df['mCO2(aq)'].tolist()[0]
                testE.composition['HCO3-'].molal = __df['mHCO3-'].tolist()[0]
                testE.composition['CO3-2'].molal = __df['mCO3-2'].tolist()[0]
                testE.env.T = T

                CIS = CarbonateIsotopeSpeciation(testE, testE0, R0, FF)
                CO2s.append(ICtools.RtodC(CIS.R_CO2_T))
                CO3s.append(ICtools.RtodC(CIS.R_CO3_T))
                HCO3s.append(ICtools.RtodC(CIS.R_HCO3_T))
                DICs.append(ICtools.RtodC(CIS.R_DIC))

            ax.plot(Tvals-273.15, DICs, ls='dotted', c=cmap(i))
            ax.plot(Tvals-273.15, CO2s, c=cmap(i))
            ax.plot(Tvals-273.15, HCO3s, ls='dashed', c=cmap(i))
            ax.plot(Tvals-273.15, CO3s, ls='dashdot', c=cmap(i))

        ax.set_xlabel(r'Temperature [$\degree$C]')
        ax.set_ylabel(r'$\delta^{13}$C [$\perthousand$]')
        ax.text(110, 60.5, r'pH at 0$\degree$C: '+str(pH), va='bottom', ha='right')
        ax.set_xlim(0,120)

    DICaxs[-1].plot(np.nan, np.nan, ls='dotted', c='k', label='$\delta^{13}$C$_{\mathregular{DIC}}$')
    DICaxs[-1].plot(np.nan, np.nan, c='k', label='$\delta^{13}$C$_{\mathregular{CO2}}$')
    DICaxs[-1].plot(np.nan, np.nan, ls='dashed', c='k', label='$\delta^{13}$C$_{\mathregular{HCO3}}$')
    DICaxs[-1].plot(np.nan, np.nan, ls='dashdot', c='k', label='$\delta^{13}$C$_{\mathregular{CO3}}$')

    DICfig.subplots_adjust(top=0.982, bottom=0.2, left=0.091, right=0.978)

    DICaxs[-1].legend(loc='lower center',ncol=2, handlelength=4, bbox_to_anchor=[0.5, -0.43])


    cb_ax = DICfig.add_axes([0.091,0.08,0.4,0.03])
    norm = mpl.colors.Normalize(vmin=Clconcs[0], vmax=Clconcs[-1])

    DICfig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cb_ax, orientation='horizontal', label=r'Ocean Cl$^-$ molality [mol kg$^{-1}$]')

    plt.savefig('DICvariation.pdf')
    plt.savefig('DICvariation.png', dpi=300)

    plt.close()

# CIS_plot(stdout=True)
