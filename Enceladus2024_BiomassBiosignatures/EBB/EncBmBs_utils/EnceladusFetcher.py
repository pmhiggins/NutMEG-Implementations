from EncBmBs_utils import DataFrameFetcher
# import DataFrameFetcher
import NutMEG as nm
from NutMEG.reactor.saved_systems.Enceladus import Enceladus as Enc

from copy import deepcopy
import numpy as np
import random

class EnceladusFetcher:
    """
    Class for generating an NutMEG.saved_organism.Enceladus object, and
    adjusting its properties for analysis.
    """

    def __init__(self):
        # this dict translates Enc reagents to their keys in the speciation dataframe
        # we don't need to put in everything from the speciation.
        self.species_amendments = {
          'H2O(l)' : 'H2O',
          'H+' : 'H+',
          'OH-' : 'OH-',
          'Na+' : 'Na+',
          'Cl-' : 'Cl-',
          'HCO3-' : 'HCO3-',
          'CO2(aq)' : 'CO2(aq)',
          'CO3--' : 'CO3-2'
          }

    def add_species_amendments(self, extra_species):
        self.species_amendments.update(extra_species)


    def update_reagents(self, rtr, df, extra_species={}):
        """ update the composition of a passed reactor-like object and return it

        df : pandas DataFrame
            DataFrame generated from an EncSpecRetriever that contains reaktoro
            chemical speciation output. This will be used to update rtr's
            composition attribute.
        extra_species : dict
            Any additional chemical species to include. By default only
            self.species_amendments will be updated based on df. If the species
            is not currently in the rtr it will be added based on entries in df.
        """
        _specs = deepcopy(self.species_amendments)
        _water = deepcopy(self.species_amendments)
        _df = deepcopy(df)

        _specs.update(extra_species)

        for k,v in _specs.items():
            try:
                g = float(_df['g'+v])
                rtr.composition[k].gamma = _df['g'+v]

                if k == 'H2O(l)':
                    a = float(_df['a'+v])
                    rtr.composition[k].activity = a
                else:
                    m = float(_df['m'+v])
                    rtr.composition[k].activity = g*m

                    # set both conc and molal to the same value, with molal act. coeff.
                    # this is because some legacy code used activity=conc for reactions
                    # so any taking place with conc should still be valid
                    rtr.composition[k].conc = m
                    rtr.composition[k].molal = m
            except KeyError:
                # the species is not present in the reactor yet
                g = float(_df['g'+v])
                m = float(_df['m'+v])
                rct = nm.reaction.reagent(k, rtr.env, phase='aq',
                  conc=m, activity=m*g, molal=m, gamma=g)
                rtr.add_reagent(rct)
        rtr.DIC = rtr.composition['HCO3-'].molal + rtr.composition['CO2(aq)'].molal + rtr.composition['CO3--'].molal
        return rtr


    def update_gases(self, rtr, CO2_bo=None, H2021_saltlevel=None, how='uniform',
      GasScale=1.0):
        """
        Update the dissolved gas concentration (H2 and CH4) in reactor-like
        object rtr by factor GasScale.

        CO2_bo : float
            molality of CO2 in Enceladus' bulk ocean for this case study
        how : string or dict
            how to calculate the bulk ocean dissolved gas concentration.
            [H2] and [CH4] are calculated as ufloats, and this attribute
            describes how to assign a float value based on their mean and
            standard deviation values (here the std dev is used as a min/max
            endmember for simplicity).
            Pass a dict  with keys 'H2' and 'CH4' each with a
            float value between +1 and -1 to manually assign a position within
            the range of values.
            Pass 'max' to maximise energy yield (max H2 min CH4)
            Pass 'min' to minimize energy yield (min H2 max CH4)
            Pass 'nom' for the nominal case ("mean")
            Pass 'uniform' for an independent random uniform selection for both
            between nom - std and nom + std
        GasScale : float (default 1.0)
            scale by which to change the the gas concentrations after their
            bulk ocean values are assigned.
        """


        if CO2_bo == None:
            # if the CO2 at ocean top is unknown, use the nominal value at a
            # given Higgins+2021-like saltlevel
            # if that isn't known, use the Higgins+2021 'nominal' ocean.
            salt_index = 'nom'
            if type(H2021_saltlevel) == type(''):
                salt_index = {'nom':0, 'high':1, 'low':2}[saltlevel]
            mol_CO2lst_oc = Enc.get_CO2_from_HTHeating(T=273.15,
              pH_0=rtr.ocean_pH, nominals=False, CO2unc=0., salts=True)[0]

            CO2_bo= mol_CO2lst_oc[salt_index]

        mH2 = rtr.calc_mol_H2(CO2_bo)
        mCH4 = rtr.calc_mol_CH4(CO2_bo)
        if type(how) == type({}):
            mH2 = mH2.n + (how['H2']*mH2.s)
            mCH4 = mCH4.n + (how['CH4']*mCH4.s)
        elif how == 'max':
            mH2 = mH2.n + mH2.s
            mCH4 = mCH4.n - mCH4.s
        elif how == 'min':
            mH2 = mH2.n - mH2.s
            mCH4 = mCH4.n + mCH4.s
        elif how == 'nom':
            mH2 = mH2.n
            mCH4 = mCH4.n
        elif how == 'uniform':
            rH2 = random.uniform(-1.,1.)
            rCH4 = random.uniform(-1.,1.)
            mH2 = mH2.n + (rH2*mH2.s)
            mCH4 = mCH4.n + (rCH4*mCH4.s)
        else:
            raise ValueError('Unknown methodology for H2 and CH4 molalities!')

        rtr.composition['Methane(aq)'].activity = mCH4*GasScale
        rtr.composition['Methane(aq)'].molality = mCH4*GasScale
        rtr.composition['Methane(aq)'].conc = mCH4*GasScale

        rtr.composition['H2(aq)'].activity = mH2*GasScale
        rtr.composition['H2(aq)'].molality = mH2*GasScale
        rtr.composition['H2(aq)'].conc = mH2*GasScale

        return rtr



    def fetch(self, T, pH_bo, saltlevel, P=1, model='pitzerPHREEQCnoGases',
      fixGases=True, spec_dr='spec_T_273-473_pH_7-12', gases='nom',
      GasScale=1.0):
        """
        Generate and return an Enceladus object representing temperature T,
        bulk ocean pH pH, Cl- concentraiton saltlevel, at pressure P bar, and
        a speciaton based on `model` and `fixGases` leaded from `spec_dr`.
        If gas concentrations are to be set and/or scaled, pass `gases` and
        `GasScale` accordingly.

        gases : string or dict
            Pass as rkt to use whatever dissolved gas concs (if any) are in the
            chemical speciation. Otherwise, refer to update_gases for the options.
            The default here is to fix the gases to the nominal plume
            mixing ratio vs CO2. GasScale will not work if gases='rkt'.

        """

        DFF = DataFrameFetcher(Ts=[T], pHs=[pH_bo], dr=spec_dr, salts=[saltlevel],
          P=P, model=model, fixGases=fixGases)

        # get a speciation df
        df = None
        try:
            df = DFF.spec()
            if df.empty:
                raise FileNotFoundError
        except FileNotFoundError:
            emsg = 'T, pH, or salt level not found in speciation.'
            raise FileNotFoundError(emsg)

        #  build a standard Enceladus object per Higgins et al 2021
        # speciation details (e.g. salt level) are not important here, because we
        # are about to re-assign them bsed on the speciation.
        _saltlevel = saltlevel
        if type(saltlevel) is type(np.float64(0.)):
                _saltlevel = 'nom'
        _Enc = Enc('E1', pH=pH_bo, T=T, CO2origin='HTHeatingSalts',
          saltlevel=_saltlevel, nominals=True)

        _Enc = self.update_reagents(_Enc, df)


        if gases != 'rkt':
            # update dissolved gas concentrations.
            # first pull out the concentration of CO2 in the bulk ocean
            # for use in getting bulk ocean dissolved gas concentrations
            DFF.Ts = [273.15]
            df = DFF.spec()
            CO2_bo = float(df['mCO2(aq)'])
            # we can leave the saltlevel as None here, because we are using
            # CO2_bo as the basis for gas concentrations
            _Enc = self.update_gases(_Enc, how=gases,
              CO2_bo=CO2_bo, H2021_saltlevel=None, GasScale=GasScale)

        return _Enc
