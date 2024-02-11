# this will help allow relative imports for code which is run within the package
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__)+'/../')

from EncTOMSystem import EncTOMSystem as ETS
from BMDRTO_utils import BMDRTO_utils
import numpy as np
import pandas as pd
from itertools import chain

import random
import multiprocessing
import math

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


class BMDRTO_MC:
    """
    Class for performing a series of Monte Carlo simulations to predict the
    habitability, steady state biomass and turnover potential of an Enceladus-
    ocean-like environment.

    Attributes
    ----------
    spec_models : list
        list of different speciation models to use. Usually this will just be one.
    Tdef_names : list
        list of different microbial maintenance power identifiers to use.
    pHvals : list or array-like
        1D list of bulk ocean pH values to use.
    Tvals : list or array-like
        1D list of temperature values to use.
    Clvals : list or array-like
        1D list of input [Cl] (determining [DIC]) values to use. In the MC,
        rather than a uniform distribution, a random entry in this list is chosen.
    GasScale : float
        Scaling for the concentrations of H2 and CH4. When GasScale=1., they are
        set based on ocean-top [CO2], and these gases ratio with that in the
        Enceladus space plume (Waite et al 2017)
    savedir : string
        name for directory to save this output in.
    cap_mr : boolean
        if True, when estimated metabolic rate is larger than the mean value
        observed in, limit it to that value. When False, allow the estimated
        metabolic rate to take its original value.
    cores : int
        Number of cores to use when performing the MC.
    spec_dr : string
        Identifier for the directory of the carbonate speciation data
    fixGases : boolean
        When True, the carbonate speciation does not interact with dissolved gases.
    """


    _recordsnames = ['H2offset', 'CH4offset', 'saltlevel', 'koffset',
      'DeltaG_M', 'nATP','BM', 'TO', 'DR', 'mmr_exceeded', 'dead']

    # Cls = np.logspace(math.log10(0.05),math.log10(0.2), num=21)
    # Cls = np.round(Cls, decimals=5) # for wide_pH, and ideally future dfs when doing 1salt

    # These are lifted from the TOM class, c.f Higgins & Cockell 2020
    maxmet_poly = [-4.07231645e-17,  1.32059394e-01, -8.10101975e+01]
    maxGR_poly = [0.06602829, -31.84363583]

    def __init__(self, spec_models, Tdef_names, pHvals, Tvals, Clvals, GasScale, savedir,
      cap_mr=False, cores=1,
      spec_dr = 'spec_T_273-473_pH_7-12',
      fixGases=False):

        self.spec_models = spec_models
        self.savedir = savedir
        self.Tdef_names = Tdef_names
        self.cap_mr = cap_mr
        self.cores = cores
        self.pHvals = pHvals
        self.Tvals = Tvals
        self.Clvals = Clvals
        self.spec_dr = spec_dr
        self.fixGases = fixGases
        self.GasScale = GasScale

        pHrange, Trange = np.meshgrid(pHvals, Tvals)
        self.pHrange = pHrange.flatten()
        self.Trange = Trange.flatten()




    @staticmethod
    def roll(count, T, pH, Tdef, Cls, model,
      fixGases,
      GasScale,
      spec_dr,
      cap_mr):
        """
        Perform a Monte Carlo simluation of size `count` across a meshgrid of T and pH.

        Values will be chosen from uniform distributions over uncertainties in H2
        concentration, CH4 concentration, ocean salinty, and metabolic rate constant.
        Other parameters will be fixed as the values passed to this function.
        """
        records = {k:v for v, k in enumerate(BMDRTO_MC._recordsnames)}

        # make a 2D array to store data in
        # rows are count long, there are len(records) of them
        data = np.zeros((len(BMDRTO_MC._recordsnames), count))

        for i in range(int(count)):
            _sl = random.choice(Cls)
            _k = random.uniform(-1.,1.)
            _H2 = random.uniform(-1.,1.)
            _CH4 = random.uniform(-1.,1.)

            data[records['H2offset']][i] = _H2
            data[records['CH4offset']][i] = _CH4
            data[records['saltlevel']][i] = _sl
            data[records['koffset']][i] = _k

            EncParams={
              'T':T, 'pH_bo':float(pH),
              'gases':{'H2':_H2, 'CH4':_CH4},
              'saltlevel':_sl,
              'CarbonatesModel':model,
              'fixGases':fixGases,
              'spec_dr':spec_dr,
              'GasScale':GasScale}
            TOMParams={
              'Tdef':Tdef, 'k_corr':_k,
              'cap_mr':cap_mr}

            _ETS = ETS(1, EncParams, TOMParams)
            _ETS.update_steadystate_BMDRTO()

            data[records['BM']][i] = math.log10(_ETS.ssBM)
            data[records['nATP']][i] = _ETS.TOM.respiration.n_ATP
            data[records['DeltaG_M']][i] = _ETS.TOM.respiration.G_A


            if _ETS.ssDR == 0.:
                data[records['dead']][i] = 1
                data[records['TO']][i] = np.nan
                data[records['DR']][i] = np.nan
            else:
                data[records['TO']][i] = math.log10(_ETS.ssTO)
                data[records['DR']][i] = math.log10(_ETS.ssDR)

                # check how the simulation compares to the max empirically
                # observed metabolic rate
                _mm = BMDRTO_MC.get_maxmet(T)

                if _mm <= _ETS.TOM.respiration.rate:
                    # note that this TOM is respiring faster than mean empirical observations
                    data[records['mmr_exceeded']][i] = 1

        df = pd.DataFrame(data.T, columns = BMDRTO_MC._recordsnames)

        return df

    @staticmethod
    def get_maxGR(T):
        """
        return the best-fit exponential growth rate of methanogens from
        Higgins & Cockell 2020
        """
        return math.exp(BMDRTO_MC.maxGR_poly[0]*(T)+BMDRTO_MC.maxGR_poly[1])

    @staticmethod
    def get_maxmet(T):
        """
        return  the best-fit exponential metabolic rate of methanogens
        from Higgins & Cockell 2020
        """
        return math.exp(BMDRTO_MC.maxmet_poly[0]*(T*T)+BMDRTO_MC.maxmet_poly[1]*(T)+BMDRTO_MC.maxmet_poly[2])

    def save_df(self, count):
        """
        Run and save a MC simuation for biomass, death rate, and turnover.
        """
        a1 = None
        num = [count,]*len(self.Trange)

        # set up the directories for saving data
        for T, pH in zip(self.Trange, self.pHrange):
            BMDRTO_utils.get_dirname(T, pH, self.savedir)


        for model in self.spec_models:
            for Tdef in self.Tdef_names:
                # to successfully pool across cores, we need a list
                # for each parameter of the roll function that matches
                # the dimensionality of Trange and pHrange.
                Tdefs = len(self.Trange)*[Tdef,]
                spec_models = len(self.Trange)*[model,]
                cap_mrs = len(self.Trange)*[self.cap_mr,]
                fixGases_s = len(self.Trange)*[self.fixGases]
                spec_drs = len(self.Trange)*[self.spec_dr]
                Cls_s = len(self.Trange)*[self.Clvals]
                GasScale_s = len(self.Trange)*[self.GasScale]

                # this generates the sets of data over Trange and pHrange
                with multiprocessing.Pool(self.cores) as pool:
                    a1 = pool.starmap(BMDRTO_MC.roll, zip(num,
                      self.Trange, self.pHrange, Tdefs, Cls_s, spec_models,
                       fixGases_s,
                       GasScale_s,
                       spec_drs,
                       cap_mrs))
                # a1 is a list of DataFrames full of output data!
                # save them in the directories we made earlier.
                for T, pH, a in zip(self.Trange, self.pHrange, a1):
                    a.to_csv(BMDRTO_utils.get_filename(model, T, pH, Tdef, self.GasScale, self.savedir, cap_mr=self.cap_mr))


    def load_df(self, model, Tdef, pH, T):
        """
        load an output DataFrame from an MC previously performed
        """
        for _v, _l in zip([model, Tdef, pH, T],
          [self.spec_models, self.Tdef_names, self.pHvals, self.Tvals]):
            if _v not in _l:
                raise ValueError(str(model)+' not found in this BMDRTO_MC object! It only contains '+str(_l))

        _df = pd.read_csv(BMDRTO_utils.get_filename(model, T, pH, Tdef, self.GasScale, self.savedir, cap_mr=self.cap_mr))
        return _df


    def load_df_as_dicts(self, model, Tdef, param='BM'):
        """
        Return a series of nested dictionaries that each are populated by
        a specific output parameter (see below) across this BMDRTO_MC object's
        Temperature and bulk_ocean pH density.

        output parameters can take the names in BMDRTO_MC._recordsnames, or
        just read from the output csv files.

        returns a dictionary. Each value starting with f_ is a float that
        characterizes the fraction of the sample that meets a specific criteria.
        All other values are lists with a sample of param
        that meets various criteria. The dictionary has the following keys:
          'full': the full dataset with no conditions
          'EL' : subset which is energy limited
          'NEL' : subset which is not energy limited (i.e., metabolic rates
            would be able to exceed mean values in exponential growth in the
            lab - something other than energy is likely to slow this down
            in situ)
          'UIH' : subset which is uninhabitable
          'f_EL' : fraction which is energy limited,
          'f_NEL' : fraction which is not energy limited,
          'f_UIH' : fraction which is uninhabitable}

        """
        # full dataset with no conditions
        p = {pH : {T : [] for T in self.Trange} for pH in self.pHrange}
        # Energy limited: neither GR is too big nor MR is too big, and not dead
        pEL = {pH : {T : [] for T in self.Trange} for pH in self.pHrange}
        # NOT energy limited: either GR is too big or MR is too big
        pNEL = {pH : {T : [] for T in self.Trange} for pH in self.pHrange}
        # uninhabitable: no growth
        pUIH = {pH : {T : [] for T in self.Trange} for pH in self.pHrange}
        # all habitable settings (EL + NEL)
        p_allhab = {pH : {T : [] for T in self.Trange} for pH in self.pHrange}

        _frac_UIH = {pH : {T : 0. for T in self.Trange} for pH in self.pHrange}
        _frac_NEL = {pH : {T : 0. for T in self.Trange} for pH in self.pHrange}
        _frac_EL = {pH : {T : 0. for T in self.Trange} for pH in self.pHrange}

        for i in range(len(self.Trange)):
            _maxGR = math.log10(BMDRTO_MC.get_maxGR(self.Trange[i]))
            _df = self.load_df(model, Tdef, self.pHrange[i], self.Trange[i])

            _df_UIH = _df[_df['dead']==1.]
            _df_EL = _df[(_df['mmr_exceeded'] == 0.) & (_df['dead'] == 0.) & (_df['DR'] < _maxGR) ]
            _df_NEL = _df[(_df['mmr_exceeded'] == 1.) | (_df['DR'] > _maxGR)]

            p[self.pHrange[i]][self.Trange[i]].extend(_df[param].tolist())
            pUIH[self.pHrange[i]][self.Trange[i]].extend(_df_UIH[param].tolist())
            pEL[self.pHrange[i]][self.Trange[i]].extend(_df_EL[param].tolist())
            pNEL[self.pHrange[i]][self.Trange[i]].extend(_df_NEL[param].tolist())

            p_allhab[self.pHrange[i]][self.Trange[i]].extend(_df_EL[param].tolist())
            p_allhab[self.pHrange[i]][self.Trange[i]].extend(_df_NEL[param].tolist())

            tot = len(_df.index)
            tot_UIH = len(_df_UIH.index)
            tot_EL = len(_df_EL.index)
            tot_NEL = len(_df_NEL.index)

            _frac_UIH[self.pHrange[i]][self.Trange[i]] = round(100 * tot_UIH / tot, 2)
            _frac_NEL[self.pHrange[i]][self.Trange[i]] = round(100 * (tot_NEL) / tot,2)
            _frac_EL[self.pHrange[i]][self.Trange[i]] = round(100 * tot_EL / tot, 2)

        return {'full':p,
          'EL' : pEL,
          'NEL' : pNEL,
          'UIH' : pUIH,
          'allhab' : p_allhab,
          'f_EL' : _frac_EL,
          'f_NEL' : _frac_NEL,
          'f_UIH' : _frac_UIH}
