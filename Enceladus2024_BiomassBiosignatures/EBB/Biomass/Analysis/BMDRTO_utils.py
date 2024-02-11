# this will help allow relative imports for code which is run within the package
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__)+'/../')

from EncTOMSystem import EncTOMSystem as ETS
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import numpy as np
import pandas as pd
from copy import deepcopy
import random
import multiprocessing
import time, math

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

class BMDRTO_utils:
    """
    Helper class with usseful methods for the Monte Carlo implementation.
    """

    def get_dirname(T, pH_bo, savedir):
        """
        Create the directory savedir in the data/biomass directory for this
        temperature and pH_bo combination
        """
        datadir = os.path.dirname(__file__)+'/../../data/biomass/'
        try:
            os.mkdir(datadir+savedir)
        except FileExistsError:
            pass
        except:
            raise
        subdir = '/'+str(int(T))+'_'+str(pH_bo)
        try:
            os.mkdir(datadir+savedir+subdir)
        except FileExistsError:
            pass
        except:
            raise

        return datadir+savedir+subdir


    def get_filename(model, T, pH_bo, Tdef, GasScale, savedir, cap_mr=False):
        """
        get the filepath for the results under the conditions passed.
        It will have the format: biomass/T_pH_bo/model_Tdef_mmr
        """
        suffix=model+'_'+Tdef
        if GasScale != 1.:
            suffix += '_GS'+str(int(GasScale))
        if cap_mr:
            suffix += '_mmr'
        suffix += '.csv'
        subdir = '/'+str(int(T))+'_'+str(pH_bo)
        return os.path.dirname(__file__)+'/../../data/biomass/'+savedir + subdir + '/' + suffix

    def get0(_):
        return 0*_


    def get_TOM_Esynth_PMs(Trange):

        """
        Get a dictionary of TOM Esynths and maintenance powers, each in the
        temperature range Trange.
        """

        dim = len(Trange)
        res = {'E_synth' : np.zeros(dim),
          'TOM' : np.zeros(dim),
          'TijhuisAnaerobe' : np.zeros(dim),
          'Lever2pc' : np.zeros(dim),
          'Lever10pc' :  np.zeros(dim)}
        for i, T in enumerate(Trange):
            EncParams={
              'T':T, 'pH_bo':8.,
              'gases':{'H2':0., 'CH4':0.},
              'saltlevel':np.float64(0.1),
              'fixGases':False}

            for Tdef in ['TOM', 'TijhuisAnaerobe', 'Lever2pc', 'Lever10pc']:
                TOMParams={
                  'Tdef':Tdef, 'k_corr':0.,
                  'cap_mr':False}

                _ETS = ETS(1, EncParams, TOMParams)
                res['E_synth'][i] = _ETS.TOM.E_synth
                res[Tdef][i] = deepcopy(_ETS.TOM.maintenance.net_dict['T'])

        return res
