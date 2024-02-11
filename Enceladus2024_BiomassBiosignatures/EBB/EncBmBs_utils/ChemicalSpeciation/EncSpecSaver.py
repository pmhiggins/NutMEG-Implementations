
import EncSpecSolver as ESSolver
import numpy as np
import pandas as pd
import sys, os
from copy import deepcopy



class EncSpecSaver:
    """
    Class for performing then saving a speciation calculation across a suite of
    ocean pH values, Temperatures, and salinities at a specific pressure,
    passed in bar.

    Dissolved gases can be either fixed (do not react), or removed entirely
    by choosing a model with 'noGases' at the end of its identifier.
    """

    def __init__(self, pHfloats, Tfloats,  salts, savedirectory,
      P=1, model='pitzerPHREEQC', fixGases=False):

        self.P_bar = P
        self.P_Pa = P*1e5
        self.pHfloats = pHfloats
        self.Tfloats = Tfloats
        self.salts = salts
        self.savedirectory = savedirectory
        self.model = model
        self.fixGases = fixGases

        self.this_dir = os.path.dirname(__file__)+'/../../data/speciation/'+savedirectory
        try:
            os.mkdir(self.this_dir)
        except:
            pass

        gf = ''
        if fixGases:
            gf = '_H2CH4fixed_'

        self.filename = 'spec_'+str(int(P))+'bar_'+model+gf+'.csv'



    def save_spec_params(self, prefix='', alt_dir=None):
        """
        Generate and save a suite of speciation calculations.
        """

        _l = []


        for k, s in enumerate(self.salts):
            for i, pH_bo in enumerate(self.pHfloats):
                for j, T in enumerate(self.Tfloats):
                    _d = ESSolver.get_spec_params(T, round(pH_bo, 1), s,
                      P=self.P_Pa,
                      model=self.model,
                      fixGases=self.fixGases)
                    _d['T'] = T
                    _d['pH_bo'] = pH_bo
                    _d['salt_lvl'] = s
                    _l.append(_d)

        df = pd.DataFrame(_l)
        if alt_dir == None:
            df.to_csv(self.this_dir+'/'+prefix+self.filename)
        else:
            df.to_csv(self.alt_dir+'/'+prefix+self.filename)



    def save_spec_params_1salt(self, salt, prefix=''):
        """
        Generate and save a suite of speciation calculations at only one salt
        concentration, with those calculations organised into directories by
        salt concentration.

        """

        saltdir = self.this_dir+'/Clconc_'+str(salt)
        try:
            os.mkdir(saltdir)
        except:
            pass

        onesalt_ESSav = deepcopy(self)
        onesalt_ESSav.salts = [salt]
        onesalt_ESSav.this_dir = saltdir
        onesalt_ESSav.save_spec_params(prefix=prefix)
