# this will help allow relative imports for code which is run within the package
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import pandas as pd
import numpy as np

import sys, os
sys.path.append(os.path.dirname(__file__)+'/../')

from EncBmBs_utils import EncSpecRetriever

class DataFrameFetcher:

    def __init__(self,
        Ts = np.linspace(273.15, 473.15, num=21),
        pHs = np.linspace(7,12, num=11),
        salts = ['low', 'nom', 'high'],
        P=1,
        model='pitzerPHREEQCnoGases',
        fixGases=True,
        dr='spec_T_273-473_pH_7-12'):

        self.Ts = Ts
        self.pHs = pHs
        self.salts = salts
        self.P = P
        self.model = model
        self.fixGases = fixGases
        self.dr = dr # name of the data directory of this dataframe


    def spec(self, prefix = ''):
        """
        Return a pandas DataFrame containing carbonate speciation data.
        """
        _df = None

        if len(self.salts)==1 and type(self.salts[0])==type(np.float64(0.)):
            _df = EncSpecRetriever.retrieve_1salt(self.salts[0], P=self.P, model=self.model, dr=self.dr,
              fixGases=self.fixGases, prefix=prefix)
            __df = _df[_df['T'].isin(self.Ts)]
            return __df[__df['pH_bo'].isin(self.pHs)]
        else:
            _df = EncSpecRetriever.retrieve(P=self.P, model=self.model, dr=self.dr,
              fixGases=self.fixGases, prefix=prefix)

            __df = _df[_df['T'].isin(self.Ts)]
            ___df = __df[__df['pH_bo'].isin(self.pHs)]

            return ___df[___df['salt_lvl'].isin(self.salts)]



    def microbial(self):

        ### code to extract a microbial dataframe goes here
        return

    def spec_and_microbial(self):

        ### code to get both speciation and microbial dataframe goes here
        return


    def spec_Z(self, Zname, prefix=''):
        """
        Return a pandas DataFrame containing carbonate speciation data, but only for
        Temperature, bulk ocean pH, salt level and parameter Zname.
        """
        new_df = None
        full_df = self.spec(prefix=prefix)

        if Zname[0] != 'a':
            # anything that is not an activity
            new_df = full_df.filter(['T', 'pH_bo', 'salt_lvl', Zname])
        elif Zname == 'aH' or Zname =='aOH' or Zname == 'aH2O':
            # these three activities are in the dataframe
            new_df = full_df.filter(['T', 'pH_bo', 'salt_lvl', Zname])
        else:
            # only activities should remain, calculate with gamma and molality
            _species = Zname[1:]
            _act  = full_df.filter(['T', 'pH_bo', 'salt_lvl'])
            _act['a'+_species] = full_df['g'+_species] * full_df['m'+_species]
            new_df = _act
        # catch any problems
        if len(new_df.columns) ==3:
            emsg = "Column name: '"+str(Zname)+"' not found"
            if type(Zname) != type(''):
                emsg += "', please pass a string.'"
            raise ValueError(emsg)

        return new_df
