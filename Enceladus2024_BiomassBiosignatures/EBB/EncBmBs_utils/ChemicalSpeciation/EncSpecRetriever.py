import pandas as pd
import os

class EncSpecRetriever:
    """
    Simple class for retrieving speciation dataframes.

    Speciation must first be performed and saved using the EncSpecSaver class using an
    environment with reaktoro v2.

    """


    def retrieve(P=1, model='pitzerPHREEQCnoGases', dr='rkt', fixGases=True,
      prefix=''):

        thisdir = os.path.dirname(__file__)+'/../../data/speciation/'+dr

        gf = ''
        if fixGases:
            gf = '_H2CH4fixed_'

        try:
            df = pd.read_csv(thisdir+'/'+prefix+'spec_'+str(int(P))+'bar_'+model+gf+'.csv')
            return df
        except FileNotFoundError:
            emsg = 'No speciation exists with this pressure, this should first be created using EncSpecSolver.save_spec_params(P) and reaktoro v2'
            raise FileNotFoundError(emsg)


    def retrieve_1salt(salt, P=1, model='pitzerPHREEQCnoGases', dr='rkt',
      fixGases=True, prefix=''):

        thisdir = os.path.dirname(__file__)+'/../../data/speciation/'+dr+'/Clconc_'+str(salt)

        gf = ''
        if fixGases:
            gf = '_H2CH4fixed_'

        try:
            df = pd.read_csv(thisdir+'/'+prefix+'spec_'+str(int(P))+'bar_'+model+gf+'.csv')
            return df
        except FileNotFoundError:
            emsg = 'No speciation exists with this pressure, this should first be created using EncSpecSolver.save_spec_params(P) and reaktoro v2'
            raise FileNotFoundError(emsg)
