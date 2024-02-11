# this will help allow relative imports for code which is run within the package
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__))

from BMDRTO_MC import BMDRTO_MC
from BMDRTO_utils import BMDRTO_utils

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as spys
import pandas as pd
import math
import csv


model = 'pitzerPHREEQCnoGases'
Tdefs = ['Lever2pc', 'TOM']
GasScales = [1.0,10.0]

pHvals =  [7.0,7.5,8.0,8.5,9.0,9.5,10.0]

Tvals = np.linspace(273.15, 393.15, num=13)

datadir = 'pH7to10'

savedir = os.path.dirname(__file__)+'/../../data/BMTO_means_CIs/pH7to10_tex'
try:
    os.mkdir(savedir)
except FileExistsError:
    pass


for GasScale in GasScales:
    _BMDRTO_MC = BMDRTO_MC([model], Tdefs, pHvals, Tvals, None, GasScale, datadir)
    for Tdef in Tdefs:

        BMdict = _BMDRTO_MC.load_df_as_dicts(model, Tdef, param='BM')
        TOdict = _BMDRTO_MC.load_df_as_dicts(model ,Tdef, param='TO')
        for pi, pH in enumerate(pHvals):
            this_dict = {'Temperature':Tvals,
              'Mean BM HAB':np.empty(len(Tvals), dtype=object),
              '95pc CI BM HAB':np.empty(len(Tvals), dtype=object),
              'Mean BM EL':np.empty(len(Tvals), dtype=object),
              '95pc CI BM EL':np.empty(len(Tvals), dtype=object),
              'Mean BM NEL':np.empty(len(Tvals), dtype=object),
              '95pc CI BM NEL':np.empty(len(Tvals), dtype=object),
              'Mean TO HAB':np.empty(len(Tvals), dtype=object),
              '95pc CI TO HAB':np.empty(len(Tvals), dtype=object),
              'Mean TO EL':np.empty(len(Tvals), dtype=object),
              '95pc CI TO EL':np.empty(len(Tvals), dtype=object),
              'Mean TO NEL':np.empty(len(Tvals), dtype=object),
              '95pc CI TO NEL':np.empty(len(Tvals), dtype=object),
              }


            for ti, T in enumerate(Tvals):

                TO_NEL = np.array(TOdict['NEL'][pHvals[pi]][Tvals[ti]])
                TO_EL = np.array(TOdict['EL'][pHvals[pi]][Tvals[ti]])
                TO_both = np.concatenate((TO_NEL, TO_EL))

                BM_NEL = np.array(BMdict['NEL'][pHvals[pi]][Tvals[ti]])
                BM_EL = np.array(BMdict['EL'][pHvals[pi]][Tvals[ti]])
                BM_both = np.concatenate((BM_NEL, BM_EL))

                for param in ['TO', 'BM']:
                # NEL.mean(), np.percentile(NEL, 2.5), np.percentile(NEL, 97.5)
                    if len(eval(param+'_NEL'))>0 or len(eval(param+'_EL'))>0:
                        this_dict['Mean '+param+' HAB'][ti] = str(eval(param+'_both').mean().round(2))
                        this_dict['95pc CI '+param+' HAB'][ti] = str(np.percentile(eval(param+'_both'), 2.5).round(2))+'-'+str(np.percentile(eval(param+'_both'), 97.5).round(2))
                    else:
                        this_dict['Mean '+param+' HAB'][ti] = '--'
                        this_dict['95pc CI '+param+' HAB'][ti] = '--'
                    if len(eval(param+'_NEL'))>0:
                        this_dict['Mean '+param+' NEL'][ti] = str(eval(param+'_NEL').mean().round(2))
                        this_dict['95pc CI '+param+' NEL'][ti] =  str(np.percentile(eval(param+'_NEL'), 2.5).round(2))+'-'+str(np.percentile(eval(param+'_NEL'), 97.5).round(2))
                    elif len(eval(param+'_NEL'))==0:
                        this_dict['Mean '+param+' NEL'][ti] = '--'
                        this_dict['95pc CI '+param+' NEL'][ti] = '--'
                    if len(eval(param+'_EL'))>0:
                        this_dict['Mean '+param+' EL'][ti] = str(eval(param+'_EL').mean().round(2))
                        this_dict['95pc CI '+param+' EL'][ti] =  str(np.percentile(eval(param+'_EL'), 2.5).round(2))+'-'+str(np.percentile(eval(param+'_EL'), 97.5).round(2))
                    elif len(eval(param+'_EL'))==0:
                        this_dict['Mean '+param+' EL'][ti] = '--'
                        this_dict['95pc CI '+param+' EL'][ti] = '--'

            df = pd.DataFrame(this_dict)

            texfn = savedir+datadir+'_pH'+str(pH)+'_'+Tdef+'.tex'
            if GasScale != 1.0:
                texfn = savedir+datadir+'_GS'+str(int(GasScale))+'_pH'+str(pH)+'_'+Tdef+'.tex'

            with open(texfn, 'w') as f:
                f.write(df.to_latex(index=False,
                  column_format='ccccccccccccc'))
