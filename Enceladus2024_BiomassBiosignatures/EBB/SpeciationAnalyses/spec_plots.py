# this will help allow relative imports for code which is run within the package
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import matplotlib.pyplot as plt
import numpy as np
import os, sys

sys.path.append(os.path.dirname(__file__)+'/../EncBmBs_utils/')

from pH_T_LineGenerator import  pH_T_LineGenerator as pHT_LG
from DataFrameFetcher import DataFrameFetcher as DFF

import multiprocessing


""" Code to generate and save speciation output with changing bulk ocean pH and T"""


def main(P, model, fixGases=False, salts_list=[0.05,0.1,0.2],
  datadir='spec_T_273-473_pH_7-12',
  savedir = 'bulk_speciation/'):
    """
    Generate and save a series of three panel figures that plot all
    chemical speciation output data with temperature and bulk ocean pH
    values beetween 7 and 12. The salinity of each panel is given
    by the salts_list parameter.

    """

    plotting_dir = os.path.dirname(os.path.abspath(__file__))+'/figures/'+savedir+model+'_'+str(P)+'bar_fixGases_'+str(fixGases)

    try:
        os.mkdir(plotting_dir)
    except:
        pass

    salts_list = [np.float64(s) for s in salts_list]

    salts = [salts_list[0]]

    # use a DataFrameFetcher to load in the speciation data
    _DFF = DFF(model=model, salts=salts, fixGases=fixGases, P=P, dr=datadir,
      pHs=np.linspace(7,12, num=6),
      Ts = np.linspace(273.15, 403.15, num=27))
    spec_df = _DFF.spec()

    # LineGenerator for plotting
    LG = pHT_LG(_DFF)

    # do not plot these
    _ignores = ['T', 'pH_bo', 'salt_lvl']

    # get all of the output parameters - plot them all!
    _params = [str(lg) for lg in spec_df.columns[1:]]

    for Z in _params:

        if Z not in _ignores:

            fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(6,15), sharey=True, sharex=True)

            for ax, salt in zip(axs, salts_list):
                _DFF.salts = [salt]

                LG = pHT_LG(_DFF)
                ax = LG.plot_spec_lines(ax, Z, salt_lvl=salt)
                ax.set_title(r'[Cl$^{-}$]: '+str(salt)+' M')
                ax.grid()

                if Z.startswith('a') or Z.startswith('m'):
                    ax.set_yscale('log')

            fig.subplots_adjust(top=0.88, hspace=0.2, bottom=0.04, left=0.2, right=0.98)
            axs[0].legend(loc='upper center', bbox_to_anchor=[0.,0.33, 1.0, 1.0], ncol=2, handlelength=6.)

            # individual titling of plots and setting y limits
            Ztitle = ''
            if Z.startswith('a'):
                Ztitle = 'activity '+Z[1:]
                if axs[0].get_ylim()[1] > 1:
                    axs[0].set_ylim(axs[0].get_ylim()[0], 1)
            elif Z.startswith('m'):
                Ztitle = 'molality '+Z[1:]
                if axs[0].get_ylim()[1] > 10:
                    axs[0].set_ylim(axs[0].get_ylim()[0], 10)
            elif Z.startswith('g'):
                if axs[0].get_ylim()[1] > 2:
                    axs[0].set_ylim(axs[0].get_ylim()[0], 2)
                Ztitle = 'activity coefficient of '+Z[1:]
            elif Z.startswith('logSI'):
                Ztitle = 'log(saturation index of '+Z[6:]+')'
                if Z[6:] == 'CO2g' and P==1:
                    axs[0].set_ylim(-7,4)
                elif Z[6:] == 'CO2g' and P==100:
                    axs[0].set_ylim(-9,0)
                if Z[6:] == 'H2g' and P==1:
                    axs[0].set_ylim(-6,5)
                elif Z[6:] == 'H2g' and P==100:
                    axs[0].set_ylim(-8,1)
            elif Z.startswith('I'):
                Ztitle = 'Ionic Strength'
                if axs[0].get_ylim()[1] > 1:
                    axs[0].set_ylim(0, 1)
            elif Z == 'pH':
                if axs[0].get_ylim()[1] > 14:
                    axs[0].set_ylim(5, 14)
            elif Z == 'pE':
                if axs[0].get_ylim()[1] > -5:
                    axs[0].set_ylim(-5, -14)
            else:
                Ztitle = Z

            plt.suptitle(Ztitle+'\n'+model+'; Pressure: '+str(P)+' bar')

            try:
                plt.savefig(plotting_dir+'/'+str(Z)+'.pdf')
            except:
                print('/figures/'+savedir+model+'_'+str(P)+'bar_fixGases_'+str(fixGases)+'                 FAILED')
            plt.close()


""" implementation """
def plot_all(cores):
    P1s = [1,1,1,1]
    P100s = [100,100,100,100]
    models = ['SUPCRTnoGases','pitzerPHREEQCnoGases', 'SUPCRT', 'pitzerPHREEQC']
    fixGases = [False, False, True, True]


    datadir= ['spec_T_273-473_pH_7-12',]*cores
    MainSalt = [0.05,0.1,0.2]
    HighSalt = [0.2,0.30314,0.4]

    for sd in ['bulk_speciation/', 'bulk_speciation_HighSalt/']:
        try:
            os.mkdir(os.path.dirname(os.path.abspath(__file__))+'/figures/'+sd)
        except:
            pass

    with multiprocessing.Pool(cores) as pool:
        for Ps in [P1s, P100s]:
            for _sl, sd in zip([MainSalt, HighSalt], ['bulk_speciation/', 'bulk_speciation_HighSalt/']):
                pool.starmap(main, zip(Ps, models, fixGases, [_sl]*cores, datadir, [sd,]*cores))

# amend number of cores as appropriate
plot_all(4)
