import matplotlib.pyplot as plt
import numpy as np
import os, sys

sys.path.append(os.path.dirname(__file__)+'/../EncBmBs_utils/')

from pH_T_LineGenerator import  pH_T_LineGenerator as pHT_LG
from DataFrameFetcher import DataFrameFetcher as DFF

import multiprocessing




def main(P, model, fixGases=False, salts_list=[0.05,0.1,0.2],
  datadir='spec_T_273-473_pH_7-12',
  savedir = 'bulk_speciation/'):

    plotting_dir = os.path.dirname(__file__)+'/figures/'+savedir+model+str(P)+'bar_'+str(fixGases)
    try:
        os.mkdir(plotting_dir)
    except:
        pass

    salts_list = [np.float64(s) for s in salts_list]

    salts = [salts_list[0]]

    _DFF_SUPCRT =  DFF(model='SUPCRT'+model, salts=salts,fixGases=fixGases, P=P, dr=datadir,
      pHs=np.linspace(7,12, num=6),
      Ts = np.linspace(283.15, 373.15, num=20))

    _DFF_pitzer =  DFF(model='pitzerPHREEQC'+model, salts=salts, fixGases=fixGases, P=P, dr=datadir,
      pHs=np.linspace(7,12, num=6),
      Ts = np.linspace(283.15, 373.15, num=20))

    s_df = _DFF_SUPCRT.spec()
    LG_SUPCRT = pHT_LG(_DFF_SUPCRT)
    LG_pitzer = pHT_LG(_DFF_pitzer)

    # ignore these parameters
    _ignores = ['T', 'pH_bo', 'salt_lvl', 'mNaCl(aq)', 'gNaCl(aq)', 'aNaCl(aq)',
      'mNaOH(aq)', 'gNaOH(aq)', 'aNaOH(aq)',
      'mH2(aq)', 'gH2(aq)', 'aH2(aq)',
      'mCH4(aq)', 'gCH4(aq)', 'aCH4(aq)']

    #otherwise, generate a comparison figure for each
    _params = [str(lg) for lg in s_df.columns[1:]]

    for Z in _params:

        if Z not in _ignores:

            fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(6,15), sharey=True, sharex=True)

            for ax, salt in zip(axs, salts_list):
                _DFF_SUPCRT.salts = [salt]
                _DFF_pitzer.salts = [salt]

                # we cannot use the LineGenerator here, so compute manually

                Sdf = _DFF_SUPCRT.spec()
                Sdf = Sdf[Sdf['salt_lvl'] == salt]

                pdf = _DFF_pitzer.spec()
                pdf = pdf[pdf['salt_lvl'] == salt]

                for _i, _pH in enumerate(_DFF_SUPCRT.pHs):
                    _pdf = pdf[pdf['pH_bo']==_pH]
                    _sdf = Sdf[Sdf['pH_bo']==_pH]

                    _pdf['ratio'] = _pdf[Z] / _sdf[Z]
                    if Z == 'pH':
                        # for pH, instead show pitzer - SUPCRT
                        _pdf['ratio'] = _pdf[Z] - _sdf[Z]
                    ax.plot(_pdf['T'], _pdf['ratio'], c=LG_SUPCRT.cmaplist[_i], label='O$\degree$C pH = '+str(_pH), ls=LG_SUPCRT.lineslist[_i])

                ax.set_xlabel('Temperature [K]')
                ax.set_ylabel(Z+' ratio (pitzer / SUPCRT)')
                if Z == 'pH':
                    ax.set_ylabel('pH (via pitzer) - pH (via SUPCRT)')
                ax.set_title('[Cl$^{-}$]: '+ str(salt))
                if Z != 'pH':
                    ax.set_yscale('log')
                ax.grid(which='both')

            # fig.set_constrained_layout(False)
            fig.subplots_adjust(top=0.88, hspace=0.2, bottom=0.04, left=0.2, right=0.98)
            axs[0].legend(loc='upper center', bbox_to_anchor=[0.,0.33, 1.0, 1.0], ncol=2, handlelength=6.)

            Ztitle = ''
            if Z.startswith('a'):
                Ztitle = 'activity '+Z[1:]
            elif Z.startswith('m'):
                Ztitle = 'molality '+Z[1:]
            elif Z.startswith('g'):
                Ztitle = 'activity coefficient of '+Z[1:]
            elif Z.startswith('logSI'):
                Ztitle = 'log(saturation index of '+Z[1:]+')'
            elif Z.startswith('I'):
                Ztitle = 'Ionic Strength'
            else:
                Ztitle = Z

            plt.suptitle(Ztitle+' model ratios \n'+model+' Pressure: '+str(P)+' bar', y=0.995)

            try:
                plt.savefig(plotting_dir+'/'+str(Z)+'.pdf')
            except:
                print('/figures/'+savedir+model+'_'+str(P)+'bar_fixGases_'+str(fixGases)+'                 FAILED')
            plt.close()

def plot_all():

    P1s = [1,1]
    P100s = [100,100]

    models = ['', 'noGases']
    fixGases = [True, False]

    datadir= ['spec_T_273-473_pH_7-12',]*2
    # savedir = 'bulk_speciation/'
    MainSalt = [0.05,0.1,0.2]
    HighSalt = [0.2,0.30314,0.4]

    for sd in ['speciation_ratios/', 'speciation_raios_HighSalt/']:
        try:
            os.mkdir(os.path.dirname(os.path.abspath(__file__))+'/figures/'+sd)
        except:
            pass

    with multiprocessing.Pool(2) as pool:
        for Ps in [P1s, P100s]:
            for _sl, sd in zip([MainSalt, HighSalt], ['speciation_ratios/', 'speciation_raios_HighSalt/']):
                pool.starmap(main, zip(Ps, models, fixGases, [_sl]*2, datadir, [sd,]*2))

plot_all()
