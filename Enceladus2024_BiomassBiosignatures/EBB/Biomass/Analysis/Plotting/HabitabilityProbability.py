# this will help allow relative imports for code which is run within the package
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__)+'/../')

from BMDRTO_MC import BMDRTO_MC
from BMDRTO_utils import BMDRTO_utils

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd

""" Plotting the habitability probability figures found in Higgins+ 2024 """



def HabProb(model, Tdef, pHvals, Tvals, Clvals, GasScale, datadir='v7b', preamble='8to9',
  labels=['a', 'b', 'c'], pHlabel=False):

    fig, axs = plt.subplots(nrows=1,ncols=3, figsize=(12,4.7))

    # to retrieve MC output
    _BMDRTO_MC = BMDRTO_MC([model], [Tdef], pHvals, Tvals, Clvals, GasScale, datadir)

    # load MC results as dictionaries
    BMdict = _BMDRTO_MC.load_df_as_dicts(model ,Tdef, param='BM')

    # extract habitability fractions
    _D = pd.DataFrame(BMdict['f_UIH']) # uninhabitable
    _D_NEL = pd.DataFrame(BMdict['f_NEL']) # not-energy limited
    _D_EL = pd.DataFrame(BMdict['f_EL']) # energy-limited

    # linestyles
    lss = ['solid', 'dotted', (0,(1,3,1,3,1,3,1,3)), (0,(3,3,1,3,1,3,1,3)), (0,(3,3,3,3,1,3,1,3)), (0,(3,3,3,3,3,3,1,3)), 'dashed']

    cmap_D = mpl.cm.autumn
    cmap_EL = mpl.cm.winter
    cmap_NEL = mpl.cm.summer

    cmapdim = int(len(pHvals)*4/3)

    for i, pH in enumerate(pHvals):
        axs[0].plot(Tvals, _D[np.float64(pH)], c=cmap_D(i/cmapdim), ls=lss[i])
        axs[1].plot(Tvals, _D_EL[np.float64(pH)], c=cmap_EL(i/cmapdim), ls=lss[i])
        axs[2].plot(Tvals, _D_NEL[np.float64(pH)], c=cmap_NEL(i/cmapdim), ls=lss[i])
        axs[2].plot(280, -0.5, c='k', ls=lss[i], label=r'pH$_\mathregular{bo}$: '+str(pH))

        if pHlabel:
            axs[0].text( 268, _D[np.float64(pH)].tolist()[0], str(pH),
              ha='right', va='center', color=cmap_D(i/cmapdim))
            axs[1].text( 268, _D_EL[np.float64(pH)].tolist()[0], str(pH),
              ha='right', va='center', color=cmap_EL(i/cmapdim))
            axs[2].text( 398, _D_NEL[np.float64(pH)].tolist()[-1], str(pH),
              ha='left', va='center', color=cmap_NEL(i/cmapdim))

            axs[0].set_xlim(250, 393)
            axs[1].set_xlim(250, 393)
            axs[2].set_xlim(273,415)


    axs[0].set_title('Probability uninhabitable')
    axs[1].set_title('Probability habitable; \n'+'energy limits growth rate')
    axs[2].set_title('Probability habitable; \n'+'energy does not limit growth rate')


    for i, ax in enumerate(axs):
        ax.set_ylim(-2,109)
        ax.set_xlabel('Temperature [K]', fontsize=16)
        ax.set_ylabel('Probability [%]', fontsize=16)
        ax.text(-0.05,1.01, labels[i], ha='center', va='bottom', fontweight='bold', fontsize=20, transform = ax.transAxes)
        ax.grid()

    fig.subplots_adjust(top=0.905,
      bottom=0.22,
      left=0.053,
      right=0.987,
      hspace=0.2,
      wspace=0.210)

    fig.legend(loc='lower center',ncol=7, handlelength=5)
    plt.savefig('Figures/'+preamble+'_Habprob_'+Tdef+'.svg')
    plt.savefig('Figures/'+preamble+'_Habprob_'+Tdef+'.pdf')
    plt.savefig('Figures/'+preamble+'_Habprob_'+Tdef+'.png', dpi=200)
    
    plt.close()



def __main__():

    # the top three panels of Figure 1
    datadir='pH8to9_test'
    model = 'pitzerPHREEQCnoGases'
    Tdef = 'Lever2pc'
    preamble = 'F1_top_8to9'
    pHvals = np.linspace(8,9,num=6)
    Tvals = np.linspace(273.15, 393.15, num=25)
    Clvals =None
    GasScale = 1.0
    HabProb(model, Tdef, pHvals, Tvals, Clvals, GasScale, datadir='pH8to9', labels=['A','B','C'], preamble=preamble,pHlabel=True)

    """ pH 7 to 10 for figure 1 lower panel and Figures S1 & S2"""

    pHvals = [7.0,7.5,8.0,8.5,9.0,9.5,10.0]
    Tvals = np.linspace(273.15, 393.15, num=13)
    datadir='pH7to10'
    Clvals =None

    # lower three panels of Figure 1
    preamble= 'F1_bottom_7to10'
    GasScale = 1.0
    HabProb(model, Tdef, pHvals, Tvals, Clvals, GasScale, datadir=datadir,
      labels=['D','E','F'], preamble=preamble,pHlabel=True)

    # the top three panels of Figure S1
    preamble= 'FS1_top_7to10'
    GasScale = 1.0
    HabProb(model, Tdef, pHvals, Tvals, Clvals, GasScale, datadir=datadir,
      labels=['A','B','C'], preamble=preamble,pHlabel=False)
    # the bottom three panels of Figure S1
    preamble= 'FS1_bottom_7to10'
    Tdef = 'TOM'
    GasScale = 1.0
    HabProb(model, Tdef, pHvals, Tvals, Clvals, GasScale, datadir=datadir,
      labels=['D','E','F'], preamble=preamble,pHlabel=False)


    # figure S2 uses the GasScaling as 10x
    preamble= 'FS2_top_7to10'
    Tdef = 'Lever2pc'
    GasScale = 10.0
    HabProb(model, Tdef, pHvals, Tvals, Clvals, GasScale, datadir=datadir,
      labels=['A','B','C'], preamble=preamble,pHlabel=False)

    preamble= 'FS2_bottom_7to10'
    Tdef = 'TOM'
    GasScale = 10.0
    HabProb(model, Tdef, pHvals, Tvals, Clvals, GasScale, datadir=datadir,
      labels=['D','E','F'], preamble=preamble,pHlabel=False)




__main__()
