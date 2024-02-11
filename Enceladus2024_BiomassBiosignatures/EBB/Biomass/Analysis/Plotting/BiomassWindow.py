# this will help allow relative imports for code which is run within the package
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__)+'/../')

from BMDRTO_MC import BMDRTO_MC
from BMDRTO_utils import BMDRTO_utils

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import math

mpl.rcParams['font.family'] = 'Liberation sans'#'sans-serif'


def BMstraightlines(model, pHvals, Tvals, Tdef, param='BM',
  GasScale=1.0, datadir='pH7to10',
  axnums=['A','B','C'], preamble = ''):

    fig, axs = plt.subplots(nrows=1,ncols=3, figsize=(12,5.5))

    # extract results from MC simulation
    _BMDRTO_MC = BMDRTO_MC([model], [Tdef], pHvals, Tvals, None, GasScale, datadir)
    BMdict = _BMDRTO_MC.load_df_as_dicts(model ,Tdef, param='BM')
    DGdict = _BMDRTO_MC.load_df_as_dicts(model ,Tdef, param='DeltaG_M')

    # mean empirical metabolic rates
    maxmet = [_BMDRTO_MC.get_maxmet(T) for T in Tvals]

    # we will need to process the default BMDRTO dictionaries to find the
    # means and max/min ranges, so set up equivalent dictionaries to monitor them
    pminNEL = {pH : {T : 0. for T in _BMDRTO_MC.Trange} for pH in _BMDRTO_MC.pHrange}
    pminEL = {pH : {T : 0. for T in _BMDRTO_MC.Trange} for pH in _BMDRTO_MC.pHrange}
    pmaxEL = {pH : {T : 0. for T in _BMDRTO_MC.Trange} for pH in _BMDRTO_MC.pHrange}

    pmeanEL = {pH : {T : 0. for T in _BMDRTO_MC.Trange} for pH in _BMDRTO_MC.pHrange}
    # mean of all habitable scenarios
    pmeanALL = {pH : {T : 0. for T in _BMDRTO_MC.Trange} for pH in _BMDRTO_MC.pHrange}

    DGmin = {pH : {T : 0. for T in _BMDRTO_MC.Trange} for pH in _BMDRTO_MC.pHrange}
    DGmax = {pH : {T : 0. for T in _BMDRTO_MC.Trange} for pH in _BMDRTO_MC.pHrange}

    PMs = BMDRTO_utils.get_TOM_Esynth_PMs(Tvals) # dict in form PMs[Tdef][Trange list]

    for i in range(len(_BMDRTO_MC.Trange)):

        try:
            pminEL[_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]] = np.array(
              BMdict['EL'][_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]]).min()

            pmaxEL[_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]] = np.array(
              BMdict['EL'][_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]]).max()
        except ValueError:
            # this typically happends because this array is empty, so set it to nan
            pminEL[_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]] =  np.nan
            pmaxEL[_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]] =  np.nan

        try:
            pminNEL[_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]] = np.array(
              BMdict['NEL'][_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]]).min()
        except ValueError:
            # this is also often because this array is empty, so set it to the min
            # of the energy limited result instead
            pminNEL[_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]] = np.array(
              BMdict['EL'][_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]]).min()

        pmeanEL[_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]] = np.array(
          BMdict['EL'][_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]]).mean()

        pmeanALL[_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]] = np.array(
          BMdict['allhab'][_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]]).mean()

        # for the DeltaG approximation, use the total range, without a
        # kinetic habitability assessment
        DGmax[_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]] = np.array(
          DGdict['full'][_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]]).max()

        DGmin[_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]] = np.array(
          DGdict['full'][_BMDRTO_MC.pHrange[i]][_BMDRTO_MC.Trange[i]]).min()


    for i, pH in enumerate(pHvals):

        # pull out temperature lists of biomass
        _p_minEL = np.array(list(pminEL[np.float64(pH)].values()))
        _p_maxEL = np.array(list(pmaxEL[np.float64(pH)].values()))
        _p_meanEL = np.array(list(pmeanEL[np.float64(pH)].values()))

        _p_minNEL = np.array(list(pminNEL[np.float64(pH)].values()))
        _p_meanALL = np.array(list(pmeanALL[np.float64(pH)].values()))

        _DG_min = np.array(list(DGmin[np.float64(pH)].values()))
        _DG_max = np.array(list(DGmax[np.float64(pH)].values()))


        axs[i].fill_between(Tvals, 10**_p_minEL, 10**_p_maxEL, alpha=0.5, facecolor='none', hatch='X', edgecolor='tab:blue', label='Range of energy limited distribution')
        axs[i].plot(Tvals, 10**_p_minNEL, alpha=1., color='tab:blue', label='Min. of not-energy limited biomass')

        axs[i].fill_between(Tvals, 0.25/(np.array(maxmet)*10), 0.25/(np.array(maxmet)*0.1), alpha=0.5, facecolor='none', hatch='O', edgecolor='tab:orange', label='Energy limitation threshold')


        axs[i].plot(Tvals, 10**_p_meanEL, color='tab:blue', label='Mean of energy limited distribution', ls='dotted')
        axs[i].plot(Tvals, 10**_p_meanALL, color='tab:blue', label='Mean of entire habitable distribution', ls='dashed')

        axs[i].fill_between(Tvals, -0.25*_DG_min/PMs['Lever2pc'], -0.25*_DG_max/PMs['Lever2pc'], alpha=0.5, facecolor='tab:red', label='Alternative biomass estimate (using \n'+'minimal maintenance power $\hat{P}^{2\%}_{M}$)')
        axs[i].fill_between(Tvals, -0.25*_DG_min/PMs['TOM'], -0.25*_DG_max/PMs['TOM'], alpha=0.5, facecolor='tab:brown', label='Alternative biomass estimate (using \n'+'maximal maintenance power $\hat{P}^{TOM}_{M}$)')

        axs[i].set_yscale('log')

        axs[i].set_ylabel(r'Biomass per mol. H$_\mathregular{2}$ consumed [cells mol$^{-1}$ s$^{-1}$]')
        axs[i].set_xlabel('Temperature [K]')
        axs[i].set_title('Bulk ocean pH: '+str(pH))
        axs[i].set_xlim(273.15,393.15)
        axs[i].set_ylim(1e3,1e29)


    fig.subplots_adjust(top=0.95,
      bottom=0.33,
      left=0.06,
      right=0.987,
      hspace=0.2,
      wspace=0.24)

    axs[1].legend(loc='lower center',ncol=3, bbox_to_anchor=(0.5,-0.53), handlelength=4)
    for a, _i in zip(axs, axnums):
        a.text(0.02,0.98, _i, ha='left', va='top', fontweight='bold', fontsize=20, transform = a.transAxes)


    plt.savefig('Figures/'+preamble+'_BiomassWindow_'+Tdef+'.pdf')
    plt.savefig('Figures/'+preamble+'_BiomassWindow_'+Tdef+'.png', dpi=300.)

    plt.close()



def __main__():
    datadir='pH7to10'
    Tdef = 'Lever2pc'
    model = 'pitzerPHREEQCnoGases'
    pHvals = [7.,8.,9.]
    Tvals = np.linspace(273.15, 393.15, num=13)
    axnums = ['A', 'B','C']
    GasScale=1.0

    BMstraightlines(model, pHvals, Tvals, Tdef, GasScale=GasScale,
      datadir=datadir, axnums=axnums, preamble='F4_')

    pHvals = [7.,8.]
    Tdef = 'TOM'
    axnums = ['A', 'B','C']
    BMstraightlines(model, pHvals, Tvals, Tdef,  GasScale=GasScale,
      datadir=datadir, axnums=axnums, preamble='FS3_top_')

    GasScale=10.

    pHvals = [7.,8.,9.]
    Tdef = 'Lever2pc'
    axnums = ['D', 'E','F']
    BMstraightlines(model, pHvals, Tvals, Tdef,  GasScale=GasScale,
      datadir=datadir, axnums=axnums, preamble='FS3_middle_')

    pHvals = [7.,8.,9.]
    Tdef = 'TOM'
    axnums = ['G', 'H','I']
    BMstraightlines(model, pHvals, Tvals, Tdef,  GasScale=GasScale,
      datadir=datadir, axnums=axnums, preamble='FS3_bottom_')





__main__()
