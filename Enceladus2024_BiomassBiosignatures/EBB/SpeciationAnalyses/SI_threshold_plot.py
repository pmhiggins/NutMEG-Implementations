# this will help allow relative imports for code which is run within the package
import math, os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d as i1d
from scipy.ndimage.filters import gaussian_filter

sys.path.append(os.path.dirname(__file__)+'/../EncBmBs_utils/')

from DataFrameFetcher import DataFrameFetcher as DFF

# globally set dimensionality of the grid we will use
EMsalts = np.logspace(math.log10(0.05),math.log10(0.4), num=31)
EMsalts = np.round(EMsalts, decimals=5)
pHs=np.linspace(7,12, num=11)
pHs = np.round(pHs, decimals=1)
Ts = np.linspace(273.15, 403.15, num=27)

datadir = 'spec_T_273-473_pH_7-12'

gas = {'logSI_H2g': 'H$_{2}$',
  'logSI_CO2g' : 'CO$_{2}$'}


def add_line(P, model, axs, param, clr='k',fixGases=True):

    salts = [round(EMsalts[0],5)]
    _DFF = DFF(model=model, salts=salts, fixGases=fixGases, P=P,
      pHs=pHs, dr=datadir,
      Ts = Ts)

    spec_df = _DFF.spec()

    pH_lim_273K = []
    pH_lim_333K = []

    T_lim_pH8 = []
    T_lim_pH9 = []

    EMgridpH, pHgrid = np.meshgrid(EMsalts, pHs)
    EMgridT, Tgrid = np.meshgrid(EMsalts, Ts)

    pH_lim_273Kgrid =[]
    pH_lim_333Kgrid = []

    T_lim_pH8grid = []
    T_lim_pH9grid = []


    for salt in EMsalts:
        _DFF.salts = [round(salt,5)]

        # get a dataframe containing T, pH_bo, salt_lvl, and the param you want
        df = _DFF.spec_Z(param)

        pH_lim_273Kgrid.append(df.query('T == 273.15')[param])
        pH_lim_333Kgrid.append(df.query('T == 333.15')[param])
        T_lim_pH8grid.append(df.query('pH_bo == 8.')[param])
        T_lim_pH9grid.append(df.query('pH_bo == 9.')[param])


    # we need the transpose of the sat index meshgrid
    pH_lim_273Kgrid = np.array(pH_lim_273Kgrid).T
    pH_lim_333Kgrid = np.array(pH_lim_333Kgrid).T
    T_lim_pH8grid = np.array(T_lim_pH8grid).T
    T_lim_pH9grid = np.array(T_lim_pH9grid).T

    # plot only the SI=0 contour line
    axs[0].contour(EMgridpH, pHgrid, pH_lim_273Kgrid, [0.], colors=[clr])
    axs[0].contour(EMgridpH, pHgrid, pH_lim_333Kgrid, [0.], colors=[clr], linestyles='dotted')

    axs[0].plot([np.nan], [np.nan], c=clr, label=gas[param]+' at 273 K,'+'\n'+'at '+str(P)+' bar pressure')
    axs[0].plot([np.nan], [np.nan], c=clr, ls='dotted', label=gas[param]+' at 333 K,'+'\n'+'at '+str(P)+' bar pressure')

    if P ==1:
        axs[1].contour(EMgridT, Tgrid, T_lim_pH8grid, [0.], colors=[clr])
        axs[1].contour(EMgridT, Tgrid, T_lim_pH9grid, [0.], colors=[clr], linestyles='dashed')

        if param=='logSI_CO2g':
            axs[1].plot([np.nan], [np.nan], c=clr, label=gas[param]+' with bulk ocean pH 8,'+'\n'+'at '+str(P)+' bar pressure')

        axs[1].plot([np.nan], [np.nan], c=clr, ls='dashed', label=gas[param]+' with bulk ocean pH 9,'+'\n'+'at '+str(P)+' bar pressure')



def __main__():

    cmap = plt.get_cmap('tab10', 10)
    cmaplist = [cmap(4), cmap(1), cmap(5)]

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(8,8))

    # add only the contour line where SI = 0 for these gases at these pressures
    # we do not need to plot CO2 at 100 bar because it is always dissolved
    add_line(1, 'pitzerPHREEQC', axs, 'logSI_H2g', clr=cmaplist[0],fixGases=True)
    add_line(1, 'pitzerPHREEQC', axs, 'logSI_CO2g', clr=cmaplist[1], fixGases=True)
    add_line(100, 'pitzerPHREEQC', axs, 'logSI_H2g', clr=cmaplist[2], fixGases=True)

    for ax in axs:
        ax.set_xlabel('Molality of Cl$^{-}$ [mol kg$^{-1}$]')
        ax.set_xlim(0.05,0.4)

    axs[0].set_ylabel('Bulk ocean pH')
    axs[1].set_ylabel('Temperature [K]')

    axs[0].set_ylim(7,10)
    axs[1].set_ylim(270,373)

    axs[0].text(0.02,1.02, 'A', ha='left', va='bottom',
      fontweight='bold', fontsize=16, transform = axs[0].transAxes)
    axs[1].text(0.02,1.02, 'B', ha='left', va='bottom',
      fontweight='bold', fontsize=16, transform = axs[1].transAxes)

    fig.suptitle("Saturation thresholds for H$_{2}$(g) and CO$_{2}$(g) in Enceladus seawater")

    plt.tight_layout()
    fig.subplots_adjust(right=0.63)

    # add legends and discretionary texts.

    leg0 = axs[0].legend(title='    Areas below plotted lines are '+'\n'+'          likely to be saturated; '+'\n'+'at 100 bar CO$_{2}$ is always unsaturated',
      bbox_to_anchor=(1.05, 0.5), loc='center left', labelspacing=0.7)

    leg0._legend_box.align = "center"

    leg1 = axs[1].legend(title='     Areas above plotted lines are '+'\n'+'           likely to be saturated; '+'\n'+' at 100 bar all gases are unsaturated; '+'\n'+ '  at pH$_{bo}$ = 8, H$_{2}$ is always saturated',
      bbox_to_anchor=(1.05, 0.5), loc='center left', labelspacing=0.7)

    leg1._legend_box.align = "center"

    plt.savefig('figures/SI_threshold.pdf')
    plt.savefig('figures/SI_threshold.png', dpi=200)

    plt.close()

__main__()
