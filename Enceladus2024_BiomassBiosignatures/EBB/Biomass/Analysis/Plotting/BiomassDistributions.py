# this will help allow relative imports for code which is run within the package
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__)+'/../')

from BMDRTO_MC import BMDRTO_MC
from BMDRTO_utils import BMDRTO_utils

import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd

from scipy import stats as spys
from itertools import chain

mpl.rcParams['font.family'] = 'Liberation sans'#'sans-serif'


def get_binspace(param):
    """ for plotting histograms and KDEs, get a reasonable list of bins. """
    bs = {'BM': np.linspace(5,25, num=61),
      'DR': np.linspace(-15,10, num=51),
      'TO': np.linspace(10,14, num=41)}
    return bs[param]

def get_label(param):
    """ For identifier param, get corresponding axes label"""
    lbls = {
    'BM': '$\log_{10}$(steady state biomass per mole of \n H$_{2}$ consumption [cells / (mol H$_{2}$ / s)])',
    'TO': '$\log_{10}$(biomass turnover per mole of \n H$_{2}$ consumption [cells / s (mol H$_{2}$ / s)])',
    'DR': '$\log_{10}$(death rate [/s])'
    }
    return lbls[param]

def get0(_):
    """ Function that returns zero, for null output of density. """
    return 0*_


def getKDE_and_boxwidth(paramlst, minx, maxx, bw_method=1.):
    """
    perform a KDE approximation on distribution paramlst, between
    minx and maxx. Adjust bin width to use for KDE using bw_method.
    return the density function from the scipy KDE fit and the list of x values.
    """

    # in case paramlst is multidimensional, flatten it
    _all = list(chain.from_iterable(paramlst))
    xs = np.arange(minx, maxx,.05)

    if len(_all) ==0:
        # no values at all, return 0 everytime.
        return get0, xs, [1, [maxx]]

    density = spys.gaussian_kde(_all, bw_method=bw_method*(len(paramlst)**(-0.2)))

    return density, xs

def plot_2D_pooled_Tdefs(ax, param='BM', model='pitzerPHREEQCnoGases',
  Tdef='Lever2pc', pH=8., GasScale=1.0, datadir='pH7to10'):

    Tvals = np.linspace(273.15, 393.15, num=7)

    # load MC output
    _BMDRTO_MC = BMDRTO_MC([model], [Tdef], [pH], Tvals, None, GasScale, datadir)
    BMdict = _BMDRTO_MC.load_df_as_dicts(model ,Tdef, param=param)

    Tbreak = Tvals[1] - Tvals[0] # temperature gap

    binspace = get_binspace(param) # bin edges
    binbreak = binspace[1] - binspace[0] # bin width

    cmap = mpl.cm.coolwarm

    for T_i, T in enumerate(Tvals):

        if T in BMdict['allhab'][pH]:

            _density, xs = getKDE_and_boxwidth([BMdict['allhab'][pH][T]], binspace[0], binspace[-1], bw_method=0.25)

            # probability of habitability
            prob_hab = (1-(BMdict['f_UIH'][pH][T]/100))
            ax.plot(xs, _density(xs), c=cmap(T_i/(len(Tvals)-1))) #*prob_not_dead

            # label probability of habitability above top of distribution.
            _idensitymax = np.argmax(_density(xs))

            ax.text(xs[_idensitymax], _density(xs[_idensitymax]),
              str(math.floor(prob_hab*100))+'%', c=cmap(T_i/(len(Tvals)-1)),
              ha='center', va='bottom')

    # labelling
    ax.text(0.02,0.98,'Bulk ocean pH: '+str(int(pH)), ha='left', va='top', transform=ax.transAxes)
    if Tdef == 'Lever2pc':
        ax.text(0.02,0.93,'Maintenance: Minimal $(\hat{P}_M^{2\%})$', ha='left', va='top', transform=ax.transAxes)
    elif Tdef == 'TOM':
        ax.text(0.02,0.93,'Maintenance: Maximal $(\hat{P}_M^{TOM})$', ha='left', va='top', transform=ax.transAxes)

    ax.set_xlabel(get_label(param))
    ax.set_ylabel('Habitable-space probability density')
    return ax


# code to generate Figure 3 in Higgins et al., 2024 JGR:Planets
def __main__():
    fig, axs = plt.subplots(nrows=3,ncols=2, figsize=(9,12), constrained_layout=True)
    axs = axs.flatten()

    axs[0] = plot_2D_pooled_Tdefs(axs[0], param='BM', model='pitzerPHREEQCnoGases',
      Tdef='Lever2pc', pH=8., GasScale=1.0)
    axs[1] = plot_2D_pooled_Tdefs(axs[1], param='TO', model='pitzerPHREEQCnoGases',
      Tdef='Lever2pc', pH=8., GasScale=1.0)
    axs[2] = plot_2D_pooled_Tdefs(axs[2], param='BM', model='pitzerPHREEQCnoGases',
      Tdef='Lever2pc', pH=9., GasScale=1.0)
    axs[3] = plot_2D_pooled_Tdefs(axs[3], param='TO', model='pitzerPHREEQCnoGases',
      Tdef='Lever2pc', pH=9., GasScale=1.0)
    axs[4] = plot_2D_pooled_Tdefs(axs[4], param='BM', model='pitzerPHREEQCnoGases',
      Tdef='TOM', pH=8., GasScale=1.0)
    axs[5] = plot_2D_pooled_Tdefs(axs[5], param='TO', model='pitzerPHREEQCnoGases',
      Tdef='TOM', pH=8., GasScale=1.0)

    axs[0].set_ylim(-0.02,0.28)
    axs[3].set_ylim(-0.2,4.1)


    for ax, L in zip(axs, ['A', 'B', 'C', 'D', 'E', 'F']):
        ax.text(0.98,0.98,L, ha='right', va='top', fontweight='bold', fontsize=18,transform=ax.transAxes)

    # make a consistent colormap at the bottom
    cmap=mpl.cm.coolwarm
    clist = [cmap(T_i/(6)) for T_i in range(7)]
    cb = fig.colorbar(
      plt.cm.ScalarMappable(norm=mpl.colors.BoundaryNorm(np.linspace(263.15, 403.15, num=8), cmap.N), cmap=cmap),
      ax=axs, shrink=0.5, label='Temperature [K]', orientation='horizontal',
      anchor=(0.5,-0.2), pad=+0.01,
      ticks=np.linspace(273, 393, num=7))

    plt.savefig('Figures/F3_BiomassTurnoverDistributions.pdf')
    plt.savefig('Figures/F3_BiomassTurnoverDistributions.svg')
    plt.savefig('Figures/F3_BiomassTurnoverDistributions.png', dpi=300)

    plt.close()


__main__()
