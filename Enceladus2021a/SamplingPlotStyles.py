import sys, os, ast, math
sys.path.append(os.path.dirname(__file__)+'../../NutMEG')
sys.path.append(os.path.dirname(__file__)+'../../NutMEG-Implementations/TOM')
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

import Sampler
import EnceladusPlotStyles as EPS
from colormapping import cmapper

"""
This file is for setting up and adding lines to plots which involve
sampling over the parameter space. Not all of these were used in the manuscript
but may come in handy for future work.
"""


# can remove this if you import enceladusplotstyles
mpl.rcParams['contour.negative_linestyle'] = 'solid'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cmr10'
mpl.rcParams['font.serif'] = 'Times.tcc'
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rc('axes', unicode_minus=False)
mpl.rcParams['hatch.linewidth'] = 2.0
# print(mpl.rcParams.keys()) # to see what can be changed
mpl.rcParams['font.size'] = 12
# mpl.rcParams['font.weight'] = 'bold'


def nominalline_T(ax, pH, c, Trange=[273,400], salttype='nom', zeroed=True, ls='-',
  fixupdate={}, fn_preamble=''):
    """ Plot a line of nominal power supplies with temperature at a given pH """

    nT, nTPS = Sampler.nominal_T(pH, salttype=salttype, Templims=Trange, zeroed=zeroed,
      fixupdate=fixupdate, fn_preamble=fn_preamble)
    ax.plot(nT, nTPS, c=c, ls=ls)
    return ax

def nominalline_pH(ax, T, c, pHrange=[7,12], salttype='nom', zeroed=True, ls='-',
  fixupdate={}, fn_preamble=''):
    """ Plot a line of nominal power supplies with pH at a given temperature """

    npH, npHPS = Sampler.nominal_pH(T, salttype=salttype, pHlims=pHrange, zeroed=zeroed,
      fixupdate=fixupdate, fn_preamble=fn_preamble)
    ax.plot(npH, npHPS, c=c, ls=ls)
    return ax




def samplebins_T(ax, pH, samplesize, cmap, Trange=[273,400], salttype='nom',
  ylims=[-30,-5], bins=[45,45],
  fixupdate={}, fn_preamble=''):
    """
    plot sampling results as a histogram with T at a given pH
    (for plotting over a nonimal line)
    """

    Tvals, PS = Sampler.sampleT(samplesize, pH=pH, salttype=salttype,
      fixupdate=fixupdate, fn_preamble=fn_preamble)

    zeroes = (np.log10(PS) != -50.)
    y = np.log10(PS)[zeroes]
    x = Tvals[zeroes]

    # ax[0].hexbin(x, y, gridsize=(int(pHnum/2), int(pHnum/8)), cmap=cmaps[i])
    ax.hist2d(x,y, bins=bins, range=[Trange, ylims], cmap=cmap, edgecolor=None)

    return ax

def samplebins_pH(ax, T, samplesize, cmap, pHrange=[7,12], salttype='nom',
  ylims=[-30,-5], bins=[45,45],
  fixupdate={}, fn_preamble=''):
    """
    plot sampling results as a histogram with pH at a given T
    (for plotting over a nonimal line)
    """

    pHvals, PS = Sampler.samplepH(samplesize, Temp=T, pHrange=pHrange, salttype=salttype,
      fixupdate=fixupdate, fn_preamble=fn_preamble)

    zeroes = (np.log10(PS) != -50.)
    y = np.log10(PS)[zeroes]
    x = pHvals[zeroes]
    # ax[0].hexbin(x, y, gridsize=(int(pHnum/2), int(pHnum/8)), cmap=cmaps[i])
    ax.hist2d(x,y, bins=bins, range=[pHrange, ylims], cmap=cmap, edgecolor=None)
    return ax


def all_varianceplot_pH(ax, samplesize, T, salttype='nom', pHrange=[8,12], fixupdate={}, fn_preamble='', cm='Blues'):
    """
    Plot the total variance in power supply in log space given by the
    parameter space as a histogram. Do so for a given T at variable pH.
    """

    nx, nPS = Sampler.nominal_pH(T, salttype=salttype, pHlims=pHrange, zeroed=False,
      fixupdate=fixupdate, fn_preamble=fn_preamble)

    nPS = np.log10(nPS)

    xvals, PS = Sampler.samplepH(samplesize, Temp=T, pHrange=pHrange, salttype=salttype,
      fixupdate=fixupdate, fn_preamble=fn_preamble)

    _logy = np.log10(PS)
    zeroes = (_logy != -50.)
    logy = _logy[zeroes]
    x = xvals[zeroes]
    logvariance = np.empty(len(logy))

    ic = 0
    for xi, logyi in zip(x, logy):
        close_nom = int((xi-pHrange[0]) / (pHrange[1]-pHrange[0]) * len(nx))
        logvariance[ic] = logyi - nPS[close_nom]
        ic+=1
    ax.hist2d(x,logvariance, bins=[45,45], range=[pHrange, [-4, 4]], cmap=cm, edgecolor=None)

    return ax


def all_varianceplot_T(ax, samplesize, pH, salttype='nom', Trange=[273,400], fixupdate={}, fn_preamble='', cm='Blues'):
    """
    Plot the total variance in power supply in log space given by the
    parameter space as a histogram. Do so for a given pH at variable T.
    """

    nx, nPS = Sampler.nominal_T(pH, salttype=salttype, Templims=Trange, zeroed=False,
      fixupdate=fixupdate, fn_preamble=fn_preamble)

    nPS = np.log10(nPS)

    xvals, PS = Sampler.sampleT(samplesize, pH=pH, Trange=Trange, salttype=salttype,
      fixupdate=fixupdate, fn_preamble=fn_preamble)

    _logy = np.log10(PS)
    zeroes = (_logy != -50.)
    logy = _logy[zeroes]
    x = xvals[zeroes]
    logvariance = np.empty(len(logy))
    xreds = []
    vreds = []

    ic = 0
    for xi, logyi in zip(x, logy):
        close_nom = int((xi-Trange[0]) / (Trange[1]-Trange[0]) * len(nx))
        logvariance[ic] = logyi - nPS[close_nom]
        if logvariance[ic] > 10:
            vreds.append(0)
            xreds.append(x[ic])
        ic+=1
    hb = ax.hist2d(x,logvariance, bins=[45,45], range=[Trange, [-4, 4]], cmap='Blues', edgecolor=None, vmin=0, vmax=samplesize/200)
    hr = ax.hist2d(xreds, vreds, bins=[45,45], range=[Trange, [-4, 4]], cmap=cmapper.r2a(), edgecolor=None, vmin=0, vmax=samplesize/100)
    return ax, hb, hr





def nominal2x3(Ts=[275, 300, 325], pHs=[8,9,10], save='figs/nominal2x3.pdf', show=False, maintenance=True, nATPchoice=1.0):
    """
    Plot the nominal power supply and endmembers of salt content for the
    temperatures in Ts with pH and the pHs in pHs with temperature.
    Can be adjusted for different choices of nATP.
    This also plots some characteristic maintenance powers
    """

    fig, ax = plt.subplots(figsize=(9,12), ncols=2, nrows=3)

    nomcols = [mpl.cm.get_cmap('autumn')(i*0.15) for i in range(0,5)]
    maincols = [mpl.cm.get_cmap('winter')(i*0.33) for i in range(0,4)]

    for i, T in enumerate(Ts):

        # ax[i][0].set_title('T = '+str(T)+' K')
        ax[i][0].text(7.2, -8, 'T = '+str(T)+' K', fontsize=14)
        ax[i][0] = nominalline_pH(ax[i][0], T, nomcols[0], ls='dotted', salttype='high', fixupdate={'H2':1., 'CH4':-1., 'nATP':nATPchoice, 'k_corr':1.0}, fn_preamble='maxed_nATP'+str(nATPchoice))
        ax[i][0] = nominalline_pH(ax[i][0], T, nomcols[1], ls='dashed', salttype='high', fixupdate={'nATP':nATPchoice}, fn_preamble='nATP_'+str(nATPchoice))
        ax[i][0] = nominalline_pH(ax[i][0], T, nomcols[2], fixupdate={'nATP':nATPchoice}, fn_preamble='nATP_'+str(nATPchoice))
        ax[i][0] = nominalline_pH(ax[i][0], T, nomcols[3], ls='dashed', salttype='low', fixupdate={'nATP':nATPchoice}, fn_preamble='nATP_'+str(nATPchoice))
        ax[i][0] = nominalline_pH(ax[i][0], T, nomcols[4], ls='dotted', salttype='low', fixupdate={'H2':-1., 'CH4':1., 'nATP':nATPchoice, 'k_corr':-1.0}, fn_preamble='mined_nATP'+str(nATPchoice))

        ax[i][0].set_xlabel('wider ocean pH')
        ax[i][0].set_ylabel('Power [W/cell]')
        if maintenance:
            EPS.add_pH_maintenance_lines(ax[i][0], T, colors=maincols)
        ax[i][0].set_xlim(7,12)

    for j, pH in enumerate(pHs):

        # ax[j][1].set_title('wider ocean pH = '+str(pH), loc='left')
        ax[j][1].text(275, -8, 'wider ocean pH = '+str(pH), fontsize=14)
        ax[j][1] = nominalline_T(ax[j][1], pH, nomcols[0], ls='dotted', salttype='high', fixupdate={'H2':1., 'CH4':-1., 'nATP':nATPchoice, 'k_corr':1.0}, fn_preamble='maxed_nATP'+str(nATPchoice))
        ax[j][1] = nominalline_T(ax[j][1], pH, nomcols[1], ls='dashed', salttype='high', fixupdate={'nATP':nATPchoice}, fn_preamble='nATP_'+str(nATPchoice))
        ax[j][1] = nominalline_T(ax[j][1], pH, nomcols[2], fixupdate={'nATP':nATPchoice}, fn_preamble='nATP_'+str(nATPchoice))
        ax[j][1] = nominalline_T(ax[j][1], pH, nomcols[3], ls='dashed', salttype='low', fixupdate={'nATP':nATPchoice}, fn_preamble='nATP_'+str(nATPchoice))
        ax[j][1] = nominalline_T(ax[j][1], pH, nomcols[4], ls='dotted', salttype='low', fixupdate={'H2':-1., 'CH4':1., 'nATP':nATPchoice, 'k_corr':-1.0}, fn_preamble='mined_nATP'+str(nATPchoice))


        ax[j][1].set_xlabel('Temperature [K]')
        ax[j][1].set_ylabel('Power [W/cell]')
        ax[j][1].yaxis.set_label_position("right")
        ax[j][1].yaxis.tick_right()
        if maintenance:
            EPS.add_maintenance_lines(ax[j][1], colors=maincols)
        ax[j][1].set_xlim(270,400)
    if maintenance:
        ax[0][1] = EPS.add_maintenance_labels(ax[0][1], colors=maincols)

    ax[0][0].plot([0],[0], c=nomcols[2], label='Power supply in a nominal-salt ocean')
    ax[0][0].plot([0],[0], c=nomcols[1], ls='dashed', label='Nominal power supply in a high-salt ocean')
    ax[0][0].plot([0],[0], c=nomcols[0], ls='dotted', label='Highest power supply in a high-salt ocean')
    ax[0][0].plot([0],[0], c=nomcols[3], ls='dashed', label='Nominal power supply in a low-salt ocean')
    ax[0][0].plot([0],[0], c=nomcols[4], ls='dotted', label='Lowest power supply in a low-salt ocean')

    ax[0][0].legend(bbox_to_anchor=(0., 1.5, 2.05, .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
    ax[0][1].legend(bbox_to_anchor=(-0.8, 1.05, 1.55, .102), loc=3,
           ncol=1, mode="expand", borderaxespad=0.)
    for a in ax.flatten():
        # a.set_yscale('log')
        a.set_ylim(-30, -5)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.05)
    plt.savefig(save)
    if show:
        plt.show()
    return fig, ax





def binsexample():
    """ Example of plotting histograms over nominal power supply lines """
    fig, ax = plt.subplots(figsize=(10,5), ncols=2, nrows=1)

    # ax[0] = samplebins_T(ax[0], 8, 1000, cmapper.k2a())
    ax[0] = samplebins_T(ax[0], 8, 1000, cmapper.b2a(), salttype='low')
    ax[0] = samplebins_T(ax[0], 8, 1000, cmapper.r2a(), salttype='high')
    ax[0] = nominalline_T(ax[0], 8, 'b', salttype='low')
    ax[0] = nominalline_T(ax[0], 8, 'r', salttype='high')

    # ax[1] = samplebins_pH(ax[1], 300, 100, cmapper.k2a())
    ax[1] = samplebins_pH(ax[1], 300, 1000, cmapper.b2a(), salttype='low')
    ax[1] = samplebins_pH(ax[1], 300, 1000, cmapper.r2a(), salttype='high')
    ax[1] = nominalline_pH(ax[1], 300, 'b', salttype='low')
    ax[1] = nominalline_pH(ax[1], 300, 'r', salttype='high')

    plt.show()






def ones_varianceplot(T=300, pH=None, salttype='nom', samplesize=100,
  Trange=[273,400], pHrange=[7,12], save=True, show=True, cm='Blues',
  fixupdate={}, fn_preamble=''):
    """
    Plot the independent variance of the secondary effects in the parameter
    space: nATP, aCH4, aH2 and k.
    """

    fig, ax = plt.subplots(figsize=(8,8), ncols=2, nrows=2)
    ax = ax.flatten()

    minx,maxx = 0, 0
    nx, nPS = 0,0

    if pH == None:
        minx, maxx = pHrange[0], pHrange[1]
        nx, nPS = Sampler.nominal_pH(T, salttype=salttype, pHlims=pHrange, zeroed=False,
          fixupdate=fixupdate, fn_preamble=fn_preamble)
    elif T == None:
        minx, maxx = Trange[0], Trange[1]
        nx, nPS = Sampler.nominal_T(pH, salttype=salttype, Templims=Trange, zeroed=False,
          fixupdate=fixupdate, fn_preamble=fn_preamble)

    nPS = np.log10(nPS)

    ind = ['CH4', 'H2' , 'nATP', 'k_corr']
    for i, v in enumerate(ind):
        xvals, PS = Sampler.independent_sample(T=T, pH=pH, to_sample=v, samplesize=samplesize, salttype=salttype,
          Trange = Trange, pHrange=pHrange)

        if pH == None:
            ax[i].set_xlabel('pH')
        elif T == None:
            ax[i].set_xlabel('Temperature [K]')

        if save or show:
            _logy = np.log10(PS)
            zeroes = (_logy != -50.)
            logy = _logy[zeroes]
            x = xvals[zeroes]
            logvariance = np.empty(len(logy))
            xreds = []

            ic = 0
            for xi, logyi in zip(x, logy):
                close_nom = int(((xi-minx) / (maxx-minx)) * len(nx))
                logvariance[ic] = logyi - nPS[close_nom]
                if logvariance[ic] > 10:
                    xreds.append(x[ic])
                ic+=1
            ax[i].hist2d(x,logvariance, bins=[45,45], range=[[minx, maxx], [-4, 4]], cmap=cm, edgecolor=None, vmax=samplesize/100)
            ax[i].hist2d(xreds,[0]*len(xreds), bins=[45,45], range=[[minx, maxx], [-4, 4]], cmap=cmapper.r2a(), edgecolor=None, vmax=samplesize/100)

            # ax[i].scatter(x,logvariance, alpha=0.1, c='b', edgecolor=None)
            # ax[i].set_xlim([minx, maxx])
            # ax[i].set_ylim([-4, 4])

            ax[i].set_title(v)
            ax[i].set_ylabel('Variance in log10 (power supply)')

    if pH == None:
        plt.suptitle('Variance with pH when T = '+str(T)+', '+salttype+' salt')
    elif T == None:
        plt.suptitle('Variance with T when pH = '+str(pH)+', '+salttype+' salt')
    plt.tight_layout()
    if save:
        plt.savefig('figs/_indy_samples/'+str(samplesize)+'_T_'+str(T)+'_pH_'+str(pH)+'_salt_'+salttype+'.pdf')
    if show:
        plt.show()
    plt.close()







###### may not be needed
def varianceexample():
    """ Plot overall variance at select pH and salt levels """
    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
    ax[0][0] = all_varianceplot_T(ax[0][0], 100, 8, salttype='high')
    ax[0][1] = all_varianceplot_T(ax[0][1], 100, 9, salttype='high')
    ax[0][1] = all_varianceplot_T(ax[0][2], 100, 10, salttype='high')

    ax[1][0] = all_varianceplot_T(ax[1][0], 100, 8, salttype='nom')
    ax[1][1] = all_varianceplot_T(ax[1][1], 100, 9, salttype='nom')
    ax[1][1] = all_varianceplot_T(ax[1][2], 100, 10, salttype='nom')

    ax[2][0] = all_varianceplot_T(ax[2][0], 100, 8, salttype='low')
    ax[2][1] = all_varianceplot_T(ax[2][1], 100, 9, salttype='low')
    ax[2][1] = all_varianceplot_T(ax[2][2], 100, 10, salttype='low')

    # ax[1][0] = all_varianceplot_pH(ax[1][0], 100, 300)
    # ax[1][1] = all_varianceplot_pH(ax[1][1], 100, 300, salttype='high')
    plt.tight_layout()
    plt.savefig('allvtest.pdf')
    plt.show()

# varianceexample()
