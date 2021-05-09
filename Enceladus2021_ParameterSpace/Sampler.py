import sys, os, ast, math
sys.path.append(os.path.dirname(__file__)+'../../NutMEG')
sys.path.append(os.path.dirname(__file__)+'../../NutMEG-Implementations/TOM')

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from copy import deepcopy

import matplotlib.colors as clr
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy.stats import gaussian_kde as kde
from matplotlib.colors import Normalize
from matplotlib import cm

import EnceladusGrids as EG
import EnceladusPlotStyles as EPS

""" This is for getting power supply across samples of the parameter space """

def plot_sigmas(sample):
    """
    Plot and show the distribution given by a sample size of `sample'.
    """
    sig = EG.getsigmas(sample, fixed={'T':330, 'pH':8.5}, dummyvals=True)

    a, b, c, d, e, f = [],[],[],[], [],[]

    for i in sig:
        a.append(i[0])
        b.append(i[1])
        c.append(i[2])
        d.append(i[3])
        e.append(i[4])
        f.append(i[5])


    fig, axs = plt.subplots(nrows=3, ncols=2)

    axs = axs.flatten()

    axs[0].scatter(a,b, alpha=0.2, marker='.')
    axs[1].scatter(a,c, alpha=0.2, marker='.')
    axs[2].scatter(a, d, alpha=0.2, marker='.')
    axs[3].scatter(a, e, alpha=0.2, marker='.')
    axs[4].scatter(a, f, alpha=0.2, marker='.')
    axs[5].scatter(e, f, alpha=0.2, marker='.')

    plt.show()

# plot_sigmas(5000)

""" PH SAMPLING METHODS """


def samplepH(size, Temp=298, pHrange=[7,12], salttype='nom',
  fixupdate={}, fn_preamble=''):
    """
    Get power supplies from a sample with fixed T at Temp and a flat
    distribution across pHrange, for salt level salttype.
    If you want to update the fix dictionary (e.g. to only sample select
    parameters), send them via fixupdate.
    """

    f = {'T':Temp, 'pH':None}
    f.update(fixupdate)

    cd_preamble = 'E21data/pH_samples/'
    thissamplename = cd_preamble+fn_preamble+'pH_'+str(size)+'_'+str(Temp)+'_'+salttype+'.npy'

    try:
        xvals, PS = np.load(thissamplename, allow_pickle=True)
        return xvals, PS
    except:
        sigmas = EG.getsigmas(size, fixed=f, pHlims=pHrange, dummyvals=True)
        PS, xvals = np.empty(size), np.empty(size)

        for i, s in enumerate(sigmas):
            PS[i] = EG.getPSfromSigmas(s, dummyvals=True, salttype=salttype)
            xvals[i] = s[-1]

        np.save(thissamplename, np.array([xvals, PS]))
        return xvals, PS

def samplepH_onlyone(size, one, Temp=298, pHrange=[7,12], salttype='nom',
  uniquefix=None, fn_preamble=''):
    """
    Get power supplies from a sample with fixed T at Temp and a flat
    distribution across pHrange, for salt level salttype.
    This specifically only samples across one other variable, named by
    the attribute one.
    """
    if one=='all':
        return samplepH(size, Temp=Temp, pHrange=pHrange, salttype=salttype)

    cd_preamble = 'E21data/pH_samples/ind_variance/'
    thissamplename = cd_preamble+fn_preamble+one+'_pH_'+str(size)+'_'+str(Temp)+'_'+salttype+'.npy'

    try:
        xvals, PS = np.load(thissamplename, allow_pickle=True)
        return xvals, PS
    except:
        fixed = {}
        if uniquefix == None:
            fixed = deepcopy(EG.default_nomdict)
        else:
            fixed=uniquefix
        fixed['T'] = Temp
        fixed['pH'] = None
        fixed[one] = None
        sigmas = EG.getsigmas(size, fixed=fixed, pHlims=pHrange, dummyvals=True)
        PS, xvals = np.empty(size), np.empty(size)

        for i, s in enumerate(sigmas):
            PS[i] = EG.getPSfromSigmas(s, dummyvals=True, salttype=salttype)
            xvals[i] = s[-1]

        np.save(thissamplename, np.array([xvals, PS]))
        return xvals, PS


def nominal_pH(Temp, pHlims=[7,12], num=100, salttype='nom', zeroed=False,
  fixupdate={}, fn_preamble=''):
    """
    Get power supplies from the nominal case with fixed T at Temp and a uniform distribution across pHrange of length num, for salt level salttype.
    You can change the nominal parameter space values by passing them in
    fixupdate. Pass zeroed as True if you want the Power supplies that NutMEG
    couldn't calculate removed, otherwise ther will be 1e-50 W/cell.
    """

    cd_preamble = 'E21data/nominalPS/'
    thisnominalname = cd_preamble+fn_preamble+'T_'+str(Temp)+'_pHs_'+str(pHlims[0])+'_'+str(pHlims[1])+'_'+salttype+'.npy'

    try:
        xvals, PS = np.load(thisnominalname, allow_pickle=True)
        if zeroed:
            zeroes = (np.log10(PS) != -50.)
            PS = np.log10(PS)[zeroes]
            xvals = xvals[zeroes]
        return xvals, PS
    except:
        pHnoms = np.linspace(pHlims[0],pHlims[1],num=num)
        fix=deepcopy(EG.default_nomdict)
        fix.update(fixupdate)
        PSnoms=[]
        for pHnom in pHnoms:
            PSnoms.append(EG.getPSfromSigmas(
              [fix['CO2'],fix['H2'],fix['CH4'],fix['nATP'],fix['k_corr'],
              Temp,pHnom],
              dummyvals=True, salttype=salttype))
        np.save(thisnominalname, np.array([pHnoms, PSnoms]))
        if zeroed:
            zeroes = (np.log10(PSnoms) != -50.)
            PSnoms = np.log10(PSnoms)[zeroes]
            pHnoms = pHnoms[zeroes]
        return pHnoms, PSnoms


######## T SAMPLING METHODS

def sampleT(size, pH=8.5, Trange=[273,400], salttype='nom',
  fixupdate={}, fn_preamble=''):
    """
    Get power supplies from a sample with fixed pH at pH and a flat
    distribution across Trange, for salt level salttype.
    If you want to update the fix dictionary (e.g. to only sample select
    parameters), send them via fixupdate.
    """

    f = {'T':None, 'pH':pH}
    f.update(fixupdate)

    cd_preamble = 'E21data/T_samples/'
    thissamplename = cd_preamble+fn_preamble+'T_'+str(size)+'_'+str(pH)+'_'+salttype+'.npy'

    try:
        xvals, PS = np.load(thissamplename, allow_pickle=True)
        return xvals, PS
    except:

        sigmas = EG.getsigmas(size, fixed=f, Tlims=Trange, dummyvals=True)
        PS, xvals = np.empty(size), np.empty(size)

        for i, s in enumerate(sigmas):
            PS[i] = EG.getPSfromSigmas(s, dummyvals=True, salttype=salttype)
            xvals[i] = s[-2]

        np.save(thissamplename, np.array([xvals, PS]))
        return xvals, PS



def sampleT_onlyone(size, one, pH=8.5, Trange=[273,400], salttype='nom',
  uniquefix=None, fn_preamble=''):
    """
    Get power supplies from a sample with fixed pH at pH and a flat
    distribution across Trange, for salt level salttype.
    This specifically only samples across one other variable, named by
    the attribute one.
    """

    if one=='all':
        return(sampleT(size, Trange=Trange, pH=pH, salttype=salttype))
    cd_preamble = 'E21data/T_samples/ind_variance/'
    thissamplename = cd_preamble+fn_preamble+one+'_T_'+str(size)+'_'+str(pH)+'_'+salttype+'.npy'


    try:
        xvals, PS = np.load(thissamplename, allow_pickle=True)
        return xvals, PS
    except:
        fixed = {}
        if uniquefix == None:
            fixed = deepcopy(EG.default_nomdict)
        else:
            fixed=uniquefix
        fixed['T'] = None
        fixed['pH'] = pH
        fixed[one] = None
        sigmas = EG.getsigmas(size, fixed=fixed, Tlims=Trange, dummyvals=True)
        PS, xvals = np.empty(size), np.empty(size)

        for i, s in enumerate(sigmas):
            PS[i] = EG.getPSfromSigmas(s, dummyvals=True, salttype=salttype)
            xvals[i] = s[-2]

        np.save(thissamplename, np.array([xvals, PS]))
        return xvals, PS


def nominal_T(pH, Templims=[273,400], num=100, salttype='nom', zeroed=False,
  fixupdate=None, fn_preamble=''):
    """
    Get power supplies from the nominal case with fixed T at Temp and a uniform distribution across pHrange of length num, for salt level salttype.
    You can change the nominal parameter space values by passing them in
    fixupdate. Pass zeroed as True if you want the Power supplies that NutMEG
    couldn't calculate removed, otherwise they will be 1e-50 W/cell.
    """

    cd_preamble = 'E21data/nominalPS/'

    thisnominalname = cd_preamble+fn_preamble+'pH_'+str(pH)+'_Temps_'+str(Templims[0])+'_'+str(Templims[1])+'_'+salttype+'.npy'

    try:
        xvals, PS = np.load(thisnominalname, allow_pickle=True)
        if zeroed:
            zeroes = (np.log10(PS) != -50.)
            PS = np.log10(PS)[zeroes]
            xvals = xvals[zeroes]
        return xvals, PS
    except:
        Tnoms = np.linspace(Templims[0], Templims[1], num=num)
        fix=deepcopy(EG.default_nomdict)
        fix.update(fixupdate)
        PSnoms=[]
        for Tnom in Tnoms:
            PSnoms.append(EG.getPSfromSigmas(
              [fix['CO2'],fix['H2'],fix['CH4'],fix['nATP'],fix['k_corr'],
              Tnom,pH],
              dummyvals=True, salttype=salttype))
        np.save(thisnominalname, np.array([Tnoms, PSnoms]))
        if zeroed:
            zeroes = (np.log10(PSnoms) != -50.)
            PSnoms = np.log10(PSnoms)[zeroes]
            Tnoms = Tnoms[zeroes]
        return Tnoms, PSnoms



"""variance methods + building samples"""


def independent_sample(T=300, pH=None, to_sample='all', samplesize=100,
  Trange = [273,400], pHrange=[8,12], salttype='nom'):
    """
    Get a sample over a specific parameter (or all of them) at given T or pH.
    """
    xvals, PS = [],[]
    if pH == None:
        xvals, PS = samplepH_onlyone(samplesize, to_sample, Temp=T, pHrange=pHrange, salttype=salttype)#, pHnum=pHnum) #random sample
    elif T == None:
        xvals, PS = sampleT_onlyone(samplesize, to_sample, Trange=Trange, pH=pH, salttype=salttype)#, pHnum=pHnum) #random sample
    return xvals, PS


def get_variancelist(size, T, pH, salttype='nom',
  ones=['all', 'CH4', 'CO2','H2','nATP', 'k_corr']):
    """ for samples where only one is varied at a time. Get a list of filenames. """
    vlist=[]
    for one in ones:
        if pH == None:
            if one=='all':
                vlist.append('E21data/pH_samples/pH_'+str(size)+'_'+str(T)+'_'+salttype+'.npy')
            else:
                vlist.append('E21data/pH_samples/ind_variance/'+one+'_pH_'+str(size)+'_'+str(T)+'_'+salttype+'.npy')
        elif T==None:
            if one == 'all':
                vlist.append('E21data/T_samples/T_'+str(size)+'_'+str(pH)+'_'+salttype+'.npy')
            else:
                vlist.append('E21data/T_samples/ind_variance/'+one+'_T_'+str(size)+'_'+str(pH)+'_'+salttype+'.npy')
    return vlist

def load_variancelist(variancelist, ones=['all', 'CH4', 'CO2','H2','nATP', 'k_corr']):
    """Load saved samples.

    return an ndarry of dimensions len(ones), where each element is in the format
    [xvals, PSvals]
    """
    xvalues = []
    PSvalues = []

    for i in range(len(ones)):
        xvals, PS = np.load(variancelist[i], allow_pickle=True)
        xvalues.append(np.array(xvals))
        PSvalues.append(np.array(PS))
    return xvalues, PSvalues


def variance_chain_sample(T, pH, stepsize=100, totsize=10000, startfrom=200,
  salttype='nom',
  ones=['all', 'CH4', 'CO2','H2','nATP', 'k_corr']):
    """
    Continue performing variance samples and updating 'master' samples.
    stepsize : size of each sample incrememnt
    totsize : goal total size of sample
    startfrom : Which file to start from.

    This is designed for sampling multiple times e.g. on less powerful PCs. If
    you have the computing power or the time, it is likely better to use the
    independent_sample method
    """
    completed = startfrom
    sofar = get_variancelist(completed, T, pH, salttype=salttype, ones=ones)
    sf_x, sf_PS = load_variancelist(sofar, ones=ones)
    while completed < totsize:
        sofar = get_variancelist(completed, T, pH, salttype=salttype, ones=ones)

        for one in ones:
            independent_sample(T=T, pH=pH, to_sample=one, samplesize=stepsize,
              Trange=[273,400], salttype=salttype, pHrange=[7,12])


        justdone = get_variancelist(stepsize, T, pH, salttype=salttype, ones=ones)
        jd_x, jd_PS = load_variancelist(justdone, ones=ones)

        combined = get_variancelist(stepsize+completed, T, pH, salttype=salttype, ones=ones)

        for i in range(len(justdone)):
            sf_x[i] = np.append(sf_x[i], jd_x[i])
            sf_PS[i] = np.append(sf_PS[i], jd_PS[i])
            np.save(combined[i], np.array([sf_x[i], sf_PS[i]]))
            os.remove(sofar[i])
            os.remove(justdone[i])
        completed += stepsize
    return(sf_x, sf_PS)









##### may not be needed after all
def uniform_samplepH(size, Temp=298, pHrange=[8,12], pHnum=50, k_corr=0.0):
    thissamplename = '2_uniform_data_samples/pH_'+str(size)+'_'+str(Temp)+'.npy'
    try:
        xvals, PS = np.load(thissamplename, allow_pickle=True)
        return xvals, PS
    except:
        PS, xvals = np.empty(size), np.empty(size)
        ix=0
        for pH in np.linspace(8,12, num=pHnum):
            sigmas = EG.getsigmas(int(size/pHnum), fixed={'T':Temp, 'pH':pH}, pHlims=pHrange, dummyvals=True)

            for s in sigmas:
                PS[ix] = EG.getPSfromSigmas(s, dummyvals=True)
                xvals[ix] = s[-1]
                ix += 1

        np.save(thissamplename, np.array([xvals, PS]))
        return xvals, PS


def makeColours( vals ):
    colours = np.zeros( (len(vals),3) )
    norm = Normalize( vmin=vals.min(), vmax=vals.max() )

    #Can put any colormap you like here.
    colours = [cm.ScalarMappable( norm=norm, cmap='jet').to_rgba( val ) for val in vals]

    return colours

def colorscatter(ax, x, y, s=10, cmap='Blues'):

    zeroes = (y != -50.)
    y = y[zeroes]
    x = x[zeroes]

    xy = np.vstack([x, y])
    z = kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    ax.scatter(x, y, c=z, s=s, cmap=cmap, edgecolor=None)

    return ax

def kdecontour(ax, x, y, s=10, cmap='Blues'):

    zeroes = (y != -50.)
    y = y[zeroes]
    x = x[zeroes]

    xgrid = np.linspace(min(x)*0.9, max(x)*1.1, num=250)
    ygrid = np.linspace(-45, -5, num=250)

    xx, yy = np.meshgrid(xgrid, ygrid)
    # xy = np.vstack([x, y])
    z = np.ndarray((len(xgrid),len(ygrid)))

    _kde = kde([x,y])
    for i,xi in enumerate(xgrid):
        for j, yi in enumerate(ygrid):
            z[j][i] = _kde.evaluate([xi, yi])
    # print(math.log10(z.max()))
    ax.contourf(xx, yy, z, cmap=cmap, levels=np.logspace(-6,math.log10(z.max())))
    return ax

def sampleplot1():
    samplesize=1000
    pHvals, PSpH = samplepH(samplesize)
    Tvals, PST = sampleT(samplesize)
    pHvals105, PSpH105 = samplepH(samplesize, Temp=350)
    Tvals350, PST350 = sampleT(samplesize, pH=10.5)

    fig, ax = plt.subplots(figsize=(10,5), ncols=2, nrows=1)
    # ax[0].hexbin(pHvals, np.log10(PSpH), gridsize=30, cmap='Blues')
    # ax[1].hexbin(Tvals, np.log10(PST), gridsize=100, cmap='Blues')


    # just density
    # densObj = kde( [Tvals, np.log10(PST)] )
    # colours = makeColours( densObj.evaluate( [Tvals, np.log10(PST)] ) )
    #
    # ax[1].scatter(Tvals, np.log10(PST), color=colours)

    # density layered
    # Calculate the point density

    # logPS = np.log10(PST)
    # zeroes = (logPS != -50.)
    # logPS = logPS[zeroes]
    # Tvals = Tvals[zeroes]
    #
    # xy = np.vstack([Tvals, logPS])
    # z = kde(xy)(xy)
    #
    # # Sort the points by density, so that the densest points are plotted last
    # idx = z.argsort()
    # x, y, z = Tvals[idx], logPS[idx], z[idx]
    #
    # ax[1].scatter(x, y, c=z, s=10, edgecolor=None)
    ax[0] = colorscatter(ax[0], pHvals, np.log10(PSpH), s=25)
    ax[1] = colorscatter(ax[1], Tvals, np.log10(PST), s=25)
    ax[0] = colorscatter(ax[0], pHvals105, np.log10(PSpH105), s=25, cmap='Reds')
    ax[1] = colorscatter(ax[1], Tvals350, np.log10(PST350), s=25, cmap='Reds')


    for a in ax:
        # a.set_yscale('log')
        # a.set_ylim(-33, -7)
        a.set_ylabel('Power Supply [W/cell]')
    # plt.colorbar()
    plt.savefig('test.pdf')
# sampleplot1()

def transparent_cmap(h, ticks=50):
    colors = [clr.hsv_to_rgb((h, a/ticks, 1))
      for a in range(ticks)]
    return mpl.colors.ListedColormap(colors)

green_to_a = {'red':[(0,0,0), (1,0,0)],
  'green':[(0,0.7,0.7),(1,0.7,0.7)],
  'blue':[(0,0,0), (1,0,0)],
  'alpha':[(0,1,1), (1,0,0)]}

blue_to_a = {'red':[(0,0,0), (1,0,0)],
  'green':[(0,0,0),(1,0,0)],
  'blue':[(0,1,1), (1,1,1)],
  'alpha':[(0,1,1), (1,0,0)]}

red_to_a = {'red':[(0,1,1), (1,1,1)],
  'green':[(0,0,0),(1,0,0)],
  'blue':[(0,0,0), (1,0,0)],
  'alpha':[(0,1,1), (1,0,0)]}

def transparent_cmap2(rgba):
    return clr.LinearSegmentedColormap('Green_trans', rgba).reversed()


def sampleplot2(Ts=[400, 350,300], pHs=[8,9,10]):
    """For testing different plot types """
    samplesize=2000
    pHnum=100
    fig, ax = plt.subplots(figsize=(10,5), ncols=2, nrows=2)

    # cmaps=[transparent_cmap(238/360), transparent_cmap(120/360), 'Reds']
    cmaps=[transparent_cmap2(blue_to_a), transparent_cmap2(green_to_a), 'Reds']

    Tresults = []
    pHresults = []
    ax[0][0].set_xlim(8,12)
    ax[1][1].set_xlim(8,12)
    ax[0][1].set_xlim(8,12)
    ax[1][0].set_xlim(8,12)
    ax[1][0].set_ylim(-45,-5)
    ax[1][1].set_ylim(-45,-5)
    for i, T in enumerate(Ts):
        pHvals, PS = uniform_samplepH(samplesize, Temp=T, pHnum=pHnum)
        # ax[0][0] = kdecontour(ax[0][0], pHvals, np.log10(PS), s=25, cmap=cmaps[i])
        # ax[0][1] = colorscatter(ax[0][1], pHvals, np.log10(PS), s=25, cmap=cmaps[i])
        # ax[1][1].scatter(pHvals, np.log10(PS), c='b', s=10, alpha=0.05, edgecolor=None)
        zeroes = (np.log10(PS) != -50.)
        y = np.log10(PS)[zeroes]
        x = pHvals[zeroes]
        ax[1][0].hexbin(x, y, gridsize=(int(pHnum/2), int(pHnum/8)), cmap=cmaps[i])
    # for i, pH in enumerate(pHs):
    #     Tvals, PS = sampleT(samplesize, pH=pH)
    #     ax[1][0] = kdecontour(ax[1][0], Tvals, np.log10(PS), s=25, cmap=cmaps[i])
        # ax[1][1] = colorscatter(ax[1][1], Tvals, np.log10(PS), s=25, cmap=cmaps[i])

    for axs in ax:
        for a in axs:
            a.set_ylabel('log10 (Power Supply [W/cell])')
    ax[0][0].set_xlabel('pH')
    ax[1][0].set_xlabel('Temperature [K]')
    ax[0][0].set_xlim(8,12)
    ax[1][1].set_xlim(8,12)
    ax[0][1].set_xlim(8,12)
    ax[1][0].set_xlim(8,12)
    ax[1][0].set_ylim(-45,-5)
    ax[1][1].set_ylim(-45,-5)
    # plt.savefig('ouioui2.pdf')
    plt.show()

# sampleplot2(Ts=[300,350],pHs=[8,9])


def sampleplot3(Ts=[400, 350,300], pHs=[8,9,10]):
    """For testing different plot types, with nominal line """
    samplesize=1500
    pHnum=100
    fig, ax = plt.subplots(figsize=(10,5), ncols=2, nrows=1)

    # cmaps=[transparent_cmap(238/360), transparent_cmap(120/360), 'Reds']
    cmaps=[transparent_cmap2(blue_to_a), transparent_cmap2(green_to_a), transparent_cmap2(red_to_a)]
    colors = ['blue', 'green','red']

    Tresults = []
    pHresults = []

    pHnoms = np.linspace(7,12,num=int(pHnum))
    Tnoms = np.linspace(273, 400, num=int(pHnum))


    for i, T in enumerate(Ts):
        pHvals, PS = samplepH(samplesize, Temp=T, pHrange=[7,12])#, pHnum=pHnum) #random sample
        PSnoms=[]
        for pHnom in pHnoms:
            PSnoms.append(EG.getPSfromSigmas([0,0,0,1,0,T,pHnom], dummyvals=True))

        zeroes = (np.log10(PS) != -50.)
        y = np.log10(PS)[zeroes]
        x = pHvals[zeroes]
        # ax[0].hexbin(x, y, gridsize=(int(pHnum/2), int(pHnum/8)), cmap=cmaps[i])
        ax[0].hist2d(x,y, bins=[45,45], range=[[7, 12], [-35, -5]], cmap=cmaps[i], edgecolor=None)

        nomzeroes = (np.log10(PSnoms) != -50.)
        nomy = np.log10(PSnoms)[nomzeroes]
        nomx = pHnoms[nomzeroes]
        ax[0].plot(nomx, nomy, c=colors[i], linewidth=2, label='T = '+ str(T))

    for i, pH in enumerate(pHs):
        Tvals, PS = sampleT(samplesize, pH=pH)

        PSnoms=[]
        for Tnom in Tnoms:
            PSnoms.append(EG.getPSfromSigmas([0,0,0,1,0,Tnom,pH], dummyvals=True))

        zeroes = (np.log10(PS) != -50.)
        y = np.log10(PS)[zeroes]
        x = Tvals[zeroes]
        # ax[0].hexbin(x, y, gridsize=(int(pHnum/2), int(pHnum/8)), cmap=cmaps[i])
        ax[1].hist2d(x,y, bins=[45,45], range=[[273, 400], [-35, -5]], cmap=cmaps[i], edgecolor=None)

        nomzeroes = (np.log10(PSnoms) != -50.)
        nomy = np.log10(PSnoms)[nomzeroes]
        nomx = Tnoms[nomzeroes]
        ax[1].plot(nomx, nomy, c=colors[i], linewidth=2, label='pH = '+str(pH))

    for a in ax:
            a.set_ylabel('log10 (Power Supply [W/cell])')
            a.set_ylim(-35,-5)
            # EPS.add_maintnenace_lines(a)
    ax[0].set_xlabel('pH')
    ax[1].set_xlabel('Temperature [K]')
    ax[0].set_xlim(7,12)
    ax[1].set_xlim(273,400)

    ax[0].legend(bbox_to_anchor=(0., 1.05, 1.0, .102), loc='center',
           ncol=3, mode="expand", borderaxespad=0., fontsize=12)
    ax[1].legend(bbox_to_anchor=(0., 1.05, 1.0, .102), loc='center',
           ncol=3, mode="expand", borderaxespad=0., fontsize=12)
    plt.tight_layout()
    plt.savefig('2_binline.pdf')
    # plt.show()

# sampleplot3(Ts=[275,325,400],pHs=[8,9,10])



def total_varianceplot(Ts = [300]):

    samplesize=2500
    pHnum=100
    fig, ax = plt.subplots(figsize=(10,5), ncols=2, nrows=1)

    # cmaps=[transparent_cmap(238/360), transparent_cmap(120/360), 'Reds']
    cmaps=[transparent_cmap2(blue_to_a), transparent_cmap2(green_to_a), transparent_cmap2(red_to_a)]
    colors = ['blue', 'green','red']

    Tresults = []
    pHresults = []
    minpH, maxpH = 7,12
    pHnoms = np.linspace(minpH,maxpH,num=int(pHnum))
    Tnoms = np.linspace(273, 400, num=int(pHnum))


    for i, T in enumerate(Ts):
        pHvals, PS = samplepH(samplesize, Temp=T, pHrange=[7,12])#, pHnum=pHnum) #random sample
        PSnoms=[]
        for pHnom in pHnoms:
            PSnoms.append(EG.getPSfromSigmas([0,0,0,1,0,T,pHnom], dummyvals=True))

        logPSnoms = np.log10(PSnoms)
        zeroes = (np.log10(PS) != -50.)
        logy = np.log10(PS)[zeroes]
        y = PS[zeroes]
        x = pHvals[zeroes]
        logvariance = np.empty(len(logy))
        variance = np.empty(len(y))
        ic = 0
        for xi, yi, logyi in zip(x, y, logy):
            close_nom = int((xi-minpH) / (maxpH-minpH) * pHnum)
            logvariance[ic] = logyi - logPSnoms[close_nom]
            variance[ic] = yi - PSnoms[close_nom]

            ic+=1

            # ax[0].hexbin(x, y, gridsize=(int(pHnum/2), int(pHnum/8)), cmap=cmaps[i])
    ax[0].hist2d(x,logvariance, bins=[45,45], range=[[8, 12], [-5, 5]], cmap=cmaps[i], edgecolor=None, label='T = 300 K')
    ax[0].set_ylabel('log10( variance in power supply [W/cell])')
    ax[1].hist2d(x,variance, bins=[45,45], range=[[8, 12], [-1e-15,1e-15]], cmap=cmaps[i], edgecolor=None, label='T = 300 K')
    ax[1].set_ylabel('Variance in power supply [Wcell]')

    ax[0].legend(bbox_to_anchor=(0., 1.05, 1.0, .102), loc='center',
           ncol=3, mode="expand", borderaxespad=0., fontsize=12)
    ax[1].legend(bbox_to_anchor=(0., 1.05, 1.0, .102), loc='center',
           ncol=3, mode="expand", borderaxespad=0., fontsize=12)
    plt.tight_layout()
    plt.savefig('totvariance_logvslin.pdf')
    # plt.show()

# total_varianceplot(Ts = [300])


def ones_varianceplot(T=300, pH=None, samplesize=100, nominals=True, save=True, show=True):

    xnum=100
    fig, ax = plt.subplots(figsize=(8,12), ncols=2, nrows=3)
    ax = ax.flatten()
    # cmaps=[transparent_cmap(238/360), transparent_cmap(120/360), 'Reds']
    cmaps=[transparent_cmap2(blue_to_a), transparent_cmap2(green_to_a), transparent_cmap2(red_to_a)]
    colors = ['blue', 'green','red']

    Tresults = []
    pHresults = []
    minpH, maxpH = 7,12
    minT, maxT = 273,400
    minx,maxx = 0, 0
    xnoms =[]

    if pH == None:
        minx, maxx = minpH, maxpH
        xnoms = np.linspace(minpH,maxpH,num=int(xnum))
    elif T== None:
        minx, maxx = minT, maxT
        xnoms = np.linspace(minT, maxT, num=int(xnum))
    if nominals:
        PSnoms=[]
        for xnom in xnoms:
            if pH == None:
                PSnoms.append(EG.getPSfromSigmas([0,0,0,1,0,T,xnom], dummyvals=True))
            elif T == None:
                PSnoms.append(EG.getPSfromSigmas([0,0,0,1,0,xnom,pH], dummyvals=True))

        logPSnoms = np.log10(PSnoms)

    ind = ['all', 'CH4', 'CO2','H2','nATP', 'k_corr']
    for i, v in enumerate(ind):
        xvals, PS = [],[]
        if pH == None:
            xvals, PS = samplepH_onlyone(samplesize, one=v, Temp=T, pHrange=[minx,maxx])#, pHnum=pHnum) #random sample
            ax[i].set_xlabel('pH')
        elif T == None:
            xvals, PS = sampleT_onlyone(samplesize, one=v, Trange=[minx,maxx], pH=pH)#, pHnum=pHnum) #random sample
            ax[i].set_xlabel('Temperature [K]')

        if save or show:
            _logy = np.log10(PS)
            zeroes = (_logy != -50.)
            logy = _logy[zeroes]
            # y = PS[zeroes]
            x = xvals[zeroes]
            logvariance = np.empty(len(logy))
            # variance = np.empty(len(y))
            ic = 0
            for xi, logyi in zip(x, logy):
                close_nom = int((xi-minx) / (maxx-minx) * xnum)
                logvariance[ic] = logyi - logPSnoms[close_nom]
                ic+=1
            ax[i].hist2d(x,logvariance, bins=[45,45], range=[[minx, maxx], [-4, 4]], cmap=cmaps[0], edgecolor=None)

            ax[i].set_title(v)
            ax[i].set_ylabel('Variance in log10 (power supply)')

    if pH == None:
        plt.suptitle('Variance with pH when T = '+str(T))
    elif T == None:
        plt.suptitle('Variance with T when pH = '+str(pH))
    plt.tight_layout()
    if save:
        plt.savefig('2_totvariance_independents_T_'+str(T)+'_pH_'+str(pH)+'.pdf')
    if show:
        plt.show()
    plt.close()


# ones_varianceplot(T=300, pH=None)




#
# for T in [275,325, 400]:
#     # variance_chain_sample(T=T, pH=None, stepsize=100, totsize=1500, startfrom=500)
#     ones_varianceplot(T=T, pH=None, samplesize=500, nominals=False, save=False, show=False)
# for pH in [8,9,10]:
#     ones_varianceplot(T=None, pH=pH, samplesize=500, nominals=False, save=False, show=False)
    # variance_chain_sample(T=None, pH=pH, stepsize=100, totsize=500, startfrom=500)
#
# for T in [275,325, 400]:
#     # variance_chain_sample(T=T, pH=None, stepsize=100, totsize=1500, startfrom=500)
#     ones_varianceplot(T=T, pH=None, samplesize=1500, nominals=True, save=True, show=False)
# for pH in [8,9,10]:
#     ones_varianceplot(T=None, pH=pH, samplesize=1500, nominals=True, save=True, show=False)
#     # variance_chain_sample(T=None, pH=pH, stepsize=100, totsize=1500, startfrom=500)

# variance_chain_sample(None, 8, stepsize=50, totsize=650, startfrom=500)

# ones_varianceplot(T=None, pH=8, samplesize=650, nominals=True, save=True, show=True)
