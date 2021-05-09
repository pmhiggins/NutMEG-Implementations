import sys, os, ast, math
sys.path.append(os.path.dirname(__file__)+'../../NutMEG')
sys.path.append(os.path.dirname(__file__)+'../../NutMEG-Implementations/TOM')
from itertools import chain
import matplotlib.pyplot as plt


import Sampler
import SamplingPlotStyles as SPS


ind = ['CH4', 'H2' , 'nATP', 'k_corr']
salt = ['nom', 'high', 'low']
Temps = [275,325,375]
pHs = [8,9,10]

sample = 1000
# for doing samples

"""
for one in ind:
    for s in salt:
        for T in Temps:
            Sampler.independent_sample(T=T, pH=None, to_sample=one,
              samplesize=sample,
              Trange = [273,400], pHrange=[7,12], salttype=s)
        for pH in pHs:
            Sampler.independent_sample(T=None, pH=pH, to_sample=one,
              samplesize=sample,
              Trange = [273,400], pHrange=[7,12], salttype=s)
        print(one, s)
"""

# plotting variance
"""
for s in salt:
    for T in Temps:
        SPS.ones_varianceplot(T=T, pH=None, salttype=s, samplesize=sample, Trange=[273,400], pHrange=[7,12], save=True, show=False, cm='Blues', fixupdate={}, fn_preamble='')
    for pH in pHs:
        SPS.ones_varianceplot(T=None, pH=pH, salttype=s, samplesize=sample, Trange=[273,400], pHrange=[7,12], save=True, show=False, cm='Blues', fixupdate={}, fn_preamble='')
"""

def varianceexample(samplesize):
    """ Plot overall variance at select pH and salt levels """
    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
    ax[0][0], hb, hr = SPS.all_varianceplot_T(ax[0][0], samplesize, 8, salttype='high')
    ax[0][1], hb, hr = SPS.all_varianceplot_T(ax[0][1], samplesize, 9, salttype='high')
    ax[0][2], hb, hr = SPS.all_varianceplot_T(ax[0][2], samplesize, 10, salttype='high')

    ax[1][0], hb, hr = SPS.all_varianceplot_T(ax[1][0], samplesize, 8, salttype='nom')
    ax[1][1], hb, hr = SPS.all_varianceplot_T(ax[1][1], samplesize, 9, salttype='nom')
    ax[1][2], hb, hr = SPS.all_varianceplot_T(ax[1][2], samplesize, 10, salttype='nom')

    ax[2][0], hb, hr = SPS.all_varianceplot_T(ax[2][0], samplesize, 8, salttype='low')
    ax[2][1], hb, hr = SPS.all_varianceplot_T(ax[2][1], samplesize, 9, salttype='low')
    ax[2][2], hb, hr = SPS.all_varianceplot_T(ax[2][2], samplesize, 10, salttype='low')

    ax[0][0].annotate('Bulk ocean pH: 8', xy=(0.5, 1.1), xytext=(0.5, 1.2), xycoords='axes fraction',
            fontsize=12, ha='center', va='center', rotation='horizontal',
            arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=2))

    ax[0][1].annotate('Bulk ocean pH: 9', xy=(0.5, 1.1), xytext=(0.5, 1.2), xycoords='axes fraction',
            fontsize=12, ha='center', va='center', rotation='horizontal',
            arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=2))

    ax[0][2].annotate('Bulk ocean pH: 10', xy=(0.5, 1.1), xytext=(0.5, 1.2), xycoords='axes fraction',
            fontsize=12, ha='center', va='center', rotation='horizontal',
            arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=2))

    ax[0][-1].annotate('high salt ocean', xy=(1.1, 0.5), xytext=(1.2, 0.5), xycoords='axes fraction',
            fontsize=12, ha='left', va='center', rotation='vertical',
            arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=2))

    ax[1][-1].annotate('nominal ocean', xy=(1.1, 0.5), xytext=(1.2, 0.5), xycoords='axes fraction',
            fontsize=12, ha='left', va='center', rotation='vertical',
            arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=2))

    ax[2][-1].annotate('low salt ocean', xy=(1.1, 0.5), xytext=(1.2, 0.5), xycoords='axes fraction',
            fontsize=12, ha='left', va='center', rotation='vertical',
            arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=2))

    cbaxes = fig.add_axes([0.1, 0.06, 0.35, 0.02])
    cbaxes2 = fig.add_axes([0.55, 0.06, 0.35, 0.02])

    plt.colorbar(hb[-1], cax=cbaxes, label='no. sample power supplies in bin', orientation='horizontal', pad=0.5, aspect=7, extend='max')

    # cmappable = ScalarMappable(norm=Normalize(0,1), cmap=cmapper.r2a())
    fig.colorbar(hr[-1], cax=cbaxes2, label='no. sample power supplies in bin (when nominal = 0)', orientation='horizontal', pad=0.5, aspect=7, extend='max')

    for a in chain(ax[0], ax[1]):
        a.get_xaxis().set_ticks([])

    for i in ax:
        for a in i[1:]:
            a.get_yaxis().set_ticks([])

    for a in ax:
        # LHS
        a[0].set_ylabel('Variance')
    for a in ax[-1]:
        # bottom row
        a.set_xlabel('Temperature [K]')

    # ax[1][0] = all_varianceplot_pH(ax[1][0], 100, 300)
    # ax[1][1] = all_varianceplot_pH(ax[1][1], 100, 300, salttype='high')
    fig.subplots_adjust(bottom=0.15, left=0.08, right=0.92, top=0.92, wspace=0.05, hspace=0.05)
    plt.savefig('figs/total_variance.pdf')
    # plt.show()

def varianceexmaple_extend(startfrom, totsize, stepsize):

    Sampler.variance_chain_sample(None, 8,
      stepsize=stepsize, totsize=totsize, startfrom=startfrom,
      salttype='nom',
      ones=['all'])
    Sampler.variance_chain_sample(None, 9,
      stepsize=stepsize, totsize=totsize, startfrom=startfrom,
      salttype='nom',
      ones=['all'])
    Sampler.variance_chain_sample(None, 10,
      stepsize=stepsize, totsize=totsize, startfrom=startfrom,
      salttype='nom',
      ones=['all'])

    Sampler.variance_chain_sample(None, 8,
      stepsize=stepsize, totsize=totsize, startfrom=startfrom,
      salttype='high',
      ones=['all'])
    Sampler.variance_chain_sample(None, 9,
      stepsize=stepsize, totsize=totsize, startfrom=startfrom,
      salttype='high',
      ones=['all'])
    Sampler.variance_chain_sample(None, 10,
      stepsize=stepsize, totsize=totsize, startfrom=startfrom,
      salttype='high',
      ones=['all'])

    Sampler.variance_chain_sample(None, 8,
      stepsize=stepsize, totsize=totsize, startfrom=startfrom,
      salttype='low',
      ones=['all'])
    Sampler.variance_chain_sample(None, 9,
      stepsize=stepsize, totsize=totsize, startfrom=startfrom,
      salttype='low',
      ones=['all'])
    Sampler.variance_chain_sample(None, 10,
      stepsize=stepsize, totsize=totsize, startfrom=startfrom,
      salttype='low',
      ones=['all'])
