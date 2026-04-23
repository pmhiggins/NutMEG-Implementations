import sys, os
# sys.path.append(os.path.dirname(__file__)+'/../../../../../NutMEG')

sys.path.append(os.path.dirname(__file__)+'/../EncBmBs_utils')
import PlotSetup

from ChiralTools import ChiralTools as CT
import matplotlib.pyplot as plt
import numpy as np
from itertools import chain
from scipy.interpolate import griddata
import pandas as pd

class ChiralPlotter:

    def ChiralChange_with_abio_component(abioL=50, T=273):

        biosig_Ls = np.linspace(50,100, num=51)

        fig, axs = plt.subplots(figsize=(5,10), nrows=3, ncols=1)
        cmap = plt.cm.get_cmap('viridis')
        axs = axs.flatten()

        sAB = np.where(biosig_Ls > abioL)[0][0]
        _abio_Ls_50 = (biosig_Ls - 100 )/ (50 - biosig_Ls)

        ks = ['kglu','kLever']
        secmin = CT.yr_to_s(10**-1.)
        secmax = CT.yr_to_s(10**8.)
        ts = [10**j for j in np.linspace(np.log10(secmin),np.log10(secmax),num=500)]
        titles = ['Slowest racemization at '+str(T-273)+r'$\degree$C', 'Fastest racemization at '+str(T-273)+r'$\degree$C']
        kts = [[],[]]


        for ax, ktype, ti in zip(axs[1:], ks, titles):

            x = []
            y  = []
            z= []
            for i, t in enumerate(ts):
                ratio =[]
                k = CT.get_k(T, ktype)
                t_thres=False
                for bs in biosig_Ls[sAB:]:
                    sourceDL = CT.get_initial_DL(k, t, CT.Lf_to_DL(bs))
                    sourceL = CT.DL_to_Lf(sourceDL)
                    ABIOtoBIO = None
                    if sourceDL == -1:
                        sourceL=100
                        ABIOtoBIO=0.
                    if sourceL==50 or sourceL==float('nan'):
                        ABIOtoBIO = 1.
                    else:
                        ABIOtoBIO = (sourceL - 100 )/ (50 - sourceL)
                    _ratio = 1/ (ABIOtoBIO+1)
                    ratio.append(_ratio)

                ax.plot(ratio, biosig_Ls[sAB:], c=cmap(i/(len(ts)-1)), label=round(np.log10(CT.s_to_yr(t)),1))

            ax.set_title(ti)

            ax.plot(1/(_abio_Ls_50[sAB:]+1), biosig_Ls[sAB:], label='Chirality at source', c='k')

            ax.set_xlim(0.,1.)
            ax.set_ylim(50,100)


            ax.set_xlabel('Fraction of amino acids that are biological at source')
            ax.set_ylabel(r'Observed chirality ($L$-form %)')


        norm=plt.Normalize(vmin=np.log10(CT.s_to_yr(min(ts))), vmax=np.log10(CT.s_to_yr(max(ts))))

        DL_abio = 1.
        f_bios=np.linspace(0,1,101)
        glu=[]
        L_thress = [51, 60]
        lss = ['-', 'dashed']
        ktypes = ['kglu', 'kLever']
        css = [['#e66101', '#5e3c99'],['#fdb863', '#b2abd2']]
        labels = [['51%; slowest racemization',
          '51%; fastest racemization'],
          ['60%; slowest racemization',
          '60%; fastest racemization']]


        for L_thres, ls, cs, lbls in zip(L_thress, lss, css, labels):
            for k, c, lbl in zip(ktypes, cs, lbls):
                t_fs = []
                for i, fb in enumerate(f_bios):
                    _DL = CT.DL_from_DLabio_fbio(DL_abio, f_bio=fb)
                    tr_g, Ratio_g, Mean_g = CT.time_to_racemic(T, k, steps=1e5, init_DL=_DL, early_cut=True)
                    t_f = tr_g[np.argmax(CT.DL_to_Lf(Ratio_g)<=L_thres)]
                    t_fs.append(t_f)

                yvs = CT.s_to_yr(np.array(t_fs))
                minS = yvs[0]
                for i, yv in enumerate(yvs):
                    if yv == minS:
                        yvs[i]=1e-50

                axs[0].plot(f_bios, yvs, c=c, ls=ls, label=lbl)

        axs[0].set_yscale('log')
        axs[0].set_ylim(1e-1, 1e9)
        axs[0].set_xlabel('Fraction of amino acids that are biological at source')
        axs[0].set_ylabel('Time until L-form threshold met [yr]')
        axs[0].set_title(r'Time at '+str(T-273)+'$\degree$C until false negative risk')
        axs[0].legend(fontsize=10)


        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        plt.subplots_adjust(left=0.15, bottom=0.14, top=0.96, right=0.9, hspace=0.4)
        cax = plt.axes([0.15, 0.05, 0.75, 0.03])

        fig.colorbar(sm, cax=cax, orientation='horizontal', label=r'$\log_{10}$ (Time [yr])')
        plt.savefig('figs/ChiralChange_with_abio_component_'+str(T)+'.pdf')
        plt.savefig('figs/ChiralChange_with_abio_component_'+str(T)+'.svg')

        plt.close()

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['dejavusans']
plt.rcParams['mathtext.fontset'] = 'dejavusans'
ChiralPlotter.ChiralChange_with_abio_component(T=273)
ChiralPlotter.ChiralChange_with_abio_component(T=333)
