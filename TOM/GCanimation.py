import methanogen_extractor as extractor
from methanogen_implementer import efficiencies
import unique_efficiencies as ueff
import sys,os, statistics, math, ast
import numpy as np
import pandas as pd
import sqlite3

sys.path.append(os.path.dirname(__file__)+'/../../../')
import NutMEG as es
import NutMEG.plotter as nutplt
import NutMEG.util.NutMEGparams as nmp

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation

this_dbpath = 'amparams_db'
es.db_helper.create_major_tables(replace=False, dbpath=this_dbpath)

""" THIS CODE IS FOR GENERATING THE SUPPLEMENTAL ANIMATION.

THERE IS SOME OVERLAP WITH EnergyNutrientLimitation.py
"""

def get_peffs(other_params={}, MPscale=1, T=300):
    """Return a list of efficincies objects for the typical optimal methanogen,
    each having different nATP yield (0.5, 1.0, 1.5 per mol CO2). Change individual
    parameters for analysis by passing other_params.
    """
    E, PT, PS, PG, S = extractor.extract_Esynths_csv('data/TOM_PT/3e-08_'+str(T))
    E05, PT05, PS05, PG05, S05 = extractor.extract_Esynths_csv('data/TOM_PT/05_3e-08_'+str(T))
    E15, PT15, PS15, PG15, S15 = extractor.extract_Esynths_csv('data/TOM_PT/15_3e-08_'+str(T))


    p_std = {'Tdef':'None',
      'n_ATP':1.0, 'MP':PT[0]*MPscale, 'Temp':T}
    p_std.update(other_params)
    peff1 = efficiencies.get_eff('averageMethanogen', Temp=T,
      paramchange=p_std)

    p_std05 = {'Tdef':'None',
      'n_ATP':0.5, 'MP':PT05[0]*MPscale, 'Temp':T}
    p_std05.update(other_params)
    peff05 = efficiencies.get_eff('averageMethanogen', Temp=T,
      paramchange=p_std05)

    p_std15 = {'Tdef':'None',
      'n_ATP':1.5, 'MP':PT15[0]*MPscale, 'Temp':T}
    p_std15.update(other_params)
    peff15 = efficiencies.get_eff('averageMethanogen', Temp=T,
      paramchange=p_std15)



    return [peff05, peff1, peff15]




def inbetweens(peffs, ax, rmv=False, stoppers={}, T=300, dt=2000, hatch=None, alph=0.6):
    """Plot growth curves for the three eddiciencies objects in peffs, with a
    fill_between plot between indexes 0 and 2, and a straight line at index 1
    ATP yield of 1 mol per mol CO2.
    """
    collist=['#2c7bb6','#abd9e9', '#fdae61', '#d7191c', 'm']
    for ic, peff in enumerate(peffs):
        L,O,S,t,pop= [0,0,0], [0,0,0],[0,0,0],[0,0,0],[0,0,0]
        for i in range(len(peff)):
            L[i], O[i], S[i] = peff[i].perform_sim(dbpath=this_dbpath, dt=dt, stoppers=stoppers)
            t2,pop2=[],[]
            if rmv:
                es.ecosystem_dbhelper.db_helper.removesimdata(S[i], dbpath=this_dbpath)
                L[i], O[i], S[i] = peff[i].perform_sim(dbpath=this_dbpath, dt=dt, stoppers=stoppers)
            t2 = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S[i], 'Time', dbpath=this_dbpath)
            pop2 = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S[i], 'no_alive_'+O[i], dbpath=this_dbpath)
            t[i], pop[i]=[],[]
            for k, p in zip(t2, pop2):
                t[i].append(k[0])
                pop[i].append(p[0])

        if t[0][-1]<t[2][-1]:
            for j in range(int((t[2][-1]-t[0][-1])/dt)):
                t[0].append(t[0][-1]+dt)
                pop[0].append(pop[0][-1])
        else:
            for j in range(int((t[0][-1]-t[2][-1])/dt)):
                t[2].append(t[2][-1]+dt)
                pop[2].append(pop[2][-1])
        try:
            if hatch!=None:
                ax.plot(t[1],pop[1], c=collist[ic], ls='--', alpha=alph)
                ax.fill_between(t[0], pop[0], pop[2], alpha=alph, edgecolor=collist[ic], facecolor='none', hatch=hatch)
            else:
                ax.plot(t[1],pop[1], c=collist[ic], alpha=alph)
                ax.fill_between(t[0], pop[0], pop[2], alpha=alph, color=collist[ic])
        except ValueError:
            inbetweens(peffs, ax, rmv=True, stoppers=stoppers, T=T, dt=dt, alph=alph)


def boxes(peffs, ax, rmv=False, stoppers={}, data='', plot=False, boxlst=[[],[]], T=300, dt=2000):
    collist=['#abd9e9', '#fdae61', '#d7191c', 'm']
    if len(boxlst[0])==0:
        collist=['#2c7bb6','#abd9e9', '#fdae61', '#d7191c', 'm']

    for ip,peff in enumerate(peffs):
        L,O,S,t,GR= [0,0,0], [0,0,0],[0,0,0],[0,0,0],[0,0,0]
        for i in range(len(peff)):
            L[i], O[i], S[i] = peff[i].perform_sim(dbpath=this_dbpath, dt=dt, stoppers=stoppers)
            t2,GR2=[],[]
            if rmv:
                es.ecosystem_dbhelper.db_helper.removesimdata(S[i], dbpath=this_dbpath)
                L[i], O[i], S[i] = peff[i].perform_sim(dbpath=this_dbpath, dt=dt, stoppers=stoppers)
            t2 = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S[i], 'Time', dbpath=this_dbpath)
            GR2 = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S[i], data, OrgID='', dbpath=this_dbpath)
            if data=='Composition':
                # print(GR2)
                GR2 = ast.literal_eval(tuple(GR2)[-1][0])['CH4(g)']
                # print(GR2)
            elif type(GR2)==type(''):
                GR2 = ast.literal_eval(GR2)[0]
                GR2 = GR2*3600
            elif type(GR2)==type([]):
                GR2=GR2[-1]
            t[i], GR[i]=[],GR2
            for k in t2:
                t[i].append(k[0])

        boxlst[0].append({
            'label' : '',
            'whislo': GR[0],    # Bottom whisker position
            'q1'    : GR[0],    # First quartile (25th percentile)
            'med'   : GR[1],    # Median         (50th percentile)
            'q3'    : GR[2],    # Third quartile (75th percentile)
            'whishi': GR[2],    # Top whisker position
            'fliers': []        # Outliers
        })
        boxlst[1].append(collist[ip])
    if not plot:
        return boxlst
    else:
        bplot = ax.bxp(boxlst[0][::-1], showfliers=False, vert=False, patch_artist=True, widths=[0.9]*len(boxlst[0][::-1]))
        for ib, patch in enumerate(bplot['boxes']):
            patch.set_facecolor(boxlst[1][::-1][ib])
            patch.set_alpha(0.6)
            patch.set_edgecolor(boxlst[1][::-1][ib])
        for ib, line in enumerate(bplot['medians']):
            line.set_color(boxlst[1][::-1][ib])
        if data=='PeakGR':
            ax.set_xlabel('Peak Growth rate [$\mathregular{hr}^{-1}$]')
        elif data=='FinBM_cells_tot':
            ax.set_xlabel('Final Biomass [cells $\mathregular{L}^{-1}$]')
        elif data=='Composition':
            ax.set_xlabel('Total CH4 produced [mol $\mathregular{L}^{-1}$]')



def maintenancepeffs(T=300):
    """return peffs with fractionally increased maintenance costs."""
    return [
      get_peffs(MPscale=1, T=T),
      get_peffs(MPscale=1.01, T=T),
      get_peffs(MPscale=1.1, T=T),
      get_peffs(MPscale=1.15, T=T)
      ]

def CO2peffs(T=300):
    """return peffs with fractionally reduced CO2 concentrations"""
    peffs = [get_peffs(T=T), get_peffs(T=T), get_peffs(T=T), get_peffs(T=T), get_peffs(T=T)]
    for peff in peffs[1]:
        peff.params['mol_CO2'] = peff.params['mol_CO2']*0.99
    for peff in peffs[2]:
        peff.params['mol_CO2'] = peff.params['mol_CO2']*0.9
    for peff in peffs[3]:
        peff.params['mol_CO2'] = peff.params['mol_CO2']*0.8
    for peff in peffs[4]:
        peff.params['mol_CO2'] = peff.params['mol_CO2']*0.6

    return peffs

def H2peffs(T=300):
    """return peffs with fractionally reduced H2 concentrations"""
    peffs = [get_peffs(T=T), get_peffs(T=T), get_peffs(T=T), get_peffs(T=T)]

    for peff in peffs[1]:
        peff.params['mol_H2'] = peff.params['mol_H2']*0.99
    for peff in peffs[2]:
        peff.params['mol_H2'] = peff.params['mol_H2']*0.95
    for peff in peffs[3]:
        peff.params['mol_H2'] = peff.params['mol_H2']*0.8
    return peffs

def Ppeffs(T=300):
    """return peffs with limiting P concentration"""
    return [get_peffs(T=T),
      get_peffs(T=T, other_params={'Pconc':1e-8}),
      get_peffs(T=T, other_params={'Pconc':1e-10})]

def CO2inpeffs(T=300):
    """return peffs with an outflow of CO2"""
    peffs = [get_peffs(T=T), get_peffs(T=T), get_peffs(T=T)]
    for peff in peffs[1]:
        peff.params['inputs'] = {'CO2(aq)':-1e-12}
    for peff in peffs[2]:
        peff.params['inputs'] = {'CO2(aq)':-1e-8}
    return peffs

def H2inpeffs(T=300):
    """return peffs with an outflow of H2"""
    peffs = [get_peffs(T=T), get_peffs(T=T), get_peffs(T=T)]
    for peff in peffs[1]:
        peff.params['inputs'] = {'H2(aq)':-1e-12}
    for peff in peffs[2]:
        peff.params['inputs'] = {'H2(aq)':-1e-8}
    return peffs

def lifespanpeffs(T=300):
    """return peffs with a lifespan"""
    peffs = [get_peffs(T=T), get_peffs(T=T), get_peffs(T=T)]

    for peff in peffs[1]:
        peff.params['lifespan'] = 2e5
    for peff in peffs[2]:
        peff.params['lifespan'] = 1e5
    return peffs

def Puptakepeffs(T=300):
    """return peffs with an altered P rate constant"""
    return [get_peffs(T=T),
      get_peffs(T=T, other_params={'uptake_consts':{'P':1e-18}, 'Pconc':0.10008}),
      get_peffs(T=T, other_params={'uptake_consts':{'P':1e-20}, 'Pconc':0.10011})]

def plot_boxes(axs, lineaxs, T, dt):

    lablist=[]
    GRboxs =[[],[]]
    TBboxs=[[],[]]
    Cboxs=[[],[]]

    plens =[]

    axs[1].set_xscale('linear')
    axs[0].set_xscale('linear')
    axs[2].set_xscale('linear')

    ############################## MAINTENANCE
    peffs = maintenancepeffs(T=T)
    plens.append(len(peffs[1:]))
    lablist.extend(['+0.1%', '+10%', '+15%'])
    GRboxs=boxes(peffs, axs[0], data='PeakGR', boxlst=GRboxs, T=T, dt=dt)
    TBboxs=boxes(peffs, axs[1], data='FinBM_cells_tot', boxlst=TBboxs, T=T, dt=dt)
    CBboxs=boxes(peffs, axs[2], data='Composition', boxlst=Cboxs, T=T, dt=dt)

    ############################## CO2 CONCENTRATION
    peffs = CO2peffs(T=T)
    plens.append(len(peffs[1:]))
    lablist.extend(['-1%', '-10%','-20%','-40%'])
    GRboxs=boxes(peffs[1:], axs[0], data='PeakGR', boxlst=GRboxs, T=T, dt=dt)
    TBboxs=boxes(peffs[1:], axs[1], data='FinBM_cells_tot', boxlst=TBboxs, T=T, dt=dt)
    CBboxs=boxes(peffs[1:], axs[2], data='Composition', boxlst=Cboxs, T=T, dt=dt)

    ############################## H2 CONCENTRATION
    peffs = H2peffs(T=T)
    plens.append(len(peffs[1:]))
    lablist.extend(['-1%','-5%','-20%'])
    GRboxs=boxes(peffs[1:], axs[0], data='PeakGR', boxlst=GRboxs)
    TBboxs=boxes(peffs[1:], axs[1], data='FinBM_cells_tot', boxlst=TBboxs)
    CBboxs=boxes(peffs[1:], axs[2], data='Composition', boxlst=Cboxs, T=T, dt=dt)

    ############################## P CONCENTRATION
    peffs = Ppeffs(T=T)
    plens.append(len(peffs[1:]))
    lablist.extend(['$10^{-8}$ M', '$10^{-10}$ M'])
    GRboxs=boxes(peffs[1:], axs[0], data='PeakGR', boxlst=GRboxs, T=T, dt=dt)
    TBboxs=boxes(peffs[1:], axs[1], data='FinBM_cells_tot', boxlst=TBboxs, T=T, dt=dt)
    CBboxs=boxes(peffs[1:], axs[2], data='Composition', boxlst=Cboxs, T=T, dt=dt)

    ############################## CO2 inflow
    peffs=CO2inpeffs(T=T)
    plens.append(len(peffs[1:]))
    lablist.extend(['$10^{-12}$ M/s', '$10^{-8}$ M/s'])
    GRboxs=boxes(peffs[1:], axs[0], data='PeakGR', boxlst=GRboxs, T=T, dt=dt)
    TBboxs=boxes(peffs[1:], axs[1], data='FinBM_cells_tot', boxlst=TBboxs, T=T, dt=dt)
    CBboxs=boxes(peffs[1:], axs[2], data='Composition', boxlst=Cboxs, T=T, dt=dt)

    ############################## H2 inflow
    peffs = H2inpeffs(T=T)
    plens.append(len(peffs[1:]))
    lablist.extend(['$10^{-12}$ M/s', '$10^{-8}$ M/s'])
    GRboxs=boxes(peffs[1:], axs[0], data='PeakGR', boxlst=GRboxs, T=T, dt=dt)
    TBboxs=boxes(peffs[1:], axs[1], data='FinBM_cells_tot', boxlst=TBboxs, T=T, dt=dt)
    CBboxs=boxes(peffs[1:], axs[2], data='Composition', boxlst=Cboxs, T=T, dt=dt)

    ############################## LIFESPAN
    peffs = lifespanpeffs(T=T)
    plens.append(len(peffs[1:]))
    lablist.extend(['200000 s','100000 s'])
    GRboxs=boxes(peffs[1:], axs[0], data='PeakGR', boxlst=GRboxs, T=T, dt=dt)
    TBboxs=boxes(peffs[1:], axs[1], data='FinBM_cells_tot', boxlst=TBboxs, T=T, dt=dt)
    CBboxs=boxes(peffs[1:], axs[2], data='Composition', boxlst=Cboxs, T=T, dt=dt)

    ############################## P UPTAKE
    peffs = Puptakepeffs(T=T)
    plens.append(len(peffs[1:]))
    lablist.extend(['$10^{-13}$ /s', '$10^{-15}$ /s' ])

    GRboxs=boxes(peffs[1:], axs[0], data='PeakGR', boxlst=GRboxs, plot=True, T=T, dt=dt)
    TBboxs=boxes(peffs[1:], axs[1], data='FinBM_cells_tot', boxlst=TBboxs, plot=True, T=T, dt=dt)
    CBboxs=boxes(peffs[1:], axs[2], data='Composition', boxlst=Cboxs, plot=True, T=T, dt=dt)

    for ax in axs[:-1]:
        ax.set_title('Temp: '+str(T)+ 'K')
        ax.set_xscale('log')


    # lines and text
    for lax in lineaxs:
        lax.set_ylim(0.5,sum(plens)+1.5)
    axs[3].set_ylim(0.5,sum(plens)+1.5)
    axs[3].set_xlim(0,0.5)
    tottt=sum(plens)+0.5
    axs[3].text(0.02,tottt+0.1,'OPTIMAL')
    axs[3].text(0.02,tottt+0.1,'OPTIMAL')
    for lax in lineaxs:
        lax.axhline(tottt, c='k', lw=2)
    axs[3].axhline(tottt, c='k', lw=2)
    names=['Maintenance', '[$\mathregular{CO}_{2}$]','[$\mathregular{H}_{2}$]',
      '[P]', '$\mathregular{CO}_{2}$ sink', '$\mathregular{H}_{2}$ sink',
      'Lifespan','RateConst P uptake']
    labcount=0
    for p,s in zip(plens, names):
        axs[3].text(0.02,tottt-0.9,s)
        axs[3].text(0.02,tottt-0.9,s)

        for testo in  range(0,p):
            axs[3].text(0.49,sum(plens)-0.4-labcount,lablist[labcount],horizontalalignment='right')
            labcount+=1

        tottt -= p
        for lax in lineaxs:
            lax.axhline(tottt, c='k', lw=2)
        axs[3].axhline(tottt, c='k', lw=2)

    for lax in lineaxs:
        lax.tick_params(colors='none', axis='y', length=0)

    axs[3].tick_params(colors='none', axis='y', zorder=1, length=0)
    axs[3].tick_params(colors='none', axis='x', length=0)

    axs[1].tick_params(axis='y', right=False, left=False)
    axs[0].tick_params(axis='y', right=False, left=False)

    return axs

def plot_growthcurves(axs, T, dt, rmv=False):

     ############################## MAINTENANCE
    peffs = maintenancepeffs(T=T)
    inbetweens(peffs, axs[0], T=T, dt=dt, rmv=rmv)
    axs[0].set_title('Maintenance')

    ############################## CO2 CONCENTRATION
    peffs = CO2peffs(T=T)
    inbetweens(peffs, axs[1], T=T, dt=dt,rmv=rmv)
    axs[1].set_title('[$\mathregular{CO}_{2}$]')

    ############################## H2 CONCENTRATION
    peffs = H2peffs(T=T)
    inbetweens(peffs, axs[2], T=T, dt=dt, rmv=rmv)
    axs[2].set_title('[$\mathregular{H}_{2}$]')

    ############################## P CONCENTRATION
    peffs = Ppeffs(T=T)
    inbetweens(peffs, axs[3], T=T, dt=dt, rmv=rmv)
    axs[3].set_title('[P]')

    ############################## CO2 inflow
    peffs=CO2inpeffs(T=T)
    inbetweens(peffs, axs[4], rmv=rmv, stoppers={'Maintenance_Fraction':{'Max':1.2, 'Min':-0.1, 'Consistency':10, 'Count':0}}, T=T, dt=dt)
    axs[4].set_title('$\mathregular{CO}_{2}$ Sink')

    ############################## H2 inflow
    peffs = H2inpeffs(T=T)
    inbetweens(peffs, axs[5], rmv=rmv, stoppers={'Maintenance_Fraction':{'Max':1.2, 'Min':-0.1, 'Consistency':10, 'Count':0}}, T=T, dt=dt)
    axs[5].set_title('$\mathregular{H}_{2}$ Sink')

    ############################## LIFESPAN
    peffs = lifespanpeffs(T=T)
    inbetweens(peffs, axs[6], rmv=rmv, T=T, dt=dt)
    axs[6].set_title('Lifespan')



    ############################## P UPTAKE
    peffs = Puptakepeffs(T=T)
    inbetweens(peffs, axs[7], rmv=rmv, T=T, dt=dt)
    axs[7].set_title('$k_\mathregular{P}$')


    for ax in [axs]:
        ax[0].set_ylabel('Biomass [cells $\mathregular{L}^{-1}$]')
        for a in ax:
            a.set_yscale('log')
            a.set_xscale('log')
            a.set_xlim(1e2, 5e7)
            a.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.8)

            a.set_xlabel('Time [s]')

def doubleanim():
    fig = plt.figure(constrained_layout=True, figsize=(16,8))
    gs = fig.add_gridspec(3,1)

    gcgs = gs[0].subgridspec(1,8)
    bgs = gs[1:].subgridspec(1,4)

    gcaxs=[]
    baraxs=[]
    for i in range(8):
        gcaxs.append(fig.add_subplot(gcgs[i]))
    for i in range(4):
        baraxs.append(fig.add_subplot(bgs[i]))

    lineaxs = [baraxs[1].twinx(),baraxs[2].twinx(), baraxs[3].twinx()]

    def doublean(num):
        T = 280+(num)
        dt =3000-(num*(1000/25))
        if T==300:
            dt=2000
        if T==325:
            dt=1000
        rmv=False


        for a in gcaxs:
            a.clear()

        for ax in baraxs:
            ax.set_xscale('linear')
            ax.clear()
        for ax in lineaxs:
            ax.clear()


        plot_growthcurves(gcaxs, int(T), int(dt), rmv=rmv)
        for i,a in enumerate(gcaxs):
            a.set_ylim(100,1e12)
            if i != 0:
                a.tick_params(labelcolor='none', axis='y')

        plot_boxes([baraxs[1],baraxs[2], baraxs[3], baraxs[0]], lineaxs, int(T), int(dt))
        baraxs[1].set_xlim(1e-4,15*3600e-5)
        baraxs[2].set_xlim(1e3,1e12)
        baraxs[3].set_xlim(2e-8,1e-3)


        plt.tight_layout()
        fig.subplots_adjust(wspace=0)

        return fig

    T=250
    dt=4000

    interval=200

    fig.dpi=80

    frames=range(51)

    ani = animation.FuncAnimation(fig, doublean, frames=len(frames), interval=interval)
    ani.save('GCan.gif', writer='imagemagick')


# don't remove any simulations
rmv=False

doubleanim()
