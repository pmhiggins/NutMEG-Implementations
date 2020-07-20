import methanogen_extractor as extractor
from methanogen_implementer import efficiencies
import unique_efficiencies as ueff
import sys,os, statistics, math, ast
import numpy as np
import pandas as pd
import sqlite3

sys.path.append(os.path.dirname(__file__)+'/../../NutMEG/')
import NutMEG as es
import NutMEG.plotter as nutplt
import NutMEG.util.NutMEGparams as nmp

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation

"""
THIS CODE WAS USED FOR FIGURES 3 AND 4: GROWTH CURVES AND TEMPERATURE CURVES

"""

this_dbpath = 'amparams_db' # path to the NutMEG database.

es.db_helper.create_major_tables(replace=False, dbpath=this_dbpath)

# Set some global parameters

# use this boolean to toggle whether simulations should be deleted and rerun
# false will use the data we already have if available.
rmv=False

# change matplotlib rcparams updating plot fonts etc.
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cmr10'
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset'] = 'cm'



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
    collist=['#2c7bb6','#abd9e9', '#fdae61', '#d7191c']
    for ic, peff in enumerate(peffs):
        L,O,S,t,pop= [0,0,0], [0,0,0],[0,0,0],[0,0,0],[0,0,0]
        for i in range(len(peff)):
            L[i], O[i], S[i] = peff[i].perform_sim(dbpath=this_dbpath, dt=dt, stoppers=stoppers)
            # print(S[i])
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

        print(T, dt)

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
    peffs = [get_peffs(T=T), get_peffs(T=T), get_peffs(T=T), get_peffs(T=T)]#, get_peffs(T=T)]
    for peff in peffs[1]:
        peff.params['mol_CO2'] = peff.params['mol_CO2']*0.99
    for peff in peffs[2]:
        peff.params['mol_CO2'] = peff.params['mol_CO2']*0.9
    # for peff in peffs[3]:
    #     peff.params['mol_CO2'] = peff.params['mol_CO2']*0.8
    for peff in peffs[3]:
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



def growthcurves(Ts=[280,300,320], dts=[3000,2000,3000-(40000/25)]):
    """plot growth curves (manuscript Figure 3)"""
    fig, axs = plt.subplots(nrows=len(Ts), ncols=4, figsize=(9,8), sharey=True)


    for i in range(len(Ts)):
        al=0.75

        ############################## CO2 CONCENTRATION
        peffs = CO2peffs(T=int(Ts[i]))
        inbetweens(peffs, axs[i][0], T=int(Ts[i]), dt=int(dts[i]), rmv=rmv, alph=al)

        ############################## H2 CONCENTRATION
        peffs = H2peffs(T=int(Ts[i]))
        inbetweens(peffs, axs[i][1], T=int(Ts[i]), dt=int(dts[i]), rmv=rmv, alph=al)

        ############################## P CONCENTRATION
        peffs = Ppeffs(T=int(Ts[i]))
        inbetweens(peffs, axs[i][2], T=int(Ts[i]), dt=int(dts[i]), rmv=rmv, alph=al)

        ############################## P UPTAKE
        peffs = Puptakepeffs(T=int(Ts[i]))
        inbetweens(peffs, axs[i][3], rmv=rmv, T=int(Ts[i]), dt=int(dts[i]), alph=al)

        # axs[0][i].set_title('T: '+str(int(Ts[i]))+' K')
    axs[0][0].set_title('[$\mathregular{CO}_{2}$]')
    axs[0][1].set_title('[$\mathregular{H}_{2}$]')
    axs[0][2].set_title('[$\mathregular{P}}$]')
    axs[0][3].set_title('$k_{P}$')

    axs[2][0].legend([n+' Optimal [$\mathregular{CO}_{2}$]' for n in ['', '-1%','-10%','-40%']], loc='upper center', bbox_to_anchor=(0.5, -0.2))
    axs[2][1].legend([n+' Optimal [$\mathregular{H}_{2}$]' for n in ['', '-1%','-5%','-20%']], loc='upper center', bbox_to_anchor=(0.5, -0.2))
    axs[2][2].legend(['[$\mathregular{P}}$] = '+n for n in ['$10^{-1}$ M (plentiful)', '$10^{-8}$ M', '$10^{-10}$ M']], loc='upper center', bbox_to_anchor=(0.5, -0.2))
    axs[2][3].legend(['$k_{P}$ = '+n+' $\mathregular{s}^{-1}$' for n in ['$10^{-10}$', '$10^{-18}$', '$10^{-20}$']], loc='upper center', bbox_to_anchor=(0.5, -0.2))

    axs[0][-1].annotate('Temperature: 280 K', xy=(1.1, 0.5), xytext=(1.2, 0.5), xycoords='axes fraction',
            fontsize=12, ha='left', va='center', rotation='vertical',
            arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=1.5))
    axs[1][-1].annotate('Temperature: 300 K', xy=(1.1, 0.5), xytext=(1.2, 0.5), xycoords='axes fraction',
            fontsize=12, ha='left', va='center', rotation='vertical',
            arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=1.5))
    axs[2][-1].annotate('Temperature: 330 K', xy=(1.1, 0.5), xytext=(1.2, 0.5), xycoords='axes fraction',
            fontsize=12, ha='left', va='center', rotation='vertical',
            arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1.0', lw=1.5))

    for ax in axs:
        ax[0].set_ylabel('Biomass [cells $\mathregular{L}^{-1}$]')
        for a in ax:

            a.set_yscale('log')
            a.set_xscale('log')
            a.set_xlim(1e4, 5e7)
            a.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.8)
            a.set_xlabel('Time [s]')

    plt.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('growthcurves.pdf')
    plt.show()

growthcurves()


def plot_Trange(energypeffs, nutrientpeffs, axs, rmv=False, stoppers={}, Tvals=range(280,331), hatch=None):
    """Function for plotting Figure 4.

    Plots one list of energy limitations and one of nutrient limitations at a
    time, passed in the form of efficiencies objects, energypeffs or nutrientpeffs.
    Pass rmv as True to delete results and perform simulations again.

    In a previous version the hatch option was used to hatch some of the plots,
    now it is used to plot the seperate axes in multiple calls to this funcion.
    """
    collist=['#2c7bb6','#abd9e9', '#fdae61', '#d7191c']
    energydf = pd.DataFrame(energypeffs)
    nutrientdf = pd.DataFrame(nutrientpeffs)
    print(energydf)
    for ic in range(len(energydf.columns)):
        # choose a good timestep. Helps with the fill_between.
        dt =3000-((num)*(1000/25))
        if T==300:
            dt=2000
        if T==325:
            dt=1000

        GR = [[],[],[]]
        BM = [[],[],[]]
        BS = [[],[],[]]

        for peff in energydf[ic]:

            L,O,S,t,pop= [0,0,0], [0,0,0],[0,0,0],[0,0,0],[0,0,0]
            for i in range(len(peff)):
                L[i], O[i], S[i] = peff[i].perform_sim(dbpath=this_dbpath, dt=dt, stoppers=stoppers)
                if rmv:
                    es.ecosystem_dbhelper.db_helper.removesimdata(S[i], dbpath=this_dbpath)
                    L[i], O[i], S[i] = peff[i].perform_sim(dbpath=this_dbpath, dt=dt, stoppers=stoppers)

                gr = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S[i], 'PeakGR', OrgID='', dbpath=this_dbpath)
                gr = ast.literal_eval(gr)[0]
                GR[i].append(gr*3600)
                bm = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S[i], 'FinBM_cells_tot', OrgID='', dbpath=this_dbpath)
                BM[i].append(bm)
                bs = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S[i], 'Composition', OrgID='', dbpath=this_dbpath)
                bs = ast.literal_eval(tuple(bs)[-1][0])['CH4(g)']
                BS[i].append(bs)


        if hatch!=None:
            axs[0][1].plot(Tvals,GR[1], c=collist[ic], alpha=0.75)#, ls='--')
            axs[1][1].plot(Tvals,BM[1], c=collist[ic], alpha=0.75)#, ls='--')
            axs[2][1].plot(Tvals,BS[1], c=collist[ic], alpha=0.75)#, ls='--')

            axs[0][1].fill_between(Tvals, GR[0], GR[2], alpha=0.75, edgecolor=collist[ic], facecolor=collist[ic])
            axs[1][1].fill_between(Tvals, BM[0], BM[2], alpha=0.75, edgecolor=collist[ic], facecolor=collist[ic])
            axs[2][1].fill_between(Tvals, BS[0], BS[2], alpha=0.75, edgecolor=collist[ic], facecolor=collist[ic])

        else:
            axs[0][0].plot(Tvals,GR[1], c=collist[ic], alpha=0.75)
            axs[1][0].plot(Tvals,BM[1], c=collist[ic], alpha=0.75)
            axs[2][0].plot(Tvals,BS[1], c=collist[ic], alpha=0.75)

            axs[0][0].fill_between(Tvals, GR[0], GR[2], alpha=0.75, edgecolor=collist[ic], facecolor=collist[ic])
            axs[1][0].fill_between(Tvals, BM[0], BM[2], alpha=0.75, edgecolor=collist[ic], facecolor=collist[ic])
            axs[2][0].fill_between(Tvals, BS[0], BS[2], alpha=0.75, edgecolor=collist[ic], facecolor=collist[ic])

    for ic in range(len(nutrientdf.columns)):

        GR = [[],[],[]]
        BM = [[],[],[]]
        BS = [[],[],[]]

        for peff in nutrientdf[ic]:

            L,O,S,t,pop= [0,0,0], [0,0,0],[0,0,0],[0,0,0],[0,0,0]
            for i in range(len(peff)):
                L[i], O[i], S[i] = peff[i].perform_sim(dbpath=this_dbpath, dt=dt, stoppers=stoppers)
                if rmv:
                    es.ecosystem_dbhelper.db_helper.removesimdata(S[i], dbpath=this_dbpath)
                    L[i], O[i], S[i] = peff[i].perform_sim(dbpath=this_dbpath, dt=dt, stoppers=stoppers)

                gr = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S[i], 'PeakGR', OrgID='', dbpath=this_dbpath)
                gr = ast.literal_eval(gr)[0]
                GR[i].append(gr*3600)
                bm = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S[i], 'FinBM_cells_tot', OrgID='', dbpath=this_dbpath)
                BM[i].append(bm)
                bs = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S[i], 'Composition', OrgID='', dbpath=this_dbpath)
                bs = ast.literal_eval(tuple(bs)[-1][0])['CH4(g)']
                BS[i].append(bs)


        if hatch!=None:
            axs[0][3].plot(Tvals,GR[1], c=collist[ic], alpha=0.6)#, ls='--')
            axs[1][3].plot(Tvals,BM[1], c=collist[ic], alpha=0.6)#, ls='--')
            axs[2][3].plot(Tvals,BS[1], c=collist[ic], alpha=0.6)#, ls='--')

            axs[0][3].fill_between(Tvals, GR[0], GR[2], alpha=0.6, edgecolor=collist[ic], facecolor=collist[ic])#, facecolor='none'), hatch=hatch)
            axs[1][3].fill_between(Tvals, BM[0], BM[2], alpha=0.6, edgecolor=collist[ic], facecolor=collist[ic])#, facecolor='none'), hatch=hatch)
            axs[2][3].fill_between(Tvals, BS[0], BS[2], alpha=0.6, edgecolor=collist[ic], facecolor=collist[ic])#, facecolor='none'), hatch=hatch)

        else:
            axs[0][2].plot(Tvals,GR[1], c=collist[ic], alpha=0.6)
            axs[1][2].plot(Tvals,BM[1], c=collist[ic], alpha=0.6)
            axs[2][2].plot(Tvals,BS[1], c=collist[ic], alpha=0.6)

            axs[0][2].fill_between(Tvals, GR[0], GR[2], alpha=0.6, edgecolor=collist[ic], facecolor=collist[ic])
            axs[1][2].fill_between(Tvals, BM[0], BM[2], alpha=0.6, edgecolor=collist[ic], facecolor=collist[ic])
            axs[2][2].fill_between(Tvals, BS[0], BS[2], alpha=0.6, edgecolor=collist[ic], facecolor=collist[ic])



def grbmbs(Tvals=range(280,331)):
    """Function for the final plot which became Figure 4 in the manuscript. """

    fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(9,8), sharey='row')

    # get typical optimal methanogen parameters
    _co2peffs=[]
    _h2peffs=[]
    _Ppeffs=[]
    _Puptakepeffs=[]
    for T in Tvals:
        _co2peffs.append(CO2peffs(T=T))
        _h2peffs.append(H2peffs(T=T))
        _Ppeffs.append(Ppeffs(T=T))
        _Puptakepeffs.append(Puptakepeffs(T=T))

    # plot the growth rates, biomass and methane
    plot_Trange(_co2peffs, _Ppeffs, axs, Tvals=Tvals)#, dt=3000)
    plot_Trange(_h2peffs, _Puptakepeffs, axs, Tvals=Tvals, hatch='x')#, dt=3000)

    # tidy up the plot
    axs[0][0].set_ylabel('Peak Growth Rate [$\mathregular{hr}^{-1}$]')
    axs[1][0].set_ylabel('Max Biomass [cells $\mathregular{L}^{-1}$]')
    axs[2][0].set_ylabel('$\mathregular{CH}_4$ Produced [mol $\mathregular{L}^{-1}$]')

    axs[0][0].set_title('[$\mathregular{CO}_2$]')
    axs[0][1].set_title('[$\mathregular{H}_2$]')
    axs[0][2].set_title('[P]')
    axs[0][3].set_title('$k_P$')

    for i, ax in enumerate(axs):
        for j, a in enumerate(ax):
            a.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.8)
            if i == 2:
                a.set_xlabel('Temperature [K]')
            else:
                a.tick_params(labelcolor='none', axis='x')
            if j!=0:
                a.tick_params(labelcolor='none', axis='y')
            a.set_yscale('log')

    axs[2][0].legend([n+' Optimal [$\mathregular{CO}_{2}$]' for n in ['', '-1%','-10%','-40%']], loc='upper center', bbox_to_anchor=(0.5, -0.2))
    axs[2][1].legend([n+' Optimal [$\mathregular{H}_{2}$]' for n in ['', '-1%','-5%','-20%']], loc='upper center', bbox_to_anchor=(0.5, -0.2))
    axs[2][2].legend(['[$\mathregular{P}}$] = '+n for n in ['$10^{-1}$ M (plentiful)', '$10^{-8}$ M', '$10^{-10}$ M']], loc='upper center', bbox_to_anchor=(0.5, -0.2))
    axs[2][3].legend(['$k_{P}$ = '+n+' $\mathregular{s}^{-1}$' for n in ['$10^{-10}$', '$10^{-18}$', '$10^{-20}$']], loc='upper center', bbox_to_anchor=(0.5, -0.2))

    plt.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('grbmbs.pdf')
    plt.show()

# grbmbs()






def k_RTPs(Trange=range(273,373)):
    """ export typical optimal methanogen parameters to a csv file"""
    Klst, MPlst, CO2lst, H2lst, CH4lst, Plst, vollst, metlst =[],[],[],[],[],[],[], []
    for T in Trange:
        pe = get_peffs(T=T)[1]
        Klst.append(pe.k_RTP)
        MPlst.append(pe.params['MP'])
        CO2lst.append(pe.params['mol_CO2'])
        H2lst.append(pe.params['mol_H2'])
        CH4lst.append(pe.params['mol_CH4'])
        Plst.append(pe.params['Pressure'])
        vollst.append((math.pi*(pe.params['radius'][0]**2)*pe.params['radius'][1]))

        R = es.reactor('r1met', workoutID=False, pH=pe.params['pH'], dbpath=this_dbpath)
        R.change_T(T)
        R.change_P(pe.params['Pressure'])

        pe.initial_conditions(R) # sets up composition, reaction

        volume, dry_mass, mass = pe.getvol()

        H = es.horde('h1met', R, pe.setup_methanogenesis(R), 500, Tdef=pe.params['Tdef'], mass=mass, dry_mass=dry_mass, volume=volume, workoutID=False, n_ATP=pe.params['n_ATP'], dbpath=this_dbpath)

        H.update_metabolic_rate()
        metlst.append(H.metabolic_rate)




    df = pd.DataFrame({'Temperature':Trange, 'k_RTP':Klst, 'Maintenance':MPlst, 'CO2':CO2lst, 'H2':H2lst, 'CH4':CH4lst, 'Pressure': Plst, 'Volume':vollst, 'MetabolicRate':metlst})
    df.to_csv('params2.csv')


# k_RTPs()
