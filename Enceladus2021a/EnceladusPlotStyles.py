import sys, os, ast, math
sys.path.append(os.path.dirname(__file__)+'../../NutMEG')
sys.path.append(os.path.dirname(__file__)+'../../NutMEG-Implementations/TOM')

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as clr
import numpy as np
import pandas as pd
from scipy import interpolate

import NutMEG as es
from NutMEG.reactor.saved_systems.Enceladus import Enceladus

import theory_emp_match as tem
import EnceladusGrids
from EnergyCalculations import getE_TOM
# import EnceladusSampler_2 as ES2


mpl.rcParams['contour.negative_linestyle'] = 'solid'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cmr10'
mpl.rcParams['font.serif'] = 'Times.tcc'
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rc('axes', unicode_minus=False)
mpl.rcParams['hatch.linewidth'] = 2.0
# print(mpl.rcParams.keys()) # to see what can be changed
mpl.rcParams['font.size'] = 14
# mpl.rcParams['font.weight'] = 'bold'



############## ADDING LINES ################

def add_pH_lines(ax, color='dimgray',
  pHnames = ['7.0','7.5','8.0','8.5','9.0','9.5','10.0','10.5','11.0','11.5','12.0'],
  unc=False):

    Tfloats = np.linspace(273.15, 473.15, num=21)
    # pHlist = np.ndarray((len(pHfloats), len(Tfloats)))

    for i, pH in enumerate(pHnames):
        df=pd.read_csv('E21data/Speciation/nominalCO2/pH'+pH+'.csv', sep=',')
        df2=pd.read_csv('E21data/Speciation/highsalt/pH'+pH+'.csv', sep=',')
        df3=pd.read_csv('E21data/Speciation/lowsalt/pH'+pH+'.csv', sep=',')
        ax.plot(df['pH'], Tfloats, c=color, linestyle='dashed', linewidth=3)
        if unc:
            ax.fill_betweenx(Tfloats, df2['pH'], df3['pH'], facecolor='y', alpha=0.3)
    # for the label
    ax.plot([0,0], [0,0], c=color, linestyle='dashed', linewidth=3, label='pH when warmed from 273.15 K')
    return ax

def add_pH_boxes(ax, CDAINMS=True, postberg=True):
    pHnames = ['7.0','7.5','8.0','8.5','9.0','9.5','10.0','10.5','11.0','11.5','12.0']
    Tfloats = np.linspace(273.15, 473.15, num=21)
    gleinlists=[]
    postlists=[]
    if postberg:
        for pH in ['8.5', '9.0']:
            df=pd.read_csv('E21data/Speciation/nominalCO2/pH'+pH+'.csv', sep=',')
            postlists.append(df['pH'])
        ax.fill_betweenx(Tfloats, postlists[0], postlists[1], color='none', edgecolor='tab:brown', hatch='X', linewidth=4, label='pH range Postberg 2009, Glein 2020')
    if CDAINMS:
        for pH in ['9.5', '11.0']:
            pd.read_csv('E21data/Speciation/nominalCO2/pH'+pH+'.csv', sep=',')
            gleinlists.append(df['pH'])
        ax.fill_betweenx(Tfloats, gleinlists[0], gleinlists[1], color='none', edgecolor='g', hatch='.',linewidth=4, label='pH range to match CDA to INMS (using Waite 2017 composition)')


def add_maintenance_lines(ax, colors=['tab:orange', 'k','y','c']):
    T1 = range(273, 375)
    PT = tem.MaintenanceRange_nATPs(Trange=T1, mCH4=3e-8, Tlst=False, fraction=False, Perform=False, dbpath='../../NutMEG-Implementations/TOM/allMtestc')

    ax.fill_between(T1, np.log10(PT[1]), np.log10(PT[2]), color=colors[0], alpha=0.6)

    T2 = range(273, 400)
    # make a TOM and an Enceladus
    _Enc = Enceladus('EncT')
    TOM = getE_TOM(_Enc)

    TE = es.applications.theory_estimates(TOM, _Enc)
    Ti, L10, L2 = [],[],[]
    for t in T2:
        td = TE.temperature_defenses(t)
        Ti.append(td['Tijhuis'])
        L10.append(td['Lever10pc'])
        L2.append(td['Lever2pc'])

    ax.fill_between(T2, np.log10(Ti), np.log10(Ti), color=colors[1], alpha=0.8, lw=3.)

    ax.fill_between(T2, np.log10(L10), np.log10(L2), color=colors[2], alpha=0.4)

    ax.fill_between([T2[0],T2[-1]], [-18, -18], [-21, -21], color=colors[3], alpha=0.4)
    return ax

def add_maintenance_labels(ax, colors=['tab:orange', 'k','g','c']):
    ax.fill_between([0,0], [1,1],[1,1], color=colors[0], alpha=0.6, label='Maintenance Power in optimum conditions (Higgins & Cockell 2020)')
    ax.fill_between([0,0], [1,1],[1,1], color=colors[1], alpha=0.6, label='Maintenance Power for anaerobes (Tijhuis et al. 1998)')
    ax.fill_between([0,0], [1,1],[1,1], color=colors[2], alpha=0.6, label='Minimal Maintenance Power (Lever et al. 2015)')
    ax.fill_between([0,0], [1,1],[1,1], color=colors[3], alpha=0.6, label='Min. subsurface power supplies (Bradley et al. 2020)')
    return ax


def add_pH_maintenance_lines(ax, T, pHrange=[7,12], colors=['tab:orange', 'k','g','c']):
    PT = tem.MaintenanceRange_nATPs(Trange=[T], mCH4=3e-8, Tlst=False, fraction=False, Perform=False, dbpath='../../NutMEG-Implementations/TOM/allMtestc')

    # pH is only 5--9 as that's as far as the data goes
    ax.fill_between([5,9], np.log10(PT[1]), np.log10(PT[2]), color=colors[0], alpha=0.4)

    _Enc = Enceladus('EncT')
    TOM = getE_TOM(_Enc)

    TE = es.applications.theory_estimates(TOM, _Enc)
    Ti, L10, L2 = [],[],[]
    td = TE.temperature_defenses(T)
    Ti.append(td['Tijhuis'])
    L10.append(td['Lever10pc'])
    L2.append(td['Lever2pc'])

    ax.fill_between(pHrange, np.log10(Ti), np.log10(Ti), color=colors[1], alpha=0.8, lw=3.)

    ax.fill_between(pHrange, np.log10(L10), np.log10(L2), color=colors[2], alpha=0.4)

    ax.fill_between(pHrange, [-18, -18], [-21, -21], color=colors[3], alpha=0.4)
    return ax


def energycolormap():
    # top = plt.cm.get_cmap('winter_r', 128)
    cyan_to_b = {'red':[(0,0,0), (1,0,0)],
      'green':[(0,1,1),(1,0,0)],
      'blue':[(0,1,1), (1,1,1)],
      'alpha':[(0,1,1), (1,1,1)]}

    blue_to_a = {'red':[(0,0,0), (1,0,0)],
      'green':[(0,0,0),(1,0,0)],
      'blue':[(0,1,1), (1,1,1)],
      'alpha':[(0,1,1), (1,0,0)]}

    top = clr.LinearSegmentedColormap('ctb', cyan_to_b)
    middle = clr.LinearSegmentedColormap('bta', blue_to_a)
    bottom = plt.cm.get_cmap('Reds', 64)

    newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       middle(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
    newcmp = clr.ListedColormap(newcolors, name='energycol')
    return newcmp



###############   Methanogenesis Energy ##############

def MethanogenesisEnergyContourPlot(ax, CO2origin='pH',
  pHrange=np.linspace(7,14, num=11), Trange = np.linspace(273,473, num=11), maincont=0, quotienttype='salty_nominal'):

    plotkey = 'Gibbs_Methanogenesis'

    Meshes = EnceladusGrids.getMesh(
    Trange, pHrange, params=[plotkey], CO2origin=CO2origin, quotienttype=quotienttype)

    XX = Meshes['pH']
    YY = Meshes['T']
    ZZ = Meshes[plotkey]

    levels = [-150000, -75000, 0, 75000]

    # contf = ax.contourf(XX, YY, ZZ[maincont], levels=np.arange(-150000, 150000, 3000), cmap=energycolormap(), vmin=-150000, vmax=150000, extend='both')
    contf = ax.contourf(XX, YY, ZZ[maincont], levels=np.arange(-160000, 80000, 2500), cmap=energycolormap(), vmin=-160000, vmax=80000, extend='both')


    contourcolor='k'#'dimgray'

    cont = ax.contour(XX, YY, ZZ[0], levels=levels, colors=[contourcolor], vmin=-150000, vmax=80000, linewidths=2.5)
    ax.clabel(cont, inline=1, fontsize=16, fmt='%d') # add label

    ax.contour(XX, YY, ZZ[1], levels=levels, linestyles='dotted', colors=[contourcolor], vmin=-150000, vmax=80000, linewidths=2.5)
    ax.contour(XX, YY, ZZ[2], levels=levels, linestyles='dotted', colors=[contourcolor], vmin=-150000, vmax=80000, linewidths=2.5)

    ax.plot([0,0], [0,0], c= contourcolor, label='Uncertainty bounds on free energy contour', linestyle='dotted', linewidth=2.5)

    return ax, contf


def make_MGEContourPlot(CO2origin='pH', save='Energyplot.pdf', show=False, pHax=False, pHbars=False, quotienttype='salty_nominal'):

    fig, ax = plt.subplots(figsize=(8,8))
    ax, contf = MethanogenesisEnergyContourPlot(ax, CO2origin=CO2origin, pHrange=np.linspace(7,12, num=11), quotienttype=quotienttype)

    fig.colorbar(contf, label='Free Energy of methanogenesis [J/mol]', orientation='horizontal', pad=0.12)
    if pHax:
        ax = add_pH_lines(ax, pHnames = ['7.0','8.0','9.0','10.0','11.0','12.0'])
    if pHbars:
        add_pH_boxes(ax, CDAINMS=False)
    ax.set_xlim(7,12)
    ax.set_xlabel('Wider ocean pH (e.g. at 273.15 K)')
    ax.set_ylabel('Temperature [K]')
    ax.set_ylim(273.15,473)
    ax.legend(bbox_to_anchor=(0., 1.05, 1., .102), loc=3,
           ncol=1, mode="expand", borderaxespad=0.)
    plt.tight_layout()
    if save != None:
        plt.savefig(save)
    if show:
        # plt.tight_layout()
        plt.show()

def MGEcomparison_plot():

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8,8))
    ax[0], contf = MethanogenesisEnergyContourPlot(ax[0],
      pHrange=np.linspace(7,12, num=11), CO2origin='pH')
    ax[1], contf2 = MethanogenesisEnergyContourPlot(ax[1],
      pHrange=np.linspace(7,12, num=11), CO2origin='HTHeating')

    fig.colorbar(contf, label='Free Energy of methanogenesis [J/mol]', orientation='vertical', pad=0.08)
    plt.show()




###############   ATP Energy  ##############


def ATPEnergyContourPlot(ax, CO2origin='pH',
  pHrange=np.linspace(7,14, num=11), Trange = np.linspace(273,473, num=11)):

    plotkey = 'ATPGibbs'

    Meshes = EnceladusGrids.getMesh(
    Trange, pHrange, params=[plotkey])

    XX = Meshes['pH']
    YY = Meshes['T']
    ZZ = Meshes[plotkey]

    contf = ax.contourf(XX, YY, ZZ[0], levels=np.arange(55000,68000,1000), cmap='cool_r', vmin=55000, vmax=68000)

    cont = ax.contour(XX, YY, ZZ[0], levels=5, colors=['k'], vmin=55000, vmax=68000)
    ax.clabel(cont, inline=1, fontsize=12, fmt='%d') # add label

    return ax, contf

def make_ATPContourPlot(CO2origin='pH', save='ConservableEnergyplot.pdf', show=False):

    fig, ax = plt.subplots(figsize=(8,8))
    ax, contf = ATPEnergyContourPlot(ax, CO2origin=CO2origin)

    fig.colorbar(contf, label='Free Energy of ATP Synthesis [J/mol ATP]', orientation='horizontal', pad=0.08)

    # ax.legend(bbox_to_anchor=(0., 1.15, 1., .102), loc=3,
    #        ncol=1, mode="expand", borderaxespad=0.)
    if save != None:
        plt.savefig(save)
    if show:
        plt.show()




############# POWER SUPPLY ################


def PSContourPlot(ax, CO2origin='pH',
  pHrange=np.linspace(7,14, num=11), Trange = np.linspace(273,473, num=21),
  cmap='PuBu', maincont=0, mesh=False, k_corr=0.0, nATP=1.0):

    plotkey = 'PowerSupply'

    Meshes = EnceladusGrids.getMesh(
    Trange, pHrange, params=[plotkey], CO2origin=CO2origin, k_corr=k_corr, nATP=nATP)

    XX = Meshes['pH']
    YY = Meshes['T']
    ZZ = Meshes[plotkey]

    contourlevels = [-25,-15,-10, -5]

    cmap = plt.cm.get_cmap(cmap, 15)

    if mesh:
        contf = ax.pcolormesh(XX, YY, np.log10(ZZ[maincont]), vmin=-25, vmax=-5, shading='nearest', cmap=cmap, edgecolor='slategray', linewidth=1)

    else:
        contf = ax.contourf(XX, YY, np.log10(ZZ[maincont]), levels=np.arange(-25,-4,1), cmap=cmap, vmin=-25, vmax=-5, extend='both')

    return ax, contf


def make_PSContourPlot(CO2origin='pH', save='Powersupply.pdf', show=False, mesh=False, Tline=True):

    fig, ax = plt.subplots(figsize=(8,8))
    ax, contf = PSContourPlot(ax, CO2origin=CO2origin, pHrange=np.linspace(7,12, num=11), mesh=mesh, cmap='BuPu')

    fig.colorbar(contf, label='log10 (Power Supply [W/cell])', orientation='horizontal', pad=0.08, extend='both')

    ax.set_xlim(6.75,12.25)
    ax.set_xlabel('Wider ocean pH (e.g. at 273.15 K)')
    ax.set_ylabel('Temperature [K]')
    ax.set_ylim(268.15,478)
    if Tline:
        ax.axhline(400, c='tab:orange', lw=9)
        txt = ax.text(12.2, 402, 'MAX. TEMPERATURE LIMIT FOR LIFE',
          horizontalalignment='right',
          verticalalignment='bottom', color='tab:orange',
          family='serif', fontsize=11, fontweight='extra bold', variant='small-caps')

    plt.tight_layout()
    if save != None:
        plt.savefig(save)
    if show:
        plt.show()



def PScomparison_plot():

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8,8),)
    ax[0], contf = PSContourPlot(ax[0],
      pHrange=np.linspace(7,12, num=11), CO2origin='pH')
    ax[1], contf2 = PSContourPlot(ax[1],
      pHrange=np.linspace(7,12, num=11), CO2origin='HTHeating')

    plt.show()



def PSunc_plot(CO2origin='HTHeating20', save='Powersupply_comp.pdf', show=False, mesh=False, Tline=True, nATP=1.0):

        fig, axs = plt.subplots(figsize=(10,6), nrows=1, ncols=2, constrained_layout=True)

        axs[0], contf = PSContourPlot(axs[0], CO2origin=CO2origin, pHrange=np.linspace(7,12, num=11), mesh=mesh, cmap='BuPu', maincont=1, k_corr=-1., nATP=nATP)
        axs[0].set_title('Lower Bound')

        axs[1], contf = PSContourPlot(axs[1], CO2origin=CO2origin, pHrange=np.linspace(7,12, num=11), mesh=mesh, cmap='BuPu', maincont=2, k_corr=1., nATP=nATP)
        axs[1].set_title('Upper Bound')


        fig.colorbar(contf, ax=[axs[0], axs[1]], label='log10 (Power Supply [W/cell])', orientation='horizontal', pad=0.001, extend='both')
        for ax in axs:
            # ax.set_xlim(7,12)
            ax.set_xlim(6.75,12.25)

            ax.set_xlabel('Wider ocean pH (e.g. at 273.15 K)')
            ax.set_ylabel('Temperature [K]')
            # ax.set_ylim(273.15,473)
            ax.set_ylim(268.15,478)

            if Tline:
                ax.axhline(400, c='tab:orange', lw=9)
                txt = ax.text(12.2, 402, 'MAX. TEMPERATURE LIMIT FOR LIFE',
                  horizontalalignment='right',
                  verticalalignment='bottom', color='tab:orange',
                  family='serif', fontsize=11, fontweight='extra bold', variant='small-caps')
        if save != None:
            plt.savefig(save)
        if show:
            plt.show()


def PShabitabilityPlot(Trange = np.linspace(273,403, num=14), pHrange=np.linspace(7,12, num=11), nATP=1.0):

    mpl.rcParams['xtick.labelsize'] = 13
    mpl.rcParams['ytick.labelsize'] = 13
    mpl.rcParams['font.size'] = 13

    nomcols = plt.get_cmap('YlOrRd', 6)
    cmaplist = [nomcols(i) for i in range(nomcols.N)]
    # force the first color entry to be grey
    cmaplist[0] = (.95, .95, .95, 1.0)

    # create the new map
    nomcols = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, nomcols.N)

    # get a mesh for the 5 cases.
    fn_preamble=''
    if nATP == 1:
        fn_preamble='E21data/PSgrids/'
    else:
        fn_preamble='E21data/PSgrids/'+str(nATP)+'_'


    # for generating and saving

    # noms = EnceladusGrids.getMesh(Trange, pHrange, params=['PowerSupply'], nATP=nATP, quotienttype='salty_endmember')
    # np.save(fn_preamble+'noms.npy', noms['PowerSupply'][0])
    # np.save(fn_preamble+'pHgrid.npy', noms['pH'])
    # np.save(fn_preamble+'Tgrid.npy', noms['T'])
    #
    # nomhigh = EnceladusGrids.getMesh(Trange, pHrange, params=['PowerSupply'], nATP=nATP, k_corr=1., quotienttype='salty_high')['PowerSupply'][0]
    # highest = EnceladusGrids.getMesh(Trange, pHrange, params=['PowerSupply'], nATP=nATP, k_corr=1., quotienttype='salty_high')['PowerSupply'][2]
    #
    # nomlow = EnceladusGrids.getMesh(Trange, pHrange, params=['PowerSupply'], nATP=nATP, k_corr=-1., quotienttype='salty_low')['PowerSupply'][0]
    # lowest = EnceladusGrids.getMesh(Trange, pHrange, params=['PowerSupply'], nATP=nATP, k_corr=-1., quotienttype='salty_low')['PowerSupply'][1]
    #
    # np.save(fn_preamble+'nomhigh.npy', nomhigh)
    # np.save(fn_preamble+'highest.npy', highest)
    # np.save(fn_preamble+'nomlow.npy', nomlow)
    # np.save(fn_preamble+'lowest.npy', lowest)


    noms = np.load(fn_preamble+'noms.npy')
    pHgrid = np.load(fn_preamble+'pHgrid.npy')
    Tgrid = np.load(fn_preamble+'Tgrid.npy')
    nomhigh = np.load(fn_preamble+'nomhigh.npy')
    highest = np.load(fn_preamble+'highest.npy')
    nomlow = np.load(fn_preamble+'nomlow.npy')
    lowest = np.load(fn_preamble+'lowest.npy')


    # get the maintenance dictionaries
    HC, Ti, L, B = EnceladusGrids.maintenancemesh(Trange, pHrange, nATP=nATP)

    Bpass = np.ndarray((len(Trange), len(pHrange)))
    Lpass = np.ndarray((len(Trange), len(pHrange)))
    Tipass = np.ndarray((len(Trange), len(pHrange)))
    HCpass = np.ndarray((len(Trange), len(pHrange)))
    # print(HCpass)

    yindex = 0
    for T in Trange:
        xindex = 0
        for pH in pHrange:
            for mp, mppass in zip([HC,Ti,L,B], [HCpass, Tipass, Lpass, Bpass]):
                mppass[yindex][xindex] = 0
                if lowest[yindex][xindex] > mp[yindex][xindex]:
                    mppass[yindex][xindex] = 1
                elif nomlow[yindex][xindex] > mp[yindex][xindex]:
                    mppass[yindex][xindex] = 2
                elif noms[yindex][xindex] > mp[yindex][xindex]:
                    mppass[yindex][xindex] = 3
                elif nomhigh[yindex][xindex] > mp[yindex][xindex]:
                    mppass[yindex][xindex] = 4
                elif highest[yindex][xindex] > mp[yindex][xindex]:
                    mppass[yindex][xindex] = 5

            xindex += 1
        yindex += 1

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10,11))
    contf = ax[0][0].pcolormesh(pHgrid, Tgrid, HCpass, vmin=0, vmax=5, shading='nearest', cmap=nomcols, edgecolor='slategray', linewidth=1)
    contf = ax[0][1].pcolormesh(pHgrid, Tgrid, Tipass, vmin=0, vmax=5, shading='nearest', cmap=nomcols, edgecolor='slategray', linewidth=1)
    contf = ax[1][0].pcolormesh(pHgrid, Tgrid, Lpass, vmin=0, vmax=5, shading='nearest', cmap=nomcols, edgecolor='slategray', linewidth=1)
    contf = ax[1][1].pcolormesh(pHgrid, Tgrid, Bpass, vmin=0, vmax=5, shading='nearest', cmap=nomcols, edgecolor='slategray', linewidth=1)

    # fig.colorbar(contf)
    ax[0][0].plot([0],[0], c=cmaplist[0], label='Not enough power available in any configuration', linewidth=10)
    ax[0][0].plot([0],[0], c=cmaplist[1], label='Enough power available in the worst-case low-salt scenario', linewidth=10)
    ax[0][0].plot([0],[0], c=cmaplist[2], label='Enough power available in the nominal low-salt scenario', linewidth=10)
    ax[0][0].plot([0],[0], c=cmaplist[3], label='Enough power available in the nominal scenario', linewidth=10)
    ax[0][0].plot([0],[0], c=cmaplist[4], label='Enough power available in the nominal high-salt scenario', linewidth=10)
    ax[0][0].plot([0],[0], c=cmaplist[5], label='Enough power available in the best-case high-salt scenario', linewidth=10)

    ax[0][0].set_title('To exceed maximum optimal maintenance \n (Higgins+ 2020)', fontsize=13)
    ax[0][1].set_title('To exceed empirical maintenance for anaerobes \n (Tijhuis+ 1993)', fontsize=13)
    ax[1][0].set_title('To exceed minimal maintenance power \n (Lever+ 2015)', fontsize=13)
    ax[1][1].set_title('To exceed minimal Earth subsurface power supplies \n (Bradley+ 2020)', fontsize=12)


    for aa in ax:
        for a in aa:
            a.set_xlim(6.75, 12.25)
            a.set_ylim(268, 408)
            a.set_xlabel('Wider ocean pH')
            a.set_ylabel('Temperature [K]')

    plt.subplots_adjust(wspace=0.2, hspace=0.33, left=0.08, right=0.98, bottom=0.05, top=0.78)

    ax[0][0].legend(bbox_to_anchor=(0., 1.18, 2.15, .102), loc=3,
           ncol=1, mode="expand", borderaxespad=0.)

    if nATP == 1.0:
        plt.savefig('figs/habgrid.pdf')
    else:
        plt.savefig('figs/habgrid_'+str(nATP)+'.pdf')
