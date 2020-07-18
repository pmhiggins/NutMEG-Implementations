import methanogen_extractor as extractor
import sys,os, statistics
sys.path.append(os.path.dirname(__file__)+'/../../NutMEG/')
import NutMEG as es
import NutMEG.util.NutMEGparams as nmp

import matplotlib.pyplot as plt
import matplotlib as mpl
# change matplotlib rcparams updating plot fonts etc.
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cmr10'
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset'] = 'cm'

def tovol(vol, lsts):
    """Convert lists passed into unit volume"""
    for l in lsts:
        i=0
        for n, v in zip(l,vol):
            l[i] = n/v
            i +=1

def plot_powers(ax, key, ilabel, files=('data/EmpiricalMethanogens/CH483.csv',), unitBM=False, theory=True, log=True, lowhigh=[], dbpath=nmp.std_dbpath):
    """Plot the maintenance power of the individual methanogens in the database """

    Edict = extractor.allmethanogens_fromcsv(filename=files[0], extra=True)
    Eyerr=[0]*len(Edict['Temp'])
    if len(lowhigh)==2:
        # errorbars are to be added
        lowErange = extractor.allmethanogens_fromcsv(filename=lowhigh[0])
        highErange = extractor.allmethanogens_fromcsv(filename=lowhigh[1])

        for k in range(len(Edict[key])):
            lowErange[key][k] = Edict[key][k] - lowErange[key][k]
            highErange[key][k] = highErange[key][k] - Edict[key][k]

        Eyerr=[lowErange[key], highErange[key]]

    if unitBM:
        tovol(Edict['vol'], [Edict[key]])

    if theory:
        Tr, L10,L2, Ti = range(273,380), [],[],[]
        meanvol = statistics.mean(Edict['vol'])
        TE = es.applications.theory_estimates.fromSim(Edict['SimID'][0], dbpath=dbpath)
        TE.org.base_volume = meanvol
        est={}
        if unitBM:
            for T in Tr:
                TE.loc.change_T(T)
                est = TE.temperature_defenses(T, per_cell=False)
                L10.append(est['Lever10pc'])
                L2.append(est['Lever2pc'])
                Ti.append(est['Tijhuis'])
        else:
            for T in Tr:
                TE.loc.change_T(T)
                est = TE.temperature_defenses(T, per_cell=True)
                L10.append(est['Lever10pc'])
                L2.append(est['Lever2pc'])
                Ti.append(est['Tijhuis'])


    if unitBM:
        ax.errorbar(Edict['Temp'],Edict[key], yerr=Eyerr, fmt='x', c='b', label=ilabel+' per L', ms=10, mew=2)
    else:
        ax.errorbar(Edict['Temp'],Edict[key], yerr=Eyerr, fmt='x', c='b', label=ilabel+ ' per cell', ms=10, mew=2)

    if theory and unitBM:
        ax.plot(Tr, Ti, c='k', label='Tijhuis maintenance')
        ax.plot(Tr, L10, c='g', label='Lever 10pc maintenance')
        ax.plot(Tr, L2, c='c', label='Lever 2pc maintenance')

    if theory and not unitBM:
        ax.plot(Tr, Ti, c='k', label='Tijhuis maintenance')
        ax.plot(Tr, L10, c='g', label='Lever 10pc maintenance')
        ax.plot(Tr, L2, c='c', label='Lever 2pc maintenance')

    if unitBM:
        ax.set_ylabel(ilabel+' Power per L biomass', fontsize=11)
        ax.set_ylim(1e-2,1e8)
    else:
        ax.set_ylabel(ilabel+' W/cell', fontsize=11)
        ax.set_ylim(1e-19,1e-8)

    ax.set_xlabel('Temperature [K]', fontsize=11)
    ax.set_xlim(273,373)
    if log:
        ax.set_yscale('log')


def MaintenanceRange_nATPs(Trange=[275,295,315,335,355,375], Perform=True, fraction=False, mCH4=1e-9, PS=True, Tlst=False, dbpath=nmp.std_dbpath):
    """Get the Predicted maintenance powers for this study, Tijhuis and Lever for the range of temperatures Trange."""
    PT, PT05, PT15 =[],[],[]
    Ti, Ti05, Ti15 =[],[],[]
    Lr10, Lr10_05, Lr10_15 = [],[],[]
    Lr2, Lr2_05, Lr2_15 = [],[],[]
    PS1, PS05, PS15 = [],[],[]

    for t in Trange:
        print(t)
        if Perform:
            mf,pt,ps,pg, s = extractor.iterateESynths([1.0], 'averageMethanogen', paramchange={'Tdef':'None', 'Temp':t, 'mol_CH4':mCH4}, save=('data/TOM_PT/'+str(mCH4)+'_'+str(t)), dbpath=dbpath)
            mf05,pt05,ps05,pg05, s05 = extractor.iterateESynths([1.0], 'averageMethanogen', paramchange={'Tdef':'None', 'Temp':t, 'mol_CH4':mCH4, 'n_ATP':0.5}, save=('data/TOM_PT/05_'+str(mCH4)+'_'+str(t)), dbpath=dbpath)
            mf15,pt15,ps15,pg15, s15 = extractor.iterateESynths([1.0], 'averageMethanogen', paramchange={'Tdef':'None', 'Temp':t, 'mol_CH4':mCH4, 'n_ATP':1.5}, save=('data/TOM_PT/15_'+str(mCH4)+'_'+str(t)), dbpath=dbpath)

        mf,pt,ps,pg,s = extractor.extract_Esynths_csv('data/TOM_PT/'+str(mCH4)+'_'+str(t))
        mf05,pt05,ps05,pg05,s05 = extractor.extract_Esynths_csv('data/TOM_PT/05_'+str(mCH4)+'_'+str(t))
        mf15,pt15,ps15,pg15,s15 = extractor.extract_Esynths_csv('data/TOM_PT/15_'+str(mCH4)+'_'+str(t))

        if fraction:
            PT.append(pt[0]/ps[0])
            PT05.append(pt05[0]/ps05[0])
            PT15.append(pt15[0]/ps15[0])
        else:
            PT.append(pt[0])
            PT05.append(pt05[0])
            PT15.append(pt15[0])
        if fraction == False and PS==True:
            PS1.append(ps[0])
            PS05.append(ps05[0])
            PS15.append(ps15[0])

        if Tlst:
            try:
                TE = es.applications.theory_estimates.fromSim(s[0], dbpath=dbpath)
            except:
                TE = es.applications.theory_estimates.fromSim(s[0], dbpath=nmp.std_dbpath)
            td = TE.temperature_defenses(t)
            ti = td['Tijhuis']
            lr = td['Lever10pc']
            lr2 = td['Lever2pc']

            # TE05 = es.applications.theory_estimates.fromSim(s05[0], dbpath=dbpath)
            # ti05 = TE05.temperature_defenses(t)['Tijhuis']
            #
            # TE15 = es.applications.theory_estimates.fromSim(s15[0], dbpath=dbpath)
            # ti15 = TE15.temperature_defenses(t)['Tijhuis']

            if fraction:
                Ti.append(ti/ps[0])
                Ti05.append(ti/ps05[0])
                Ti15.append(ti/ps15[0])

                Lr10.append(lr/ps[0])
                Lr10_05.append(lr/ps05[0])
                Lr10_15.append(lr/ps15[0])

                Lr2.append(lr2/ps[0])
                Lr2_05.append(lr2/ps05[0])
                Lr2_15.append(lr2/ps15[0])

            else:
                Ti.append(ti)
                Ti05.append(ti)
                Ti15.append(ti)

                Lr10.append(lr)
                Lr10_05.append(lr)
                Lr10_15.append(lr)

                Lr2.append(lr2)
                Lr2_05.append(lr2)
                Lr2_15.append(lr2)

    if Tlst:
        if PS:
            return [PT,PT05,PT15, PS1, PS05, PS15], [Ti, Ti05, Ti15], [Lr10, Lr10_05, Lr10_15], [Lr2, Lr2_05, Lr2_15]
        else:
            return [PT,PT05,PT15], [Ti, Ti05, Ti15], [Lr10, Lr10_05, Lr10_15], [Lr2, Lr2_05, Lr2_15]

    else:
        return [PT, PT05, PT15]








# es.db_helper.create_major_tables(replace=False, dbpath='allMtestc')

# extractor.allmethanogens_tocsv(filename='data/EmpiricalMethanogens/CH483.csv', nATP=1.0, ESfrac=1.0, mol_CH4=3e-8, dbpath='allMtestc')
# extractor.allmethanogens_tocsv(filename='data/EmpiricalMethanogens/05CH483.csv', nATP=0.5, ESfrac=1.0, mol_CH4=3e-8, dbpath='allMtestc')
# extractor.allmethanogens_tocsv(filename='data/EmpiricalMethanogens/15CH483.csv', nATP=1.5, ESfrac=1.0, mol_CH4=3e-8, dbpath='allMtestc')


fig, ax = plt.subplots(figsize=(6,6))
# plot_powers(axs[0], 'MF', 'Fraction', files=('data/EmpiricalMethanogens/CH483.csv',), unitBM=False, theory=False, dbpath='allMtest')#, lowhigh=['data/EmpiricalMethanogens/05CH483.csv', 'data/EmpiricalMethanogens/15CH483.csv'])
plot_powers(ax, 'PThrottle', 'Maintenance Power', files=('data/EmpiricalMethanogens/CH483.csv',), unitBM=False, theory=False, dbpath='allMtestc', lowhigh=['data/EmpiricalMethanogens/05CH483.csv', 'data/EmpiricalMethanogens/15CH483.csv'])


T = range(273,375)#[275,285,295,305,315,325,335,345,355,365,375]
# PT83, Ti, L10, L2 = TijhuisRange_nATPs(Trange=T, mCH4=3e-8, Tlst=True, fraction=False, Perform=False, dbpath=nmp.std_dbpath)
# PT83 = TijhuisRange_nATPs(Trange=[273,274], mCH4=3e-8,fraction=False, Perform=True, dbpath='allMtestc')
PT, Ti, L10, L2 = MaintenanceRange_nATPs(Trange=T, mCH4=3e-8, Tlst=True, fraction=False, Perform=False, dbpath='allMtestc')



ax.plot(T,PT[0], c='tab:orange', label="`Typical' Methanogen")
ax.fill_between(T, PT[1], PT[2], color='tab:orange', alpha=0.6)
# axs[i].plot(T,PT[3], c='r')
# axs[i].fill_between(T, PT[4], PT[5], color='r', alpha=0.6)
ax.plot(T,Ti[0], c='k', label='Empirical Estimate: Tijhuis et al. (1993)')
ax.fill_between(T, Ti[1], Ti[2], color='k', alpha=0.6)
ax.plot(T,L10[0], c='g', label='Theory Minimum: Lever et al. (2015) - 10% racemization')
ax.fill_between(T, L10[1], L10[2], color='g', alpha=0.6)
ax.plot(T,L2[0], c='c', label='Theory Minimum: Lever et al. (2015) - 2% racemization')
ax.fill_between(T, L2[1], L2[2], color='c', alpha=0.6)

ax.set_yscale('log')

plt.legend(loc='upper center', fontsize=11, bbox_to_anchor=(0.5,-.1))
plt.tight_layout()
# plt.subplots_adjust(hspace=0.0)
plt.savefig('MP.pdf')
plt.show()




""" THIS WAS USED TO CHECK OUT THE SUPPLY, THROTTLE AND FRACTIONS. MAY BE USEFUL IN THE FUTURE"""

# fig, axs = plt.subplots(nrows=3)
#
# T = [275,285,295,305,315,325,335,345,355,365,375]
# Edict = extractor.allmethanogens_fromcsv(filename='/data/EmpiricalMethanogens/CH483.csv', extra=True)
# PT83f, Tiaf = TijhuisRange_nATPs(Trange=T, mCH4=3e-8, Tlst=True, fraction=True, Perform=False, dbpath=nmp.std_dbpath)
# PT83, Tia = TijhuisRange_nATPs(Trange=T, mCH4=3e-8, Tlst=True, fraction=False, Perform=False, dbpath=nmp.std_dbpath)
#
# methTemps=Edict['Temp']
# methSupp=Edict['PSupply']
# methThrot=Edict['PThrottle']
# methFrac=[]
#
# for s, t in zip(methSupp, methThrot):
#     methFrac.append(t/s)
#
#
# avgSupp=[]
# avgThrot=[]
# avgFrac=[]
#
# for pt, ptf in zip(PT83[0], PT83f[0]):
#     try:
#         avgSupp.append(pt/ptf)
#     except:
#         avgSupp.append(1e-12)
#     avgThrot.append(pt)
#     avgFrac.append(ptf)
#
#
# axs[0].scatter(methTemps,methSupp)
# axs[0].plot(T,avgSupp)
# axs[0].set_ylabel('Supply')
# axs[0].set_yscale('log')
#
# axs[1].scatter(methTemps,methThrot)
# axs[1].plot(T,avgThrot)
# axs[1].set_ylabel('Throttle')
# axs[1].set_yscale('log')
#
# axs[2].scatter(methTemps,methFrac)
# axs[2].plot(T,avgFrac)
# axs[2].set_ylabel('Fraction')
#
# plt.show()
