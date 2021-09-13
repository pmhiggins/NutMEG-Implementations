import sys, os, math
sys.path.append(os.path.dirname(__file__)+'../../NutMEG')
import ThesisSetup

import NutMEG as nm
from NutMEG.reactor.saved_systems.Enceladus import Enceladus
from NutMEG.culture.saved_organisms.TypicalOptimalMethanogen import TypicalOptimalMethanogen as TOM

import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib as mpl
import numpy as np

def SAget(v):
    return (4*math.pi*((v*3)**2))**(1/3)

def pH_flux_power(org, memb_pot, PermH, PermOH, v):
    fluxH, fluxOH, pow = [],[],[]
    org.memb_pot = memb_pot
    org.PermH = PermH
    org.PermOH = PermOH
    org.base_volume = v
    org.surfacearea = SAget(v)
    pHs = np.linspace(0,14, num=5000)
    for i, pH in enumerate(pHs):
        org.locale.update_pH(pH, _from='pH')
        pHad = nm.pHadaptations(org)
        org.maintenance.get_P_pH()
        pow.append(org.maintenance.net_dict['pH'])
        fluxH.append(abs(pHad._getfluxH()))
        fluxOH.append(abs(pHad._getfluxOH()))
    # return a dictionary of the pH values, fluxes, and power demand
    return {'pH':pHs,
      'fluxH':fluxH,
      'fluxOH':fluxOH,
      'pow':pow}



def fluxplot():
    # make a reactor object
    rtr = Enceladus('EncT', nominals=True, T=298)
    # make an organism, for this example, use the preset methanogen.
    org = TOM(rtr)
    org.maintenance.Tdef='None'
    org.maintenance.pHdef='FluxPerm'

    def_vals = pH_flux_power(org, -1e-3, 1e-10, 1e-10, org.base_volume)
    high_phi = pH_flux_power(org, -1e-1, 1e-10, 1e-10, org.base_volume)
    high_perms = pH_flux_power(org, -1e-3, 1e-5, 1e-5, org.base_volume)
    bigger = pH_flux_power(org, -1e-3, 1e-10, 1e-10, org.base_volume*100)

    fig, ax = plt.subplots(figsize=(8,4.5), ncols=2)

    cols = ['k', 'tab:blue', 'cyan', 'tab:purple']
    alphas = [1, 0.8, 0.4, 0.6]

    ax[0].plot(def_vals['pH'], def_vals['fluxH'], c=cols[0],
      label=r'NutMEG default values:     $\bar{P}_H = 10^{-10},\ \bar{P}_{OH}=10^{-10},\ \Delta\Psi = -10^{-3},\ v = 3.44\times10^{-18}$', alpha=alphas[0])
    ax[0].plot(def_vals['pH'], def_vals['fluxOH'], c=cols[0], ls='dashed', alpha=alphas[0])
    ax[1].plot(def_vals['pH'], def_vals['pow'], c=cols[0], alpha=alphas[0])

    ax[0].plot(high_phi['pH'], high_phi['fluxH'], c=cols[1],
      label=r'Increase potential:             $\bar{P}_H = 10^{-10},\ \bar{P}_{OH}=10^{-10},\ \Delta\Psi = -0.1,\ v = 3.44\times10^{-18}$', alpha=alphas[1])
    ax[0].plot(high_phi['pH'], high_phi['fluxOH'], c=cols[1], ls='dashed', alpha=alphas[1])
    ax[1].plot(high_phi['pH'], high_phi['pow'], c=cols[1], alpha=alphas[1])

    ax[0].plot(high_perms['pH'], high_perms['fluxH'], c=cols[2],
      label=r'Increase permeabilities:      $\bar{P}_H = 10^{-5},\ \bar{P}_{OH}=10^{-5},\ \Delta\Psi = -10^{-3},\ v = 3.44\times10^{-18}$', alpha=alphas[2])
    ax[0].plot(high_perms['pH'], high_perms['fluxOH'], c=cols[2], ls='dashed', alpha=alphas[2])
    ax[1].plot(high_perms['pH'], high_perms['pow'], c=cols[2], alpha=alphas[2])

    ax[0].plot(bigger['pH'], bigger['fluxH'], c=cols[3],
      label=r'Cell size:                         $\bar{P}_H = 10^{-10},\ \bar{P}_{OH}=10^{-10},\ \Delta\Psi = -10^{-3},\ v = 3.44\times10^{-16}$', alpha=alphas[3])
    ax[0].plot(bigger['pH'], bigger['fluxOH'], c=cols[3], ls='dashed', alpha=alphas[3])
    ax[1].plot(bigger['pH'], bigger['pow'], c=cols[3], ls='dashed', alpha=alphas[3])

    ax[0].set_title('Goldman-style fluxes for H$^{+}$ and OH$^{-}$')
    ax[1].set_title('Corresponding power demands')
    ax[0].set_ylabel('Flux of ions [mol (s cell$^{-1}$)]')
    ax[1].set_ylabel('Power demand [W cell$^{-1}$]')
    for a in ax:
        a.set_yscale('log')
        a.set_xlabel('external pH')

    fig.subplots_adjust(wspace=0.25, top=0.65, right=0.995, left=0.1)
    ax[0].legend(bbox_to_anchor=(0., 1.12, 2.25, 1.4), loc=3,
           ncol=1, mode="expand", borderaxespad=0.)
    # plt.show()
    plt.savefig('fluxplot.pdf')

fluxplot()


def complex_fluxes():
    rtr = Enceladus('EncT', nominals=True, T=298)
    # make an organism, for this example, use the preset methanogen.
    org = TOM(rtr)
    org.maintenance.Tdef='None'
    org.maintenance.pHdef='FluxPerm'
    phidicts, phidicts2, pHdicts, pHdicts2 = [], [], [], []
    phis= [1, 1e-1, -1e-3, -1e-1, -1]
    pHints= [5,6,7,8,9]

    nomcols = plt.get_cmap('coolwarm', len(phis)+1)
    phicmaplist = [nomcols(i) for i in range(nomcols.N)][1:]
    cmap = clr.LinearSegmentedColormap.from_list(
    'Custom cmap', phicmaplist, len(phis)+1)

    nomcols2 = plt.get_cmap('PRGn', len(pHints)+3)
    pHcmaplist = [nomcols2(i) for i in range(nomcols2.N)][1:-2]
    pHcmaplist = ['#984ea3', '#fb9a99', '#fdcdac', '#b2df8a', '#33a02c']
    cmap2 = clr.LinearSegmentedColormap.from_list(
    'Custom cmap2', pHcmaplist, len(pHints)+1)

    for i, phi in enumerate(phis):
        phidicts.append(pH_flux_power(org, phi, 1e-10, 1e-10, org.base_volume))
        phidicts2.append(pH_flux_power(org, phi, 1e-10, 1e-13, org.base_volume))

    for i, pHint in enumerate(pHints):
        org.pH_interior = pHint
        pHdicts.append(pH_flux_power(org, -1e-3, 1e-10, 1e-10, org.base_volume))
        pHdicts2.append(pH_flux_power(org, -1e-3, 1e-10, 1e-13, org.base_volume))

    fig, ax = plt.subplots(figsize=(7,8), ncols=2, nrows=3)

    for i in range(len(phidicts)):
        ax[0][0].plot(phidicts[i]['pH'], phidicts[i]['fluxH'], c=phicmaplist[i], linewidth=3)
        ax[0][0].plot(phidicts[i]['pH'], phidicts[i]['fluxOH'], ls='dashed',  c=phicmaplist[i], linewidth=3)
        ax[1][0].plot(phidicts[i]['pH'], phidicts[i]['pow'], c=phicmaplist[i], linewidth=3)
        ax[2][0].plot(phidicts2[i]['pH'], phidicts2[i]['pow'], c=phicmaplist[i], linewidth=3)
    for i in range(len(pHdicts)):
        ax[0][1].plot(pHdicts[i]['pH'], pHdicts[i]['fluxH'], c=pHcmaplist[i], linewidth=3)
        ax[0][1].plot(pHdicts[i]['pH'], pHdicts[i]['fluxOH'], ls='dashed',  c=pHcmaplist[i], linewidth=3)
        ax[1][1].plot(pHdicts[i]['pH'], pHdicts[i]['pow'], c=pHcmaplist[i], linewidth=3)
        ax[2][1].plot(pHdicts2[i]['pH'], pHdicts2[i]['pow'], c=pHcmaplist[i], linewidth=3)

    for ai in ax:
        for a in ai:
            a.set_yscale('log')
            a.set_xlim(0,14)
            a.get_xaxis().set_ticks([1,3,5,7,9,11,13])
            a.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.8)

    for a in ax[-1]:
        a.set_xlabel('external pH')
        a.set_ylim(1e-27,1e-12)
    for a in ax[1]:
        a.set_ylim(1e-27,1e-12)
        a.get_xaxis().set_ticklabels([])
    for a in ax[0]:
        a.set_ylim(1e-31, 1e-16)
        a.get_xaxis().set_ticklabels([])
    for a in [ax[0][1], ax[1][1], ax[2][1]]:
        a.get_yaxis().set_ticklabels([])

    ax[0][0].set_ylabel('Flux of ions [mol (s cell$^{-1}$)]')
    ax[1][0].set_ylabel('Power demand [W cell$^{-1}$]')
    ax[2][0].set_ylabel('Power demand [W cell$^{-1}$]')

    ax[0][0].text(0.42,0.9,'(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0][0].transAxes)
    ax[0][1].text(0.42,0.9, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[0][1].transAxes)
    ax[1][0].text(0.42,0.9, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[1][0].transAxes)
    ax[1][1].text(0.42,0.9, '(d)', horizontalalignment='center', verticalalignment='center', transform=ax[1][1].transAxes)
    ax[2][0].text(0.42,0.9, '(e)', horizontalalignment='center', verticalalignment='center', transform=ax[2][0].transAxes)
    ax[2][1].text(0.42,0.9, '(f)', horizontalalignment='center', verticalalignment='center', transform=ax[2][1].transAxes)


    cbaxes = fig.add_axes([0.125, 0.91, 0.425, 0.02])
    cbaxes2 = fig.add_axes([0.57, 0.91, 0.425, 0.02])

    # plt.colorbar(clr.BoundaryNorm(phis, nomcols.N), cax=cbaxes, label='phi', orientation='horizontal', pad=0.5, aspect=7)

    mpl.colorbar.ColorbarBase(cbaxes, cmap=cmap, norm=clr.BoundaryNorm([0.5,1.5,2.5,3.5,4.5], nomcols.N-1),
    spacing='proportional', ticks=[0.5,1.5,2.5,3.5,4.5], boundaries=[0,1,2,3,4,5], format='%1i', orientation='horizontal')
    cbaxes.set_xticklabels(['1','10$^{-1}$', '$-10^{-3}$', '$-10^{-1}$', '-1'])
    cbaxes.set_xlabel('Membrane potential [V]')
    cbaxes.xaxis.set_label_position('top')
    cbaxes.xaxis.tick_top()

    mpl.colorbar.ColorbarBase(cbaxes2, cmap=cmap2, norm=clr.BoundaryNorm(pHints, nomcols2.N-3),
    spacing='proportional', ticks=pHints, boundaries=[4.5,5.5,6.5,7.5,8.5,9.5], format='%1i', orientation='horizontal')
    cbaxes2.set_xlabel('Internal pH')
    cbaxes2.xaxis.set_label_position('top')
    cbaxes2.xaxis.tick_top()

    # cmappable = ScalarMappable(norm=Normalize(0,1), cmap=cmapper.r2a())
    # fig.colorbar(clr.BoundaryNorm(pHints, nomcols.N), cax=cbaxes2, label='pH', orientation='horizontal', pad=0.5, aspect=7)

    fig.subplots_adjust(wspace=0.05, hspace=0.07, right=0.995, left=0.125, bottom=0.062, top=0.9)
    plt.savefig('complex_fluxes.pdf')
complex_fluxes()







def pow_perm():

    # make a reactor object
    rtr = Enceladus('EncT', nominals=True, T=298)
    # make an organism, for this example, use the preset methanogen.
    org = TOM(rtr)
    org.maintenance.Tdef='None'
    org.maintenance.pHdef='FluxPerm'

    Perm4 = pH_flux_power(org, 1e-3, 1e-6, 1e-6, org.base_volume)
    Perm8 = pH_flux_power(org, 1e-3, 1e-11, 1e-11, org.base_volume)
    org.locale.env.T = 400
    Perm4_350 = pH_flux_power(org, 1e-3, 1e-6, 1e-6, org.base_volume)
    Perm8_350 = pH_flux_power(org, 1e-3, 1e-11, 1e-11, org.base_volume)

    fig, ax = plt.subplots(figsize=(6,3.25), ncols=2, sharey=True)
    ax[0].fill_between(Perm4['pH'], Perm4['pow'], Perm8['pow'])
    ax[0].set_title('Variable Permeability')
    ax[1].fill_between(Perm8['pH'], Perm8_350['pow'], Perm8['pow'])
    ax[1].fill_between(Perm4['pH'], Perm4_350['pow'], Perm4['pow'], facecolor='g')

    ax[1].set_title('Variable Temperature')
    for a in ax:
        a.set_yscale('log')
        a.set_xlabel('pH')
        a.set_xlim(0,14)
        a.set_ylim(1e-27,1e-9)
        a.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.8)
    ax[0].set_ylabel('Power Demand [W cell$^{-1}$]')

    ax[0].text(10.,10**(-13.25), r'$\bar{P} = 10^{-6}$', c='C0', rotation=40)
    ax[0].text(10.,10**(-20.5), r'$\bar{P} = 10^{-11}$', c='C0', rotation=40)
    ax[1].text(10.,10**(-13), r'$\bar{P} = 10^{-6}$', c='g', rotation=40)
    ax[1].text(10.,10**(-20.5), r'$\bar{P} = 10^{-11}$', c='C0', rotation=40)

    fig.subplots_adjust(wspace=0., right=0.995, top=0.905, bottom=0.14)
    # plt.show()

    plt.savefig('pH_pow_perm.pdf')

pow_perm()
