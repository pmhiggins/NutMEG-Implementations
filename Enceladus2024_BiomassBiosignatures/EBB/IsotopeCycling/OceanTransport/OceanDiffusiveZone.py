from C_1D import C_1D
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc, erfcinv
import mpl_toolkits.axisartist as AA
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
import math
from matplotlib.offsetbox import AnchoredText

mpl.rcParams['lines.linewidth'] = 3



cmap = plt.get_cmap('tab20', 20)

D_CO2_l = C_1D.get_D_CO2(0., 0.01)
D_CO2_h = C_1D.get_D13(D_CO2_l, -0.87) # -0.87 is the epsilon for CO2 diffusion

D_CH4_l = C_1D.get_D_CH4(0.)
D_CH4_h = C_1D.get_D13(D_CH4_l, -3.29) # -3.29 is the epsilon for CH4 diffusion

D_turbulent = 1e-4


Ls = [1]
t_tilde_range = np.logspace(-7,5,num=1000) #-3 for no turbulent


def get_C_tilde(x_tilde, t_tilde, L=20000, D=0.001, t_ref=1e4):
    alpha = t_ref * D / (L**2)
    return erfc(x_tilde / (2* np.sqrt(alpha * t_tilde)))

def get_t(x_tilde, C_tilde, L=20000, D=0.001):
    return ((x_tilde*L)**2)/(4.*D*(erfcinv(C_tilde)**2))


fig, axs = plt.subplots(nrows=2, figsize=(7.5,6.5))


for L in Ls:

    t_diff = L*L / D_CO2_l


    C_tilde_CO2_h = get_C_tilde(1., t_tilde_range, L=L, D=D_CO2_h, t_ref=t_diff)
    C_tilde_CO2_l = get_C_tilde(1., t_tilde_range, L=L, D=D_CO2_l, t_ref=t_diff)

    C_tilde_CH4_h = get_C_tilde(1., t_tilde_range, L=L, D=D_CH4_h, t_ref=t_diff)
    C_tilde_CH4_l = get_C_tilde(1., t_tilde_range, L=L, D=D_CH4_l, t_ref=t_diff)

    C_tilde_turbulent = get_C_tilde(1., t_tilde_range, L=L, D=D_turbulent, t_ref=t_diff)
    C_tilde_turbulent_CO2_h = get_C_tilde(1., t_tilde_range, L=L, D=D_turbulent+D_CO2_h, t_ref=t_diff)
    C_tilde_turbulent_CO2_l = get_C_tilde(1., t_tilde_range, L=L, D=D_turbulent+D_CO2_l, t_ref=t_diff)

    C_tilde_turbulent_CH4_h = get_C_tilde(1., t_tilde_range, L=L, D=D_turbulent+D_CH4_h, t_ref=t_diff)
    C_tilde_turbulent_CH4_l = get_C_tilde(1., t_tilde_range, L=L, D=D_turbulent+D_CH4_l, t_ref=t_diff)


    t_range = t_tilde_range*t_diff/ (3600*24*365)

    print(D_CO2_l)
    axs[0].plot(t_range  / (L**2), C_tilde_CO2_l, label=r'CO$_2$; $D_l=5.85\times10^{-10}$ m$^{2}$s$^{-1}$', c=cmap(2), alpha=0.8)
    axs[0].plot(t_range  / (L**2), C_tilde_CH4_l, label=r'CH$_4$; $D_l = 8.75\times10^{-10}$ m$^{2}$s$^{-1}$', c=cmap(8))
    axs[0].plot(t_range  / (L**2), C_tilde_turbulent, label=r'Turbulent $D_e=1\times10^{-4}$ m$^{2}$s$^{-1}$', c=cmap(13))

    axs[0].set_ylabel(r"Molality at $x_3$ / Molality at $x_2$") # ; $\tilde{C}_{x_{4}}
    axs[0].set_xlabel(r"Time / $(L_2)^2$ [yr m$^{-2}$]")

    axs[1].plot(t_range / (L**2), 1000 * abs((C_tilde_CO2_h / C_tilde_CO2_l) - 1), label=r'CO$_2$ molecular diffusion;'+'\n'+r'$D_{CO2, l}=5.85\times10^{-10}$ m$^{2}$s$^{-1}$', c=cmap(2), alpha=0.8)
    axs[1].plot(t_range / (L**2), 1000 * abs((C_tilde_CH4_h / C_tilde_CH4_l) - 1), label=r'CH$_4$ molecular diffusion;'+'\n'+r'$D_{CH4, l}=8.75\times10^{-10}$ m$^{2}$s$^{-1}$', c=cmap(8))

    axs[1].plot(t_range / (L**2), 1000 * abs((C_tilde_turbulent_CO2_h / C_tilde_turbulent_CO2_l) - 1), label=r'Turbulent diffusion +'+'\n'+'CO$_{2}$ molecular diffusion;'+'\n'+'$D_t + D_{CO2, j}$', c=cmap(13))
    axs[1].plot(t_range / (L**2), 1000 * abs((C_tilde_turbulent_CH4_h / C_tilde_turbulent_CH4_l) - 1), label=r'Turbulent diffusion +'+'\n'+'CH$_{4}$ molecular diffusion;'+'\n'+'$D_t + D_{CH4, j}$', c=cmap(13), ls='dashed')

    axs[1].set_ylabel(r"Approximate $-\Delta\delta$C$_{x2/x3}$ $[\perthousand]$")
    axs[1].set_xlabel(r"Time / $(L_2)^2$ [yr m$^{-2}$]")

for ax in axs:
    ax.set_xlim(t_range[0],t_range[-1])
    ax.set_xscale('log')
    ax.axvline(1e9 / (20000*20000), c='tab:red', ls='dashed', label='1Gyr at $L_2$ = 20 km')
    ax.axvline(1e9 / (1000*1000), c='tab:green', ls='dashdot', label='1Gyr at $L_2$ = 1 km')



axs[1].set_yscale('log')
axs[1].axhline(1., c='slategray', lw=3, alpha=0.5)

plt.tight_layout()
plt.subplots_adjust(right=0.65)
axs[1].legend(loc='center right', bbox_to_anchor=(1.65, 1.1), labelspacing = 2)

plt.savefig('PureDiffusion.png', dpi=300)
plt.savefig('PureDiffusion.pdf')
plt.close()




fig, axs = plt.subplots(nrows=2, figsize=(5,8))

fig = plt.figure(figsize=(7.5,4))

gs1 = GridSpec(1, 2, width_ratios=[8, 1], wspace=0.00)
axs = [fig.add_subplot(gs1[0]), fig.add_subplot(gs1[1])]


Ls = np.logspace(1,math.log10(20000))

t_10_CO2_l = get_t(1., 0.1, L=Ls, D=D_CO2_l)
t_10_CO2_h = get_t(1., 0.1, L=Ls, D=D_CO2_h)
t_10_CH4_l = get_t(1., 0.1, L=Ls, D=D_CH4_l)
t_10_CH4_h = get_t(1., 0.1, L=Ls, D=D_turbulent)

t_10_turbulent = get_t(1., 0.1, L=Ls, D=1e-4)

t_90_CO2_l = get_t(1., 0.9, L=Ls, D=D_CO2_l)
t_90_CO2_h = get_t(1., 0.9, L=Ls, D=D_CO2_h)
t_90_CH4_l = get_t(1., 0.9, L=Ls, D=D_CH4_l)
t_90_CH4_h = get_t(1., 0.9, L=Ls, D=D_CH4_h)

t_90_turbulent = get_t(1., 0.9, L=Ls, D=D_turbulent)

axs[0].plot(Ls, t_10_CO2_l/(3600*24*365), label=r'CO$_2$ molecular diffusion;'+'\n'+r'$\tilde{C}_{x3}=0.1$', c=cmap(2), alpha=0.8)
axs[0].plot(Ls, t_10_CH4_l/(3600*24*365), label=r'CH$_4$ molecular diffusion;'+'\n'+r'$\tilde{C}_{x3}=0.1$', c=cmap(8))
axs[0].plot(Ls, t_10_turbulent/(3600*24*365), label=r'Turbulent diffusion only;'+'\n'+r'$\tilde{C}_{x3}=0.1$', c=cmap(13))


axs[0].plot(Ls, t_90_CO2_l/(3600*24*365), label=r'CO$_2$ molecular diffusion;'+'\n'+r'$\tilde{C}_{x3}=0.9$', c=cmap(2), ls='dashed', alpha=0.8)
axs[0].plot(Ls, t_90_CH4_l/(3600*24*365), label=r'CH$_{4}$ molecular diffusion;'+'\n'+r'$\tilde{C}_{x3}=0.9$', c=cmap(8), ls='dashed')
axs[0].plot(Ls, t_90_turbulent/(3600*24*365), label=r'Turbulent diffusion only;'+'\n'+r'$\tilde{C}_{x3}=0.9$', c=cmap(13), ls='dashed')


axs[0].set_ylabel(r'Time required to reach $\tilde{C}_{x3}$ threshold [yr]')
axs[0].set_xlabel(r'Depth of diffusive layer $L_2=x_3 - x_2$ [m]')


axs[1].plot(Ls, 1000 * abs((get_C_tilde(1.0, 1.0, L=Ls, D=D_CO2_h, t_ref=t_10_CO2_l) / get_C_tilde(1.0, 1.0, L=Ls, D=D_CO2_l, t_ref=t_10_CO2_l)) - 1), label=r'CO$_2$, $\tilde{C}^{L}=0.1$', c=cmap(2), alpha=0.8)
axs[1].plot(Ls, 1000 * abs((get_C_tilde(1.0, 1.0, L=Ls, D=D_CH4_h, t_ref=t_10_CO2_l) / get_C_tilde(1.0, 1.0, L=Ls, D=D_CH4_l, t_ref=t_10_CO2_l)) - 1), label=r'CH$_4$, $\tilde{C}^{L}=0.1$', c=cmap(8))
axs[1].plot(Ls, 1000 * abs((get_C_tilde(1.0, 1.0, L=Ls, D=D_CO2_h, t_ref=t_90_CO2_l) / get_C_tilde(1.0, 1.0, L=Ls, D=D_CO2_l, t_ref=t_90_CO2_l)) - 1), label=r'CO$_2$, $\tilde{C}^{L}=0.9$', c=cmap(2), ls='dashed', alpha=0.8)
axs[1].plot(Ls, 1000 * abs((get_C_tilde(1.0, 1.0, L=Ls, D=D_CH4_h, t_ref=t_90_CO2_l) / get_C_tilde(1.0, 1.0, L=Ls, D=D_CH4_l, t_ref=t_90_CO2_l)) - 1), label=r'CH$_4$, $\tilde{C}^{L}=0.9$', c=cmap(8), ls='dashed')

axs[1].yaxis.set_label_position("right")
axs[1].yaxis.tick_right()
axs[1].set_ylabel(r"Approximate $\Delta\delta$C$_{x2/x3}$ $[\perthousand]$")
axs[1].set_ylim(1e-2,10)

for ax in axs:
    ax.set_yscale('log')
    ax.set_xscale('log')

axs[1].set_xticks([])
axs[1].set_xticklabels([])
axs[0].legend(loc='center right', bbox_to_anchor=(2.05,0.5), labelspacing=2)
plt.tight_layout()
plt.subplots_adjust(right=0.6)
plt.savefig('PureDiffusionthresholds.png', dpi=300)
plt.savefig('PureDiffusionthresholds.pdf', dpi=300)
plt.close()
