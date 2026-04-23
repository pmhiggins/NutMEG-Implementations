import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from copy import deepcopy
from itertools import chain
import numpy as np
import math
from ChiralTools import ChiralTools as CT

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['font.size'] = 7
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['lines.markersize'] = 3
plt.rcParams['lines.linewidth'] =2.
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica']
plt.rcParams['mathtext.fontset'] = 'dejavusans'

plt.rcParams['errorbar.capsize'] = 8
plt.rcParams['lines.linewidth'] = 4
ls = ['dashed', '-', 'dotted']


fig, ax = plt.subplots(figsize=(5.5,2.5), nrows=1)
axs = [ax]
for ax in axs:
    ax.set_xlim(1,9)
axs[0].set_ylim(40,110)

ax = axs[0]


def get_m_eb(ls):
    """ get mean and errorbars from a list for plotting """
    m = (min(ls)+max(ls))/2
    eb = abs(m-min(ls))
    return m, eb

# endmember properties for vertical lines
vert_dict = {
    2 : {'T':273, 'ktype':'kglu', 'Lthres':51, 'ls':'-'},
    4 : {'T':273, 'ktype':'kLever', 'Lthres':51, 'ls':'-'},
    6 : {'T':333, 'ktype':'kglu', 'Lthres':51, 'ls':'-'},
    8 : {'T':333, 'ktype':'kLever', 'Lthres':51, 'ls':'-'},
}

_map = 'viridis'
cm = plt.get_cmap(_map)

tlims = [10**-1, 10**8] # time endmembers
norm = plt.Normalize(vmin=np.log10(tlims[0]), vmax=np.log10(tlims[1]))

for k,v in vert_dict.items():

    tra, Ra, Ma = CT.time_to_racemic(v['T'], ktype=v['ktype'])
    Ra = CT.DL_to_Lf(Ra)
    t51a = np.argmax(Ra<v['Lthres'])
    t99 = np.argmin(Ra>99)
    these_ts = tra[t99:t51a]
    these_ts = CT.s_to_yr(these_ts)
    y = Ra[t99:t51a]
    x = [k,]*len(y)

    norm_these_ts = norm(np.log10(these_ts))

    for i in range(len(y)-1):
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=cm(norm_these_ts[i]), ls=v['ls'])


    # ax.hlines(100, xmin=k-0.3, xmax=k+0.3, color=cm(norm_these_ts[0]))
    ax.scatter([k], 100, color=cm(norm_these_ts[0]), marker='D', s=40., edgecolors='k', zorder=10)

    ax.text(k, 102, r'<10$^{'+str(round(math.log10(these_ts[0]), 1))+'}$ yr', ha='center', va='bottom')
    ax.arrow(k, 53, 0, -0.1,  # Start at (2,4), extends 1 unit in x and 2 in y
          head_width=0.5, head_length=3, overhang=0.8,
          fc=cm(norm_these_ts[-1]), ec=cm(norm_these_ts[-1]), length_includes_head=False, lw=4)
    ax.text(k, 47, r'>10$^{'+str(round(math.log10(these_ts[-1]),1))+'}$ yr', ha='center', va='top')

ax.axhline(50, c='k', ls='dotted')

for _ax in axs:
    _ax.axvline(3, c='k', lw=2)
    _ax.axvline(5, c='k', lw=2)
    _ax.axvline(7, c='k', lw=2)

    _ax.tick_params(axis='x', bottom=False, labelbottom=False)


ax.set_ylabel(r'Biotic amino acid $L$-form %')
ax.set_yticks([50,60,70,80,90,100])

sm = plt.cm.ScalarMappable(cmap=_map, norm=norm)
plt.subplots_adjust(left=0.15, right=0.8, bottom=0.1, top=0.9)
cax = plt.axes([0.85, 0.1, 0.03, 0.8])

fig.colorbar(sm, cax=cax, orientation='vertical', label=r'$\log_{10}$ ( Time [yr] since $L$-form % = 100)')

plt.savefig('figs/Chiral_summary.pdf')
plt.close()
