import methanogen_extractor as extractor
import sys,os
sys.path.append(os.path.dirname(__file__)+'../../NutMEG/')
import NutMEG as es
import NutMEG.util.NutMEGparams as nmp

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.colors as clr
import matplotlib as mpl

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cmr10'
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset'] = 'cm'

thisdbpath=nmp.std_dbpath

Trange=np.linspace(275,375, num=21)
nATPs=np.linspace(0.5,1.5, num=11)

# for performing the simulation, remove this docstring.
"""
for T in Trange:
    for n in nATPs:
        es,pt,ps,pg,s = extractor.iterateESynths(E, 'averageMethanogen', paramchange={'Tdef':'None', 'Temp':T, 'n_ATP':n, 'mol_CH4':3e-8}, save='data/EsynthMesh/'+str(3e-8)+'_'+str(n)+'_'+str(int(T))+'.csv', dbpath=thisdbpath)
"""

# mesh color setup
dn, dT = nATPs[1]-nATPs[0], Trange[1]-Trange[0]

# ncd, Tcd = np.mgrid[slice(nATPs[0], nATPs[-1] + dn, dn),
#                 slice(Trange[0], Trange[-1] + dT, dT)]
ncdd, Tcdd = np.mgrid[slice(nATPs[0]-(dn/2), nATPs[-1] + (dn/2), dn),
                slice(Trange[0]-(dT/2), Trange[-1] + (dT/2), dT)]

Esreq = np.zeros_like(ncdd)


# mesh color: synthesis
#colorblind friendly :
clst = [(0.012, 0.386, 0.4), (0.352, 0.767, 0.684), (0.913, 0.883, 0.733), (0.939, 0.625, 0.262), (0.863,0.231,0.112)]

cmap = clr.LinearSegmentedColormap.from_list('custom cbf', clst, N=256)


for i in range(len(nATPs)):#ncdd.shape[0]):
    for j in range(len(Trange)):#ncdd.shape[1]):
        ES, PT, PS, PG, S = extractor.extract_Esynths_csv('data/EsynthMesh/'+str(3e-8)+'_'+str(nATPs[i])+'_'+str(int(Trange[j]))+'.csv')
        print(ES)
        m,c = np.polyfit(ES[:-1],PT[:-1], deg=1)
        Esreq[i,j] = -c/m


plt.pcolormesh(ncdd,Tcdd,Esreq, norm=LogNorm(vmin=Esreq[:-1].min(), vmax=Esreq[:-1].max()), cmap=cmap)
# add in ,shading='gouraud' above to smooth plot.
plt.xlabel(r'ATP per mol $\mathregular{CO}_2}$', fontsize=14)
plt.ylabel('Temperature [K]', fontsize=14)
cb = plt.colorbar(label='Synthesis scaler needed')
cb.set_label(label='Synthesis scaler needed', size=14)
plt.tight_layout()

# plt.savefig('Emesh.eps')
plt.savefig('Emesh.pdf')





# mesh color: constant basal maintenance Energy
"""
Esreq = Esreq[:-1, :-1]
Trange=np.linspace(280,375, num=20)

for i in range(len(nATPs)):#ncdd.shape[0]):
    for j in range(len(Trange)):#ncdd.shape[1]):
        ES, PT, PS, PG, S = extractor.extract_Esynths_csv('data/EsynthMesh/'+str(3e-8)+'_'+str(nATPs[i])+'_'+str(int(Trange[j]))+'.csv')

        m,c = np.polyfit(ES[:-1],PT[:-1], deg=1)
        Esreq[i,j] = m+c
        # print(nATPs[i],Trange[j],m+c)

plt.pcolor(ncdd,Tcdd,Esreq, norm=LogNorm(vmin=Esreq[:-1].min(), vmax=Esreq[:-1].max()), cmap='gist_earth_r')
plt.xlabel(r'ATP per mol $\mathregulardata2}$')
plt.ylabel('Temperature [K]')
plt.colorbar(label='Throttle Power Required')
plt.show()
"""
