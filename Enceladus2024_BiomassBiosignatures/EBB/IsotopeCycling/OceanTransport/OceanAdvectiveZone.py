import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from C_1D import C_1D
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd



class OceanAdvectiveSim:
    """
    Class for defining and performing a simple isotope-specific 1D ocean
    transport model that considers advection, diffusion, and lateral loss
    processes.
    """

    def __init__(self, L_km, Q_kg_s, vfunc, v_b2t, Afunc, A_b2t, Atop, c_CO2_top, c_CH4_top, f_loss, f_loss_range, dir=''):

        self.L_km = L_km
        self.Q_kg_s = Q_kg_s
        self.CH4_flux = self.Q_kg_s*c_CH4_top
        self.CO2_flux = self.Q_kg_s*c_CO2_top


        self.vfunc_str = vfunc
        if vfunc == 'v_func_s':
            self.vfunc = self.v_func_s
        elif vfunc == 'v_func_linear':
            self.vfunc = self.v_func_linear
        self.v_b2t = v_b2t

        self.Afunc_str = Afunc
        if Afunc == 'Area_lateral':
            self.Afunc = self.Area_lateral
        elif Afunc == 'Area_const':
            self.Afunc = self.Area_const
        self.A_b2t = A_b2t
        self.Atop = Atop
        self.c_CO2_top = c_CO2_top
        self.c_CH4_top = c_CH4_top
        self.f_loss = f_loss
        self.f_loss_range = f_loss_range

        self.set_fn(dir = dir)

    def set_fn(self, dir=''):

        self.profile_fn = dir
        for i in [self.L_km, self.Q_kg_s, self.vfunc_str, self.v_b2t, self.Afunc_str, self.A_b2t, self.Atop, self.c_CO2_top, self.c_CH4_top, self.f_loss, self.f_loss_range[0], self.f_loss_range[1]]:
            self.profile_fn += str(i)+'_'

        self.profile_fn = self.profile_fn[:-1] +'.csv'

    @staticmethod
    def v_func_s(xs, v0_to_vt):
        vy_max = 1  # Starting at upper asymptote y=1
        vy_min = vy_max/(v0_to_vt)  # Lower asymptote
        k = 0.0001  # Steepness factor, negative for downward slope
        x0 = 0.55*xs[-1]  # Midpoint, where the function starts dropping noticeably
        vs_dimless = (vy_max - vy_min) / (1 + np.exp(k * (xs - x0))) + vy_min

        return vs_dimless

    @staticmethod
    def v_func_linear(xs, v0_to_vt):
        vs_dimless = np.array([1. for i in xs])
        return vs_dimless

    @staticmethod
    def Area_lateral(xs, top, bottom):
        # Area_top_est = 100*10. * 500*10000. # 100 m x 500 km; in dm.
        # Area_bot_est = Area_top_est *1e-2
        A_sqrt = np.linspace(math.sqrt(bottom), math.sqrt(top), num=len(xs))
        Area_est = A_sqrt**2
        return Area_est

    @staticmethod
    def Area_const(xs,top,bottom):
        if top!=bottom:
            raise ValueError('constant area has different top and bottom!')
        return np.array([top for i in xs])

    @staticmethod
    def R_to_dC(R):
        return ((R/0.0112372)-1)*1000




    def implement(self, R_CO2_out, R_CH4_out, save=True):

        CO2 = C_1D(
          C_1D.get_D_CO2(0., 0.01)* 100, # bulk D in dm3
          -0.87, # epsilon for diffusion
          'out', self.CO2_flux, self.c_CO2_top, self.f_loss,
          transitprops={'v_func':self.vfunc, 'v_b2t':self.v_b2t, 'A_func':self.Afunc, 'A_b2t':self.A_b2t, 'A_top':self.Atop},
          xprops = {'dx':0.01, 'L_km': self.L_km},
          fluxprops = {'j_loss_xrange':self.f_loss_range})
          # 0.0106753399


        CH4 = C_1D(
          C_1D.get_D_CH4(0.)* 100, # bulk D in dm3
          -3.29, # epsilon for diffusion
          'out', self.CH4_flux, self.c_CH4_top, self.f_loss,
          transitprops={'v_func':self.vfunc, 'v_b2t':self.v_b2t, 'A_func':self.Afunc, 'A_b2t':self.A_b2t, 'A_top':self.Atop},
          xprops = {'dx':0.01, 'L_km': self.L_km},
          fluxprops = {'j_loss_xrange':self.f_loss_range})

        C_CO2_init = CO2.guess_C_profile_advective()
        C_CO2_l_results = []
        C_CO2_h_results = []

        C_CH4_init = CH4.guess_C_profile_advective()
        C_CH4_l_results = []
        C_CH4_h_results = []

        max_D = max(CO2.D_bulk, CO2.D_heavy, CH4.D_bulk, CH4.D_heavy)

        # vs[0] will always be the largest velocity
        default_dt = (CO2.xprops['dx_dm']**2) / (CO2.vs[0]*CO2.xprops['dx_dm']+ 2*(max_D))
        adv_timescale = CO2.xprops['L_dm'] / CO2.vs[0]


        dts = [default_dt, adv_timescale/100]

        for dt in dts:

            _CO2 = deepcopy(CO2)
            _CH4 = deepcopy(CH4)

            ### CO2
            C_CO2_l_init = deepcopy(C_CO2_init)
            C_CO2_h_init = deepcopy(C_CO2_init)*R_CO2_out

            _CO2.reset_default_fluxes(R=1.)
            C_CO2_l_SS = _CO2.integrate_depth_diff_only(C_CO2_l_init, restype='fin', dt=dt)
            C_CO2_l_results.append(C_CO2_l_SS)

            _CO2.reset_default_fluxes(R=R_CO2_out)
            C_CO2_h_SS = _CO2.integrate_depth_diff_only(C_CO2_h_init, D=CO2.D_heavy, restype='fin', dt=dt)
            C_CO2_h_results.append(C_CO2_h_SS)

            ### CH4
            C_CH4_l_init = deepcopy(C_CH4_init)
            C_CH4_h_init = deepcopy(C_CH4_init)*R_CH4_out

            _CH4.reset_default_fluxes()
            C_CH4_l_SS = _CH4.integrate_depth_diff_only(C_CH4_l_init, restype='fin', dt=dt)
            C_CH4_l_results.append(C_CH4_l_SS)

            _CH4.reset_default_fluxes(R=R_CH4_out)
            C_CH4_h_SS = _CH4.integrate_depth_diff_only(C_CH4_h_init, D=CH4.D_heavy, restype='fin', dt=dt)
            C_CH4_h_results.append(C_CH4_h_SS)

        this_df = pd.DataFrame.from_dict({'x':CO2.xprops['xs_dm']})

        this_df['CO2_bulk_start'] = C_CO2_init
        this_df['CO2_bulk_end_default_dt'] = C_CO2_l_results[0]
        this_df['CO2_bulk_end_adv_dt'] = C_CO2_l_results[1]

        this_df['CO2_h_end_default_dt'] = C_CO2_h_results[0]
        this_df['CO2_h_end_adv_dt'] = C_CO2_h_results[1]

        this_df['CH4_bulk_start'] = C_CH4_init
        this_df['CH4_bulk_end_default_dt'] = C_CH4_l_results[0]
        this_df['CH4_bulk_end_adv_dt'] = C_CH4_l_results[1]

        this_df['CH4_h_end_default_dt'] = C_CH4_h_results[0]
        this_df['CH4_h_end_adv_dt'] = C_CH4_h_results[1]

        if save:
            this_df.to_csv(self.profile_fn)
        return this_df



    def get_DdC(self, species, dt_types=['default_dt', 'adv_dt'], df=None):

        _out = df
        if type(_out) == type(None):
            _out = pd.read_csv(self.profile_fn)

        res = {}
        for dt_type in dt_types:

            bulk_top = np.array(_out[species+'_bulk_end_'+dt_type])[-1]
            bulk_bottom = np.array(_out[species+'_bulk_end_'+dt_type])[0]
            C13_top = np.array(_out[species+'_h_end_'+dt_type])[-1]
            C13_bottom = np.array(_out[species+'_h_end_'+dt_type])[0]

            if any(x < 0 for x in (bulk_top, bulk_bottom, C13_top, C13_bottom)):
                print("At least one variable is negative.")
                res[t_frac] = np.nan
            else:
                try:
                    R_top = C13_top / bulk_top
                    R_bottom = C13_bottom / bulk_bottom

                    dC_top = OceanAdvectiveSim.R_to_dC(R_top)
                    dC_bottom = OceanAdvectiveSim.R_to_dC(R_bottom)

                    res[dt_type] = dC_bottom - dC_top
                except:
                    print('Problem')
                    res[dt_type] = np.nan

        return res


def Advective_iterate(fn='AdvectiveParameterSpace'):

    dir = os.path.dirname(__file__)+'/../../data/OceanTransport/'


    Area_top_big = 100*10. * 500*10000. # 100 m x 500 km; in dm. # from N+I 2016
    Area_top_small = 1*10. * 100*10000. # 1 m x 100 km; in dm. hypothetical minimum

    iterables = {'Area_funcs':['Area_lateral', 'Area_const'],
      'Area_top':[Area_top_big, (Area_top_big + Area_top_small) / 2, Area_top_small],
      'Area_bottom_top_ratio':[1., 0.5, 0.001],
      'vfuncs':['v_func_s', 'v_func_linear'],
      'v0_to_vt':[1., 10., 100.],
      'f_loss' : [0.05, 0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95],
      'f_loss_range': [(0.4,0.6), (0.5,0.9)],
      'c_CO2_top' : [1e-3],#, 5e-3, 1e-4], # ~ pH 8,9 from Higgins et al 2024 - can expand to salinity-related unc if this has a measurable effect
      'c_CH4_top' : [1e-3 * 0.002/0.0055],#, 5e-3 * 0.002/0.0055, 1e-4 * 0.002/0.0055], # nominals from waite et al 2017
      'J_top_H2O' : [0.1,1,10,100,1000],
      'L':[20.,40.,60.]}

    pspace = {'Afunc':[],
      'Area_top':[],
      'A0_to_At':[],
      'vfunc':[],
      'v0_to_vt':[],
      'f_loss' : [],
      'f_loss_range': [],
      'c_CO2_top' : [],
      'c_CH4_top' : [],
      'J_top_H2O' : [],
      'L':[],
      'DdC (bottom-top) CO2 default_dt' : [],
      'DdC (bottom-top) CO2 adv_dt' : [],
      'DdC (bottom-top) CH4 default_dt' : [],
      'DdC (bottom-top) CH4 adv_dt' :[]}

    Cs_list = []

    i= 0
    n=1
    for L in iterables['L']:
        for Q in iterables['J_top_H2O']:
            for c_CO2, c_CH4 in zip(iterables['c_CO2_top'], iterables['c_CH4_top']):
                for f_loss_range in iterables['f_loss_range']:
                    for f_loss in iterables['f_loss']:
                        for v_b2t in iterables['v0_to_vt']:
                            for v_func in iterables['vfuncs']:
                                if not (v_b2t == 1. and v_func == 'v_func_s'):
                                    for A_b2t in iterables['Area_bottom_top_ratio']:
                                        for A_top in iterables['Area_top']:
                                            for A_func in iterables['Area_funcs']:
                                                if not (A_b2t==1. and A_func == 'Area_lateral'):
                                                    if not (A_b2t!=1. and A_func == 'Area_const'):
                                                        i += 1
                                                        # if i/100. == 1.*n:
                                                        #     print(i)
                                                        #     n += 1
                                                        _Cs = OceanAdvectiveSim(L, Q, v_func, v_b2t, A_func, A_b2t, A_top, c_CO2, c_CH4, f_loss, f_loss_range, dir=dir+'/')
                                                        df = _Cs.implement(0.0106753399, 0.0106753399, save=False)

                                                        pspace['L'].append(L)
                                                        pspace['J_top_H2O'].append(Q)
                                                        pspace['c_CO2_top'].append(c_CO2)
                                                        pspace['c_CH4_top'].append(c_CO2)
                                                        pspace['f_loss_range'].append(f_loss_range)
                                                        pspace['f_loss'].append(f_loss)
                                                        pspace['v0_to_vt'].append(v_b2t)
                                                        pspace['vfunc'].append(v_func)
                                                        pspace['A0_to_At'].append(A_b2t)
                                                        pspace['Area_top'].append(A_top)
                                                        pspace['Afunc'].append(A_func)

                                                        DdC_dict_CO2 = _Cs.get_DdC('CO2', df=df)
                                                        DdC_dict_CH4 = _Cs.get_DdC('CH4', df=df)


                                                        pspace['DdC (bottom-top) CO2 default_dt'].append(DdC_dict_CO2['default_dt'])
                                                        pspace['DdC (bottom-top) CO2 adv_dt'].append(DdC_dict_CO2['adv_dt'])

                                                        pspace['DdC (bottom-top) CH4 default_dt'].append(DdC_dict_CH4['default_dt'])
                                                        pspace['DdC (bottom-top) CH4 adv_dt'].append(DdC_dict_CH4['adv_dt'])
    pspace = pd.DataFrame(pspace)
    pspace.to_csv(dir+fn+'-summary.csv')



def Pe_vs_f_loss_contours(fn='AdvectiveParameterSpace', vmin=-3, vmax=0.):

    dir = os.path.dirname(__file__)+'/../../data/OceanTransport/'
    df = pd.read_csv(dir+fn+'-summary.csv')

    df['log(vt)'] = np.log10(df['J_top_H2O']*df['c_CO2_top']/(df['Area_top']*df['c_CO2_top']))
    df['log(Pe)'] = np.log10(10000 * df['L']) + df['log(vt)'] - np.log10(C_1D.get_D_CO2(0.,0.01)*100)

    Enc_df = df.where(df['J_top_H2O'] >=100.)

    fig, axs = plt.subplots(ncols=2,nrows=1, figsize=(8,4.5))

    this_axs = 0
    cont=None

    yparam = 'log(Pe)'
    xparam = 'f_loss'
    combo = [xparam, yparam]
    fancyspecies = ['CO$_{2}$', 'CH$_{4}$']

    for i, species in enumerate(['CO2', 'CH4']):

        dt = 'default_dt'
        grouped = df.groupby(combo)['DdC (bottom-top) '+species+' '+dt].apply(list)  # Group delta_C values

        x,y,z = [],[],[]
        for group, values in grouped.items():
            x.append(group[0])
            y.append(group[1])
            z.append(np.nanmax(np.abs(values))+1e-50)

        levels = np.linspace(vmin, vmax,num=13)
        cont = axs[i].tricontourf(x,y,np.log10(z), levels=levels, cmap='viridis', vmin=vmin, vmax=vmax, extend='both')

        axs[i].set_ylabel(yparam)
        axs[i].set_xlabel('Fraction of fluid flow lost on transit $f_{loss}$')
        axs[i].set_title(fancyspecies[i])
        axs[i].axhline(Enc_df['log(Pe)'].min(), c='r', lw=3)
        axs[i].text(0.95, 5., 'Min. Enceladus log(Pe)', c='r', va='bottom', ha='right', fontweight='bold')


    plt.tight_layout()
    fig.subplots_adjust(bottom=0.3)
    cax = fig.add_axes([0.1, 0.12, 0.9, 0.05])
    fig.colorbar(cont, cax=cax, orientation='horizontal')
    cax.set_xlabel(r'log$_{10}$[ - $\Delta \delta^{13}$C$_{x1/x2}$ $[\perthousand]$] = log$_{10}$[ $-(\delta^{13}$C$_{x1}$ - $\delta^{13}$C$_{x2})$]')

    plt.savefig('Pe_vs_f_loss_contours_default.png', dpi=300)
    plt.savefig('Pe_vs_f_loss_contours_default.pdf')

    plt.close()




Advective_iterate()
Pe_vs_f_loss_contours()
