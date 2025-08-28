import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl

from copy import deepcopy

"""

"""


def R_to_dC(R):
    return ((R/0.0112372)-1)*1000

class C_1D:
    """
    Class for performing a 1D simulation of advective/diffusive processes.
    """

    def __init__(self,
      D_bulk, D_epsilon, in_out, flux, conc, f_loss,
      transitprops={'v_func':None, 'v_b2t':100., 'A_func':None, 'A_b2t':0.01, 'Atop':5e9},
      xprops = {'dx':0.01, 'L_km': 60},
      fluxprops = {'j_loss_xrange':(0.4,0.6)}):


        """
        transitprops: physical constraints on the ocean profile
            {'v_func': function to return dimensionless velocity profile. Must take x, bottom_to_top ratio
            'v_b2t' : float, ratio between advective velocity at bottom vs top
            'A_func' : function to return area profile. Must take x, area at top, area at bottom.
            'A_b2t' : float, ratio between plume area at bottom vs top.
            'Atop' : float, Area of ocean intersecting with escaping region in dm2
            }
        xprops: for x dimensionality
            {'dx': dimensionless step in x, e.g., 0.01,
            'L' : total distance to model
            }
        fluxprops : for flux-related quantities
            { 'f_ice' : fraction of fluid flux to top of ocean that is lost before mass transfer into plume
            { 'f_loss' : fraction of fluid flux lost during ocean transit}

        """


        self.xprops = xprops
        self.xprops['xdim'] = int(1 + (1/self.xprops['dx'])) # number of entries in xs
        self.xprops['L_km'] = xprops['L_km']
        self.xprops['L_dm'] = self.xprops['L_km'] * 10000
        self.xprops['xs_dm'] = np.linspace(0., self.xprops['L_dm'], num = self.xprops['xdim'])
        self.xprops['dx_dm'] = self.xprops['xs_dm'][1]-self.xprops['xs_dm'][0]


        # area trend
        self.Areas = transitprops['A_func'](self.xprops['xs_dm'], transitprops['A_top'], transitprops['A_top']*transitprops['A_b2t'])

        # dimensionless velocity profile (max 1, at ocean floor)
        # note x axis scale is and should be between 0 and L, not 0 and 1.
        self.vs_dimless = transitprops['v_func'](self.xprops['xs_dm'], transitprops['v_b2t'])


        self.transitprops = transitprops
        self.fluxprops = fluxprops

        self.in_out = in_out
        if in_out == 'out':
            self.fluxprops['C_out'] = conc
            self.fluxprops['J_out'] = flux
            self.fluxprops['J_in'] = flux/(1.0-f_loss)
            self.v_top = self.fluxprops['J_out'] / (self.Areas[-1] * self.fluxprops['C_out'])

        elif in_out == 'in':
            self.fluxprops['C_in'] = conc
            self.fluxprops['J_out'] = flux - (f_loss*flux)
            self.fluxprops['J_in'] = flux
            self.v_top = (self.fluxprops['J_out'] / (self.Areas[0] * self.fluxprops['C_in'] * self.transitprops['v_b2t']))

        else:
            raise ValueError('option '+str(in_out)+' not valid as in_out! Must be in or out!')


        # set dimensionalised velocity profile
        # we multiply by v_top / dimless v_top to ensure that the final entry
        # in the vs array matches v_top.
        self.vs = self.vs_dimless * self.v_top / self.vs_dimless[-1]


        # finalise ocean concentrations at top/bottom
        self.fluxprops['C_in'] = self.fluxprops['J_in'] / (self.Areas[0] * self.vs[0])
        self.fluxprops['C_out'] = self.fluxprops['J_out'] / (self.Areas[-1] * self.vs[-1])

        # set misc properties
        self.fluxprops['f_loss'] = f_loss
        self.fluxprops['J_loss'] = self.fluxprops['J_in'] * self.fluxprops['f_loss']
        self.D_bulk = D_bulk
        self.D_epsilon = D_epsilon
        self.D_heavy = C_1D.get_D13(self.D_bulk, self.D_epsilon)

        # setup j_loss
        self.fluxprops['j_loss_bulk'], self.fluxprops['cumulative_J_loss'] = self.setup_stepwise_j_loss_bulk()


    def reset_default_fluxes(self, R=1.0):
        """
        Reset the molar fluxes to original values the object was created with.
        Useful for resetting after a simulation.
        Do not use this in a loop when R != 1.
        """

        if self.in_out == 'out':
            self.fluxprops['J_out'] = R * self.fluxprops['J_out']
            self.fluxprops['J_in'] = self.fluxprops['J_out']/(1.0-self.fluxprops['f_loss'])
        elif self.in_out == 'in':
            self.fluxprops['J_in'] = R * self.fluxprops['J_in']
            self.fluxprops['J_out'] = self.fluxprops['J_in']*(1.0 - self.fluxprops['f_loss'])

        self.fluxprops['J_loss'] = self.fluxprops['J_in'] * self.fluxprops['f_loss']
        self.fluxprops['j_loss_bulk'], self.fluxprops['cumulative_J_loss'] = self.setup_stepwise_j_loss_bulk()


    @staticmethod
    def get_D13(D_bulk, D_epsilon):

        aD  = D_epsilon/1000+1 # convert to frac. factor
        return D_bulk * aD

    @staticmethod
    def get_D_CO2(T, Salinity):

        A1=18.157948
        A2=0.05736
        A3=0.068700205361
        A4=0.000387610239534634
        A5=0.820561458
        A6=1.46331515077

        T=T+273.15; # Celcius input to Kelvin
        _DCO2=-A1*math.exp(-A2*Salinity)+(A3*T)-(A4*(Salinity**(A5))*(T**(A6)))
        _DCO2=_DCO2*1e-9 # m2 / s
        return _DCO2

    @staticmethod
    def get_D_CH4(T):

        T=T+273.15 # Celcius input to Kelvin

        A=-1670
        B=5.98

        _DCH4 = (A/T)+B
        _DCH4 = math.exp(_DCH4)
        _DCH4 = _DCH4*1e-5 # cm2/s
        _DCH4 = _DCH4*1e-4 # m2/s
        return _DCH4


    def setup_stepwise_j_loss_bulk(self, debug=False):

        j_loss = np.zeros(self.xprops['xdim'])
        cumulative_J_loss = np.zeros(self.xprops['xdim'])
        _sumJloss = 0.

        j_loss_start = int(self.fluxprops['j_loss_xrange'][0]*self.xprops['xdim'])
        j_loss_stop = int(self.fluxprops['j_loss_xrange'][1]*self.xprops['xdim'])

        J_loss_step = self.fluxprops['J_loss'] / (j_loss_stop-j_loss_start) # mol/s to be lost in each dx - note this is not changing with area!


        for _ in range(j_loss_start, j_loss_stop):

            j_loss[_] = J_loss_step / (self.Areas[_]*self.xprops['dx_dm'])

            _sumJloss += J_loss_step
            cumulative_J_loss[_] = _sumJloss

        cumulative_J_loss[_:] = _sumJloss

        if debug:
            # check that your total moles lost is correct:
            J_loss_check = np.sum(j_loss*self.Areas*self.xprops['dx_dm'])

        return j_loss, cumulative_J_loss




    def guess_C_profile_advective(self):
        """
        Guess a steady state concentration profile dominated by advective transport
        """
        return (self.fluxprops['J_in'] - self.fluxprops['cumulative_J_loss']) / (self.vs * self.Areas)




    def integrate_depth_diff_only(self, Cs_init, D=None, dt=None, dt_max=100, t_max=3e15, restype = 'fin'):
        """

        Currently does NOT consider and effects of changing area. If that
        affects your velocities (and concentrations!?), you should factor that
        into the vs parameter, or update this method.
        """

        dx = self.xprops['dx_dm']

        _C = np.copy(Cs_init)

        if D == None:
            D = self.D_bulk

        if dt==None:
            # set dt as if advection were important
            dt = (dx**2) / (abs(self.vs[0])*dx + 2*(D))
            # dt = (dx**2) / (2*(D))


        n10=0.05 # for dubugging to save a Cx array each 10% of tmax
        # dt_max = 1000 # maximum number of timesteps
        res=[Cs_init]


        for _k, n in enumerate(np.linspace(0.,1, num=int(dt_max))):

            C_new = np.copy(_C)  # To store the updated values

            for i in range(1, self.xprops['xdim']-2):
                # Advection (central difference)
                # dCdx = (_C[i] - _C[i-1]) / (dx)

                # Diffusion (central difference)
                dCdx2 = (_C[i+1] - 2*_C[i] + _C[i-1]) / ((dx)**2)

                C_new[i] = (D*dCdx2)*(dt) + _C[i]

                if C_new[i] < 0.:
                    C_new[i:] = 0.
                    if restype == 'lst':
                        res.append(np.copy(C_new))
                        return res
                    elif restype == 'fin':
                        return C_new

            # Enforce boundary conditions
            # C_new = self.post_integrate_clarity(C_new, _C)[1]
            C_new[-1] = _C[-1]
            C_new[0] = _C[0]

            # let advection take over and dominate the new profile,
            # then correcting C_top based on the changing "loss" flux due to diffusion
            C_new = self.adv_sweep(C_new, _C, debug=False)

            if np.allclose(C_new, _C, rtol=1e-11, atol=1e-11):

                # print ('equal at n='+str(_k))
                if restype == 'lst':
                    res.append(np.copy(C_new))
                    return res
                elif restype == 'fin':
                    return C_new

            # Update the concentration array
            _C = C_new

        # print('dt_max reached; returning final C profile NOT at SS')

        if restype == 'lst':
            res.append(C_new)
            return res
        elif restype == 'fin':
            return C_new


    def adv_sweep(self, C_new, C_old, debug=False):

        self.flux_correct(C_new, C_old, debug=debug)

        C_corrected = self.guess_C_profile_advective()

        return C_corrected

    def flux_correct(self, C_new, C_old, debug=False):

        cumulative_J_loss = np.zeros(self.xprops['xdim'])
        _sumJloss = 0.

        for _ in range(self.xprops['xdim']):

            _sumJloss += self.fluxprops['j_loss_bulk'][_]*self.Areas[_]*self.xprops['dx_dm']*C_new[_]/C_old[_]
            cumulative_J_loss[_] = _sumJloss

        if debug:
            print('old: ', self.fluxprops['cumulative_J_loss'][-1])
            print('new: ', cumulative_J_loss[-1])

        # update flux at top to relfect changes in moles lost during j_loss
        self.fluxprops['cumulative_J_loss'] = cumulative_J_loss
        if self.in_out == 'out':
            self.fluxprops['J_in'] =  _sumJloss + self.fluxprops['J_out']
        elif self.in_out == 'in':
            self.fluxprops['J_out'] = self.fluxprops['J_in'] - _sumJloss
