import math
import numpy as np

class ChiralTools:

    """ Helper functions for chirality and racemization calculations """

    kdict = {'kala': (24.91, 15863),
      'kasp' : (23.34, 14725),
      'kglu' : (22.56, 15675),
      'kisoleu' : (27.24, 16913),
      'kleu' : (24.3, 15750),
      'kphenal' : (19.82, 13440),
      'kval' : (23.61, 15754)}

    @staticmethod
    def get_k(T, ktype):
        """
        Get the racemization rate constant in [/s] of species
        of type ktype at temperature T.

        Options are: kala, kasp, kglu, kisoleu, kleu, kphenal, kval
        and kLever. Endmembers can be provided by kglu (min) and kLever (max)
        """
        kf = lambda pre, post, T : np.exp(pre-(post/T))
        if ktype=='kLever':
            return (0.00012*np.exp(0.10174*(T-273)))/(365*24*3600)
        else:
            return kf(ChiralTools.kdict[ktype][0], ChiralTools.kdict[ktype][1], T)

    @staticmethod
    def get_Qdl(DL):
        """ Useful factor for DL calcuations """
        return (1+DL)/(1-DL)

    @staticmethod
    def get_DLt(k, t, DL):
        """
        Get the D:L ratio after t seconds with an original value of DL,
        using rate constant k

        cf B. A. Cohen and C. F. Chyba. “Racemization of Meteoritic Amino
        Acids”. Icarus 145.1 (2000), pp. 272–281
        """
        if DL==1:
            return 1.
        _exponent = 2*k*t

        if _exponent > np.log(np.finfo('d').max):
            # if the time is too long such that it'll overflow the exponential,
            # send back a DL of 1.
            return 1.

        eo = np.exp(_exponent)

        top = (eo*ChiralTools.get_Qdl(DL)) - 1
        bottom = (eo*ChiralTools.get_Qdl(DL)) + 1

        return top/bottom

    @staticmethod
    def DL_to_Lf(DL):
        """ Convert D:L ratio to L-form % """
        return 100 / (1+DL)

    @staticmethod
    def DL_to_ee(DL):
        """ Convert D:L ratio to enantiometric excess """
        return math.abs((DL-1)/(DL+1))

    @staticmethod
    def Lf_to_DL(Lf):
        """ Convert L-form % to D:L ratio """
        return (100/Lf)-1.

    @staticmethod
    def Lf_to_ee(Lf):
        """ Convert L-form % to enantiometric excess """
        return ChiralTools.DL_to_ee(ChiralTools.Lf_to_DL(Lf))

    @staticmethod
    def ee_to_Lf(ee):
        """ Convert enantiometric excess to L-form % """
        return 0.5*(1-ee)*100.

    @staticmethod
    def ee_to_DL(ee):
        """ Convert enantiometric excess to DL """
        return ChiralTools.Lf_to_DL(ChiralTools.ee_to_Lf(ee))

    @staticmethod
    def DL_from_DLabio_fbio(DL_abio, f_bio=1.0):
        return (1-f_bio)/(f_bio+(1/DL_abio))

    @staticmethod
    def s_to_yr(s):
        return s*3.17098e-8

    @staticmethod
    def yr_to_s(yr):
        return yr/3.17098e-8

    @staticmethod
    def get_initial_DL(k, t_end, DL_end):
        """
        return the initial D:L ratio for an AA that has been racemizing for
        t_end seconds with rate constant k
        """
        R_t = math.exp(-2*k*t_end) * (1 + DL_end) / (1 - DL_end)
        DL_t = (R_t -1 ) / (1 + R_t)

        return DL_t

    @staticmethod
    def time_to_racemic(T, ktype, steps=1e5, init_DL=0., early_cut=False):

        k = ChiralTools.get_k(T, ktype)

        dt = math.pow(10, 12 - (math.log10(k) +15))

        tr = np.linspace(dt,dt*steps,int(steps))

        Ratios = np.zeros(len(tr))
        Means = np.zeros(len(tr))

        Ratios[:] = np.nan
        Ratios[0] = init_DL
        Means[0] = init_DL

        for ti, t in enumerate(tr[:-1]):
            Ratios[ti+1] = ChiralTools.get_DLt(k, tr[ti+1]-tr[ti], Ratios[ti])
            Means[ti+1] = np.mean(Ratios[:ti+1])
            if early_cut and (np.round(Ratios[ti+1], 4) == 1.0):
                break

        return tr, Ratios, Means
