import NutMEG as nm
from NutMEG.culture.saved_organisms.TypicalOptimalMethanogen import TypicalOptimalMethanogen as TOM

import numpy as np

class TOMFetcher:
    """
    Class for generating a NutMEG.saved_organism.TypicalOptimalMethanogen
    object, adjusting and extracting its properties for analysis.
    """

    @staticmethod
    def fetch(rtr, n_ATP='best', k_corr=0.0, Tdef='TOM', cap_mr=True):
        """
        Return a TOM instance in a reactor-like rtr.

        n_ATP : 'best' or float
            ATP yield per mole of CO2 metabolised. 'best' will calculate the
            most power efficient yield in the host reactor.
        k_corr : float
            Order-of-magnitude adjustment to microbial metabolic rate constant.
        Tdef : string
            Identifier for maintenance power
        cap_mr : boolean
            Whether the metabolic rate should be capped to the empirical values
            that were fitted in Higgins & Cockell (2020)
        """

        if n_ATP == 'best':
            _T = TOM(rtr, workoutID=False, fromdata=False, n_ATP=1.,
              paramchange={'TOMpressure':182000, 'Tdef':Tdef}, k_corr=k_corr)
            # reassign n_ATP for the TOM we return to the most efficient value
            n_ATP = TOMFetcher.get_best_nATP(rtr, _T, k_corr)

        if cap_mr:
            return TOM(rtr, workoutID=False, fromdata=False, n_ATP=n_ATP,
              paramchange={'TOMpressure':182000, 'Tdef':Tdef}, k_corr=k_corr)
        else:
            return TOM(rtr, workoutID=False, fromdata=False, n_ATP=n_ATP,
              paramchange={'TOMpressure':182000, 'Tdef':Tdef, 'maxmet':1e20},
              k_corr=k_corr)



    @staticmethod
    def get_best_nATP(rtr, org, k_corr):

        """
        Find the molar yield of ATP that maximises power supply for a given
        organism-like object org in a given reactor rtr.

        The organism's respiration should already be setup for conditions in rtr.

        """

        best_nATP = 0.0
        best_PS = 0.0

        #free energies in terms of / RT for ease
        GA = org.respiration.G_A/(8.31*rtr.env.T)
        GP = org.respiration.G_P/(8.31*rtr.env.T)
        GCarr = np.linspace(0.25*GP,-GA, num=1000)
        PSarr = GCarr-(GCarr*np.exp((GA+GCarr)))

        best_GC = GCarr[PSarr.argmax()]*(8.31*rtr.env.T)
        return best_GC / org.respiration.G_P


    @staticmethod
    def MGparams(_TOM):
        """
        Get instantaneous Quotient, Free energy, ATP energy, Power Supply, and
        maximum microbial metabolic rate constant k from a TOM
        """
        _TOM.respiration.net_pathway.update_molar_gibbs_from_quotient()
        Q = _TOM.respiration.net_pathway.quotient
        G = _TOM.respiration.net_pathway.molar_gibbs
        GP = _TOM.respiration.G_P

        PS = _TOM.get_supplied_power()
        mk = ((_TOM.max_metabolic_rate /(_TOM.locale.composition['CO2(aq)'].activity*(_TOM.locale.composition['H2(aq)'].activity**4))) /
        (2**((TOM.locale.env.T-298)/10))) - _TOM.respiration.net_pathway.rate_constant_RTP

        return Q, G, GP, PS, mk

    @staticmethod
    def getPS(_TOM):
        """ Get the power supply to a TOM"""
        _TOM.respiration.net_pathway.update_molar_gibbs_from_quotient()
        PS = _TOM.get_supplied_power()
        return PS
