from NutMEG.culture.saved_organisms.TypicalOptimalMethanogen import TypicalOptimalMethanogen as TOM

def getE_TOM(Enc, n_ATP=1.0, k_corr=0.0):
    """ Get a TOM instance in Enceladus-like Enc. """
    return TOM(Enc, workoutID=False, fromdata=False, n_ATP=n_ATP,
      paramchange={'TOMpressure':182000}, k_corr=k_corr)


def MGparams(QE, k_corr=0., nATP=1.0):
    """
    Get instantaneous Quotient, Free energy, ATP energy, Power Supply, and
    maximum k from a TOM in reactor QE.
    """
    MG = getE_TOM(QE, n_ATP=nATP, k_corr=k_corr)
    MG.respiration.net_pathway.update_molar_gibbs_from_quotient()
    Q = MG.respiration.net_pathway.quotient
    G = MG.respiration.net_pathway.molar_gibbs
    GP = MG.respiration.G_P

    PS = MG.get_supplied_power()
    mk = ((MG.max_metabolic_rate /(QE.composition['CO2(aq)'].activity*(QE.composition['H2(aq)'].activity**4))) /
    (2**((QE.env.T-298)/10))) - MG.respiration.net_pathway.rate_constant_RTP

    return Q, G, GP, PS, mk

def getPS(Enc, nATP, k_corr):
    """ Get the power supply to a TOM in reactor Enc."""
    MG = getE_TOM(Enc, n_ATP=nATP, k_corr=k_corr)
    MG.respiration.net_pathway.update_molar_gibbs_from_quotient()
    PS = MG.get_supplied_power()
    return PS
