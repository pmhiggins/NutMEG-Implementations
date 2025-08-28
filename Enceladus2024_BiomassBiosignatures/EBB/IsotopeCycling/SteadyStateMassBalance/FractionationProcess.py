import math
import numpy as np

class FractionationProcess(object):

    def __init__(self, name, f, R=None, a=None, R_S=None, equalities={}, **kwargs):
        """
        f is the flux across this process
        isfwd is a bool for direction: True for forward (production) False for reverse (consumption).
        R is the fractionation associated with it
        R_S is the fractionation of the source
        a is the equilibrium alpha of the process
        """
        self.name = name
        self.f = f
        self.R = R

        if R!= None and R_S!=None:
            self.a = R_S / R
            self.R_S = R_S
        elif a!=None and R!=None:
            self.R_S = a*R
            self.a=a
        elif a!=None and R_S!=None:
            self.R = R_S / a
            self.a = a
            self.R_S = R_S
        else:
            self.a = np.nan
            self.R_S = np.nan

        self.equalities = equalities
