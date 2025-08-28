import warnings
import numpy as np
import os, sys
sys.path.append(os.path.dirname(__file__)+'/../../../')

class MethaneBioFracFactorFixed:

    def __init__(self, a):
        self.a = a

class MethaneBioFracFactorUniform:

    def __init__(self, alims=[1.02,1.08]):
        self.a = random.uniform(alims)


class MethaneBioFracFactorGropp:

    def __init__(self, CO2conc, H2conc):
        """
        Fractionation factor for hydrogenotropic methanogenesis from Gropp et al. 2022
        at 60 C and [CH4]= 1e-5. interpolated from a calculation as a function
        of DeltaG.
        """

        # See MethanogenInterpolation/GroppFits.py for the working behind this.
        # This file is a pickled function, returned by scipy.interpolate.CloughTocher2DInterpolator
        # That class has been depreciated in newer versions of scipy, but works on version 1.3

        interp_func = np.load(os.path.dirname(__file__)+'/../../data/MethanogenFractionation/GroppInterp.npy', allow_pickle=True)[()]

        self.a = interp_func(
          [[np.log10(H2conc)]],[[np.log10(CO2conc)]])[0][0]


        if H2conc < 1e-9 or H2conc > 1e-1:
            warnings.warn('H2 concentration is outside our confidence in fractionation factor')
        if CO2conc < 1e-4 or CO2conc > 1e-1:
            warnings.warn('CO2 concentration is outside our confidence in fractionation factor')
