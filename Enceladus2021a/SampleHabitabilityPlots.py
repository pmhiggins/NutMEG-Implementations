import sys, os, ast
sys.path.append(os.path.dirname(__file__)+'../../NutMEG')
sys.path.append(os.path.dirname(__file__)+'../../NutMEG-Implementations/TOM')

import EnceladusPlotStyles as EPS
import SamplingPlotStyles as SPS
import GetSamples as GS
from colormapping import cmapper

"""nominal power supply + maintenance plots"""

# SPS.nominal2x3(save='figs/nominal2x3_n10.pdf')
# SPS.nominal2x3(save='figs/nominal2x3_n05.pdf', nATPchoice=0.5)
# SPS.nominal2x3(save='figs/nominal2x3_n15.pdf', nATPchoice=1.5)
# SPS.nominal2x3(save='suppfigs/nominal2x3_n025.pdf', nATPchoice=0.25)
# SPS.nominal2x3(save='suppfigs/nominal2x3_n200.pdf', nATPchoice=2.0)


"""individual variance plots"""
# SPS.ones_varianceplot(salttype='nom', show=False, cm=cmapper.b2a(), samplesize=1000)
# SPS.ones_varianceplot(salttype='high', show=False, cm=cmapper.b2a(), samplesize=1000)
# SPS.ones_varianceplot(salttype='low', show=False, cm=cmapper.b2a(), samplesize=1000)


"""Total variance plots"""
## use the below to extend the sample by a set amount as needed
# GS.varianceexmaple_extend(2000, 5000, 1000)
## use the below to plot the total variance
GS.varianceexample(2000)


def ind_v_plots():
    """supplemental data plots contribution of the independent variables)"""
    SPS.ones_varianceplot(T=275, salttype='nom', show=False, cm=cmapper.b2a(), samplesize=1000)
    SPS.ones_varianceplot(T=275, salttype='high', show=False, cm=cmapper.b2a(), samplesize=1000)
    SPS.ones_varianceplot(T=275, salttype='low', show=False, cm=cmapper.b2a(), samplesize=1000)

    SPS.ones_varianceplot(T=325, salttype='nom', show=False, cm=cmapper.b2a(), samplesize=1000)
    SPS.ones_varianceplot(T=325, salttype='high', show=False, cm=cmapper.b2a(), samplesize=1000)
    SPS.ones_varianceplot(T=325, salttype='low', show=False, cm=cmapper.b2a(), samplesize=1000)

    SPS.ones_varianceplot(T=375, salttype='nom', show=False, cm=cmapper.b2a(), samplesize=1000)
    SPS.ones_varianceplot(T=375, salttype='high', show=False, cm=cmapper.b2a(), samplesize=1000)
    SPS.ones_varianceplot(T=375, salttype='low', show=False, cm=cmapper.b2a(), samplesize=1000)


    SPS.ones_varianceplot(T=None, pH=8, salttype='nom', show=False, cm=cmapper.b2a(), samplesize=1000)
    SPS.ones_varianceplot(T=None, pH=8, salttype='high', show=False, cm=cmapper.b2a(), samplesize=1000)
    SPS.ones_varianceplot(T=None, pH=8, salttype='low', show=False, cm=cmapper.b2a(), samplesize=1000)

    SPS.ones_varianceplot(T=None, pH=9, salttype='nom', show=False, cm=cmapper.b2a(), samplesize=1000)
    SPS.ones_varianceplot(T=None, pH=9, salttype='high', show=False, cm=cmapper.b2a(), samplesize=1000)
    SPS.ones_varianceplot(T=None, pH=9, salttype='low', show=False, cm=cmapper.b2a(), samplesize=1000)

    SPS.ones_varianceplot(T=None, pH=10, salttype='nom', show=False, cm=cmapper.b2a(), samplesize=1000)
    SPS.ones_varianceplot(T=None, pH=10, salttype='high', show=False, cm=cmapper.b2a(), samplesize=1000)
    SPS.ones_varianceplot(T=None, pH=10, salttype='low', show=False, cm=cmapper.b2a(), samplesize=1000)

# ind_v_plots()
