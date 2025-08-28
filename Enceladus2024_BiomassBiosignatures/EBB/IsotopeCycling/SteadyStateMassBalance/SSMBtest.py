import sys, os; sys.path.append(os.path.dirname(__file__)+'/../../../../../NutMEG')

from EnceladusSSMBdf import CO2_SSMBdf
from EnceladusSSMBdf import CH4_SSMBdf


import pandas as pd

sys.path.append(os.path.dirname(__file__)+'/../../')
from EncBmBs_utils import IsotopeConversions as ICtools

class SimpleMethaneBox:

    def __init__(self, R_CO2_out, f_CH4, a, R_in_CH4=0., j_in_CH4=0.):

        self.R_CO2_out = R_CO2_out
        self.f_CH4 = f_CH4
        self.a = a
        self.R_in_CH4 = R_in_CH4
        self.j_in_CH4 = j_in_CH4

        # First calculate for CO2. Solve for the input, as the output is passed
        CO2_input = {'f':0., 'R':0.}
        CO2_loss = [{'f':-self.f_CH4, 'a':self.a, 'R_S':0.}]
        CO2_output = {'f':-1+self.f_CH4, 'R':self.R_CO2_out}

        self.CO2_SB = CO2_SSMBdf(CO2_input, CO2_output, CO2_loss, 'in')
        self.CO2_SB.compute_unknown_element(self.CO2_SB.solvefor)

        # First calculate for CO2. Solve for the output, based on the biotic CH4 fraction.
        CH4_input = {'f':self.j_in_CH4, 'R':self.R_in_CH4}
        CH4_inside = [{'f':self.f_CH4, 'a':self.a, 'R_S':self.CO2_SB.df['R']['CO2_in']}]
        CH4_output = {'f':-self.j_in_CH4-self.f_CH4, 'R':0.}

        self.CH4_SB = CH4_SSMBdf(CH4_input, CH4_output, CH4_inside, 'out')

        self.CH4_SB.compute_unknown_element(self.CH4_SB.solvefor)

"""
SMB = SimpleMethaneBox(ICtools.dCtoR(60), 0.001, 1.01)

print(SMB.CO2_SB.df)
print(SMB.CH4_SB.df)
print(ICtools.RtodC(SMB.CO2_SB.df['R']['CO2_in']))
print(ICtools.RtodC(SMB.CH4_SB.df['R']['CH4_inside']))
"""
