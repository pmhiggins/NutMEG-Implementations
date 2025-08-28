import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import math
import pandas as pd
import numpy as np
from FractionationProcess import FractionationProcess as FP
from SteadyStateMassBalance_df import SteadyStateMassBalance_df as SSMBdf
from typing import Iterable



class SSMB_builder(SSMBdf):
    """
    Streamlined class for building an SteadyStateMassBalance_df object
    """

    def __init__(self, input, output, inside, solvefor, name='C', mechnames=('in', 'out', ['inside'])):

        if isinstance(inside, dict):
            inside = [inside]

        if solvefor == 'in':
            input['f'] = -output['f']
            for ins in inside:
                input['f'] -= ins['f']
        elif solvefor == 'out':
            output['f'] = -input['f']
            for ins in inside:
                output['f'] -= ins['f']

        _in = FP(name+'_in',
          input.pop('f', 0.),
          R=input.pop('R',0.),
          equalities=input.pop('equalities', {}))

        _out = FP(name+'_out',
          output.pop('f', 0.),
          R=output.pop('R',0.),
          equalities=output.pop('equalities', {}))

        FPs = [_in, _out]

        for i, ins in enumerate(inside):
            _ins = FP(name+'_'+mechnames[2][i],
              ins.pop('f',0.),
              R=ins.pop('R', None),
              a=ins.pop('a', None),
              R_S=ins.pop('R_S',None),
              equalities=ins.pop('equalities', {})) #'a':2

            FPs.append(_ins)

        FPtup = tuple(FPs)

        SSMBdf.__init__(self, FPtup)
        self.solvefor=name+'_'+solvefor


def builder_tester():
    input = {'f':1.0, 'R':0.001}
    loss = [{'f':-0.5, 'a':1.05, 'R_S':0.001}, {'f':-0.1, 'a':1.5, 'R_S':0.001}]
    output = {'f':0., 'R':0.}

    SB = SSMB_builder(input, output, loss, 'out', name='C', mechnames=('in','out','catabolism','anabolism'))
    print(SB.df)
    SB.compute_unknown_element(SB.solvefor)
    print(SB.df)


class CO2_SSMBdf(SSMBdf):

    def __init__(self, input, output, inside, solvefor, name='CO2', mechnames=('in', 'out', ['inside'])):
        if 'equalities' not in input:
            input['equalities'] = {}
            input['equalities']['R'] = 1

        for ins in inside:
            if 'equalities' not in ins:
                ins['equalities'] = {}
                ins['equalities']['R_S'] = 1

        SSMB_builder.__init__(self, input, output, inside, solvefor,
          name='CO2', mechnames=mechnames)



class CH4_SSMBdf(SSMBdf):

    def __init__(self, input, output, inside, solvefor, name='CH4', mechnames=('in', 'out', ['inside'])):

        SSMB_builder.__init__(self, input, output, inside, solvefor,
          name=name, mechnames=mechnames)
