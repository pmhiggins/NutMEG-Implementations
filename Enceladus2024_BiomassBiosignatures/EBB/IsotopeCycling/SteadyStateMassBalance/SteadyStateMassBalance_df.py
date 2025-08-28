import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import math
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
from FractionationProcess import FractionationProcess as FP

from copy import deepcopy
import sympy as sp

class SteadyStateMassBalance_df(object):
    """
    Class for solving and monitoring the steady state flux in a system.
    For an arbitrary number of processes, monitors the isotope ratios of
    fluxes entering and leaving the system and maintains them in the
    _dict attribute.
    """

    def __init__(self, FP_it):
        _dict = {'f':[], 'R':[], 'a':[], 'R_S':[]}
        _eq = {'f':[], 'R':[], 'a':[], 'R_S':[], 'fR':[]}
        for _FP in FP_it:
            _dict['f'].append(_FP.f)
            _dict['R'].append(_FP.R)
            _dict['a'].append(_FP.a)
            _dict['R_S'].append(_FP.R_S)

            _eq['f'].append(0)
            _eq['R'].append(0)
            _eq['a'].append(0)
            _eq['R_S'].append(0)
            _eq['fR'].append(0)

            for eq, lvl in _FP.equalities.items():
                _eq[eq][-1] = lvl

        if not math.isclose(sum(_dict['f']), 0, abs_tol=5e-15):

            raise ValueError('Bulk fluxes are not conserved! The net is '+str(sum(_dict['f'])))

        self.names = [_FP.name for _FP in FP_it]
        self.cols = ['f', 'R', 'a', 'R_S', 'fR']

        self.df = pd.DataFrame(_dict, index=self.names)
        self.eq = pd.DataFrame(_eq, index=self.names)

        self.df['fR'] = self.df['f']*self.df['R']



    def to_FP_list(self):
        FP_lst = []
        for index, row in self.df.iterrows():
            _eq =  [n for n in list(self.eq) if self.eq[n][index]==True]
            _FP = FP(index, row['f'], row['R'], a=row['a'], R_S=row['R_S'], equalities=_eq)
            FP_lst.append(_FP)
        return FP_lst

    def to_FP_dict(self):
        FP_dict = {}
        for index, row in self.df.iterrows():
            _eq =  [n for n in list(self.eq) if self.eq[n][index]==True]
            _FP = FP(index, row['f'], row['R'], a=row['a'], R_S=row['R_S'], equalities=_eq)
            FP_dict[index] = _FP
        return FP_dict



    def compute_unknown_element(self, _name, _param='R'):
        """
        Solve simultaneous equations to compute _param in the steady state system.
        """

        _df = deepcopy(self.df)

        eq_count = self.eq.max().max()# no. of different equalities


        x = sp.Symbol('x') # let x be the number we are trying to find
        eq_syms = sp.symbols('y0:%d'%eq_count)


        _df.loc[_name, _param] = x

        used_eq_syms = [] # not all the syms will be used, so store the ones that are here
        eq_vals = [] # store their values too, for substitution

        for i, sym in enumerate(eq_syms):

            ues = None
            for cname in _df.columns:

                _keys = _df.loc[self.eq[cname]== i+1].index.values
                if self.eq[_param][_name] == i+1:
                    sym = x

                for _k in _keys:
                    ues = _df[cname][_k]
                    _df.loc[_k, cname] = sym
                    if cname == 'R_S':
                        _df.loc[_k, 'R'] = sym / _df['a'][_k]
                    if cname == 'a':
                        _df[_k, 'R'] = _df['R_S'][_k] / sym

                if sym == x:
                    ues = None


            if ues != None:
                used_eq_syms.append(sym)
                eq_vals.append(ues)
            else:
                used_eq_syms.append(x)
                eq_vals.append(None)


        _df['fR'] = _df['f']*_df['R']
        _df['fF'] = _df['f']*_df['R'] / (1 + _df['R'])

        expr = _df['fF'].sum() # sum as a sympy expression
        _subs = []
        for sym, val in zip(used_eq_syms, eq_vals):
            if sym != x:
                _subs.append(tuple((sym, val)))

        subbed_expr = expr.subs(_subs)
        x_sol =  sp.solve(subbed_expr)
        s2 = np.array(x_sol, dtype=complex) # convert to a numpy array as complex no.s

        x_val = s2[s2.real > 0.][0].real # choose the real, positive solution

        # now resub in all the variables to the df
        for i, sym in enumerate(eq_syms):

            for cname in self.df.columns:

                _keys = _df.loc[self.eq[cname]== i+1].index.values
                if self.eq[_param][_name] == i+1:
                    sym = x_val
                else:
                    sym = eq_vals[i]

                for _k in _keys:

                    _df.loc[_k, cname] = sym
                    if cname == 'R_S':
                        _df.loc[_k, 'R'] = sym / _df['a'][_k]
                    if cname == 'a':
                        _df.loc[_k, 'R'] = _df['R_S'][_k] / sym

        # if there are no equalities, we only need to resub in x
        # if len(eq_syms)==0:
        #     _df[_param][_name] = x_val
        # there was a problem where if the sym solving for was not an eq,
        # but there were other eqs in the system, it was not subbing properly.
        _df[_param][_name] = x_val

        _df['fR'] = _df['f']*_df['R']
        _df = _df.drop('fF', axis=1)

        self.df = _df

        return x_val
