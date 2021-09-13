import sys, os, ast, math
sys.path.append(os.path.dirname(__file__)+'../../NutMEG')

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clr

import numpy as np
import pandas as pd
import NutMEG as es # is this even used??


mpl.rcParams['contour.negative_linestyle'] = 'solid'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cmr10'
mpl.rcParams['font.serif'] = 'Times.tcc'
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['lines.linewidth'] = 4

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rc('axes', unicode_minus=False)
mpl.rcParams['hatch.linewidth'] = 2.0
# print(mpl.rcParams.keys()) # to see what can be changed
mpl.rcParams['font.size'] = 12
# mpl.rcParams['font.weight'] = 'bold'
