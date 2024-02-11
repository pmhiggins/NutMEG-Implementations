import matplotlib as mpl
from collections import OrderedDict

class PlotSetup:
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    mpl.rcParams['font.family'] = 'stixgeneral'#'sans-serif'
    mpl.rcParams['font.sans-serif'] = 'cmr10'
    mpl.rcParams['font.serif'] = 'Times.tcc'
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['lines.linewidth'] = 4

    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams['text.usetex'] = False
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    mpl.rc('axes', unicode_minus=False)
    mpl.rcParams['hatch.linewidth'] = 2.0
    # print(mpl.rcParams.keys()) # to see what can be changed
    mpl.rcParams['font.size'] = 12
    # mpl.rcParams['font.weight'] = 'bold'

    mpl_linestyles = OrderedDict(
        [('solid',               (0, ())),
         ('dotted',              (0, (1, 3))),
         ('densely dotted',      (0, (1, 1))),

         ('dashed',              (0, (5, 5))),
         ('densely dashed',      (0, (5, 1))),

         ('dashdotted',          (0, (3, 5, 1, 5))),
         ('densely dashdotted',  (0, (3, 1, 1, 1))),

         ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
         ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1))),

         ('loosely dotted',      (0, (1, 10))),
         ('loosely dashed',      (0, (5, 10))),
         ('loosely dashdotted',  (0, (3, 10, 1, 10))),
         ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10)))
         ])
