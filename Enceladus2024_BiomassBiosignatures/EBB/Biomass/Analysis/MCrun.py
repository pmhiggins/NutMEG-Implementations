import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from BMDRTO_MC import BMDRTO_MC
from BMDRTO_utils import BMDRTO_utils
import numpy as np
import math


"""
Between pH 8 and 9, for main text habitability assesment
requires spec_dr = 'spec_T_273-473_pH_8-9'
"""

Tdefs = ['Lever2pc', 'TOM']
model = 'pitzerPHREEQCnoGases'
spec_dr = 'spec_T_273-473_pH_8-9'
pHfloats = np.linspace(8.,9.,num=6)
pHvals = np.round(pHfloats, decimals=1)
Tvals = np.linspace(273.15, 393.15, num=25)
Clvals = np.logspace(math.log10(0.05),math.log10(0.2), num=21)
Clvals = np.round(Clvals, decimals=5)
GasScale = 1.0


MC = BMDRTO_MC([model], Tdefs, pHvals, Tvals, Clvals,
  GasScale, 'pH8to9_test', cap_mr=False, spec_dr=spec_dr, cores=6)
# this could take a long time. consider lowering the number a lot for testing.
MC.save_df(25000)


GasScale = 10.0
MC = BMDRTO_MC([model], Tdefs, pHvals, Tvals, Clvals,
  GasScale, 'pH8to9_test', cap_mr=False, spec_dr=spec_dr, cores=6)

# this could take a long time. consider lowering the number a lot for testing.
MC.save_df(25000)



"""
Between pH 7 and 10, for supplemental plots and tables
requires spec_dr = 'spec_T_273-473_pH_7-12'
"""

# Tdefs = ['Lever2pc']#, 'TOM']
# model = 'pitzerPHREEQCnoGases'
# spec_dr = 'spec_T_273-473_pH_7-12'
# pHvals = [7.0,7.5,8.0,8.5,9.0,9.5,10.0]
# Tvals = np.linspace(273.15, 393.15, num=13)
# Clvals = np.logspace(math.log10(0.05),math.log10(0.2), num=21)
# Clvals = np.round(Clvals, decimals=5)
# GasScale = 1.0
#
# MC = BMDRTO_MC([model], Tdefs, pHvals, Tvals, Clvals,
#   GasScale, 'pH7to10_test', cap_mr=False, spec_dr=spec_dr, cores=6)
# this could take a long time. consider lowering the number a lot for testing.
# MC.save_df(25000)
#
#
# GasScale = 10.0
# MC = BMDRTO_MC([model], Tdefs, pHvals, Tvals, Clvals,
#   GasScale, 'pH7to10_test', cap_mr=False, spec_dr=spec_dr, cores=6)
# this could take a long time. consider lowering the number a lot for testing.
# MC.save_df(25000)
