import sys, os, math
from statistics import mean
from copy import deepcopy
import pandas as pd
import numpy as np

#location of csv file describing methanogens
_m_file = os.path.dirname(__file__)+'/data/methanogens.csv'

# Use M. jannaschii as the default organism
_defaultparams = {'Tdef':'Tijhuis',
  'MP':0,
  'radius':0.5e-6,
  'CH4rate':0.4e-12/3600,
  'shape':'sphere',
  'n_ATP':1.0,
  'Pconc':0.1,
  'uptake_consts':{'P':1e-10},
  'mol_CH4':1e-9,
  'pH':6.8,
  'lifespan':float('inf'),
  'inputs':{},
  'Temp':355.15,
  'Pressure':200000}

def M_jannaschii_efficiencies(paramchange={}):
    Mj_uniqueparams = {'Tdef':'Tijhuis',
      'MP':0,
      'radius':0.5e-6,
      'CH4rate':0.4e-12/3600,
      'shape':'sphere',
      'n_ATP':1.0,
      'Pconc':0.1,
      'uptake_consts':{'P':1e-10},
      'mol_CH4':1e-9,
      'pH':6.8,
      'lifespan':float('inf'),
      'inputs':{},
      'Temp':355.15,
      'Pressure':200000}
    org_uniqueparams = deepcopy(_defaultparams)
    org_uniqueparams.update(Mj_uniqueparams)
    org_uniqueparams.update(paramchange)
    # growth rate from Ver Eecke 2013
    return org_uniqueparams, 1.19/3600


def CH4Rate_from_GrowthRate(growthrate, CH4conc=1e-8):
    """Estimate the rate of CH4 production using the growth
    rate. ref Powell 1983."""
    return CH4conc*growthrate*(math.exp(growthrate) - 1)


def avg_org_params(paramchange={}):
    """get the 'average' methanogen parameters from the methanogens csv."""
    df = pd.read_csv(_m_file, header=0)
    widths, lengths, Ts, pHs, CH4s, Pressures, GRs = [],[],[],[],[],[],[]
    for index, row in df.iterrows():
        try:
            widths.append((float(row['Min. cell width'])+float(row['Max. cell width']))/2)
            lengths.append((float(row['Min. cell length'])+float(row['Max. cell length']))/2)
            Ts.append(273.15+((float(row['Min. optimal growth temp.'])+float(row['Max. optimal growth temp.']))/2))
            pHs.append((float(row['Min. optimal growth pH'])+float(row['Max. optimal growth pH']))/2)
            Pressures.append(float(row['Pressure'])*1000)
            GRs.append(row['Growth rate']/3600)
            CH4s.append(CH4Rate_from_GrowthRate(row['Growth rate']/3600, CH4conc=paramchange.get('mol_CH4', 3e-8)))
        except Exception as e:
            print(str(e))
            continue

    CH4s = np.polyfit(Ts, np.log(CH4s), 1)
    GRs = np.polyfit(Ts, np.log(GRs), 1)
    return (mean(widths), mean(lengths), Ts, mean(pHs), CH4s, mean(Pressures), GRs)



def avg_org_efficiencies(Temp, paramchange={}):
    """Get the parameters of the TOM at the passed temperature. """
    width, length, T, pH, CH4, Pressure, GR = avg_org_params(paramchange)
    GRT = math.exp(GR[0]*(Temp)+GR[1])
    CH4T = math.exp(CH4[0]*(Temp)+CH4[1])
    avg_org_uniqueparams = {'Tdef':'Tijhuis',
      'MP':0,
      'radius':[width*0.5e-6,length*1e-6],
      'CH4rate':CH4T,
      'shape':'rod',
      'n_ATP':1.0,
      'Pconc':0.1,
      'uptake_consts':{'P':1e-10},
      'mol_CH4':paramchange.get('mol_CH4', 3e-8),
      'pH':pH,
      'lifespan':float('inf'),
      'inputs':{},
      'Temp':Temp,
      'Pressure':Pressure}

    org_uniqueparams = deepcopy(_defaultparams)
    org_uniqueparams.update(avg_org_uniqueparams)
    org_uniqueparams.update(paramchange)
    return org_uniqueparams, GRT
