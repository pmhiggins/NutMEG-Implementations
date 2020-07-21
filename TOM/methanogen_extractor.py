import sys, os, ast, re
sys.path.append(os.path.dirname(__file__)+'../../NutMEG/')
from methanogen_implementer import efficiencies
import unique_efficiencies as ueff
import math, statistics

import pandas as pd
from pandas import DataFrame, read_csv
import NutMEG as es
import NutMEG.util.NutMEGparams as nmp


def allmethanogens_tocsv(filename='output.csv', nATP=1.0, ESfrac=1.0, scale=1.0, fixed_size=None, mol_CH4=1e-9, dbpath=nmp.std_dbpath):
    """Save the output from efficiencies.get_efficiency to a csv file for
    each of the empirical methanogens.
    """
    resdict={}

    file = os.path.dirname(__file__)+'/data/methanogens.csv'
    df = pd.read_csv(file, header=0)
    reslst=[]
    for index, row in df.iterrows():
        try:
            width = (float(row['Min. cell width'])+float(row['Max. cell width']))/2
            length = (float(row['Min. cell length'])+float(row['Max. cell length']))/2
            T = (float(row['Min. optimal growth temp.'])+float(row['Max. optimal growth temp.']))/2
            pH = (float(row['Min. optimal growth pH'])+float(row['Max. optimal growth pH']))/2


            eff = efficiencies(orgname=row['Name'].replace(' ',''), radius=[width*0.5e-6,length*1e-6], shape='rod', Temp=273+T, pH=pH, n_ATP=nATP, CH4rate=ueff.CH4Rate_from_GrowthRate(row['Growth rate']/3600, CH4conc=mol_CH4), target={'GrowthRate':row['Growth rate']/3600}, Pressure=float(row['Pressure'])*1000, mol_CH4=mol_CH4)

            if fixed_size != None:
                # correct the organism's size to fised_size
                eff.params['shape'] = 'sphere'
                eff.params['radius'] = fixed_size
                v = (4/3)*math.pi*(fixed_size**3)
                eff.k_RTP = (eff.params['CH4rate']/(eff.params['mol_CO2']*(eff.params['mol_H2']**4)))/((2**((eff.params['Temp']-298)/10)))

            thislst=[]
            thislst.append(row['Name'].replace(' ','_'))

            MF, PThrottle, PSupply, PGrowth, S, O = eff.get_efficiency(startnum=500, ESfrac=1.0, scale=scale, dbpath=dbpath)
            thislst.append(MF)
            thislst.append(PThrottle)
            thislst.append(PSupply)
            thislst.append(PGrowth)
            thislst.append(S)
            thislst.append(O)
            thislst.append(273+T)
            thislst.append(pH)
            thislst.append(eff.getvol()[0])

            reslst.append(thislst)

        except Exception as e:
            raise e
            continue


    resdf = pd.DataFrame(reslst, columns=['name', 'MF', 'PThrottle', 'PSupply', 'PGrowth', 'SimID', 'OrgID', 'T', 'pH', 'vol'])

    resdf.to_csv(filename)

def allmethanogens_fromcsv(filename='output.csv', extra=False):
    """Extract the efficiencies.get_efficiency data for the empirical
     methanogens from a .csv file generated using allmethanogens_tocsv above.
     """
    resdf2 = pd.read_csv(filename)
    ret = {'MF':[], 'PThrottle':[], 'PSupply':[], 'PGrowth':[], 'SimID':[], 'OrgID':[], 'Temp':[], 'vol':[], 'source':[], 'nut':[]}
#E, PT, PS, PG, SimID, OrgID, Temp, v, source, nut = [], [], [], [], [], [], [], [], [], []
    for index, row in resdf2.iterrows():
        ret['MF'].append(row['MF'])
        ret['PThrottle'].append(row['PThrottle'])
        ret['PSupply'].append(row['PSupply'])
        ret['PGrowth'].append(row['PGrowth'])
        ret['SimID'].append(row['SimID'])
        ret['OrgID'].append(row['OrgID'])
        ret['Temp'].append(row['T'])
        if extra:
            est = get_extra_params(index)
            ret['vol'].append(row['vol'])
            ret['source'].append(est['GR source'])
            ret['nut'].append(est['nutrients'])

    return ret


def get_extra_params(index):
    """helper function for allmethanogens_fromcsv. Extract the growth rate
    source and nutrients elements from the original methanogen csv. """
    file = os.path.dirname(__file__)+'/data/methanogens.csv'
    df = pd.read_csv(file, header=0)
    ret = {}
    ret['GR source'] = df['GR source'][index]
    ret['nutrients'] = df['Extra nutrients?'][index]

    return ret




def iterateESynths(ESlist, orgname, paramchange={}, save=None, dbpath=nmp.std_dbpath):
    """Using the Esynthfrac argument for eff.get_efficiency, work out what maintenance
    requirement would be needed if synthesis were more/less expensive for the scalers in
    ESlist.
    If save is passed as a str, save to a csv there.
    """
    E, PT, PG, PS, S = [], [], [], [], []
    for e in ESlist:
        eff = efficiencies.get_eff(orgname, Temp=paramchange.get('Temp', 298), paramchange=paramchange)

        ef, pt, ps, pg, s, o = eff.get_efficiency(startnum=500, ESfrac=e, scale=1.0, dbpath=dbpath)
        PT.append(pt)
        PG.append(pg)
        PS.append(ps)
        E.append(e)
        S.append(s)
        if ef == 0.:
            break
    if save != None:
        res=[]
        for e, ps, pt, pg, s in zip(E, PS, PT, PG, S):
            res.append([e, ps, pt, pg, s])
        resdf = pd.DataFrame(res, columns=['E', 'PS', 'PT', 'PG', 'SimID'])
        resdf.to_csv(save)
    return E, PT, PG, PS, S


def extract_Esynths_csv(filename):
    """Extract and return saved data from iterateESynths."""
    resdf2 = pd.read_csv(filename)
    E, PT, PS, PG, S = [], [], [], [],[]
    for index, row in resdf2.iterrows():
        E.append(row['E'])
        PT.append(row['PT'])
        PS.append(row['PS'])
        PG.append(row['PG'])
        S.append(row['SimID'])
    return E,PT,PS,PG,S
