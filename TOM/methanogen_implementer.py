import sys, os, ast, math
sys.path.append(os.path.dirname(__file__)+'/../../NutMEG/')
import matplotlib.pyplot as plt
import unique_efficiencies

import NutMEG as es
import NutMEG.util.NutMEGparams as nmp
import NutMEG.reaction as reaction
from NutMEG.applications.NutMEmatcher import NutMEmatcher
import NutMEG.plotter as nutplt
from NutMEG.reactor.saved_systems.VenusDrop import VenusDrop

dbpath= es.util.NutMEGparams.std_dbpath
es.db_helper.create_major_tables(replace=False)

defaultparams = {'Tdef':'None',
  'mol_CO2':0.003,
  'mol_H2':0.0014,
  'MP':0,
  'radius':0.8e-6,
  'k_RTP':((0.4e-12)/(0.003*(0.0014**4)))/(3600*(2**((355.16-298)/10))),
  'n_ATP':1.0,
  'Pconc':0.1,
  'uptake_consts':{'P':1e-10},
  'mol_CH4':1e-9,
  'pH':6.8,
  'lifespan':float('inf'),
  'inputs':{},
  'Temp':355.15}

class efficiencies:
    """Simple class to get the efficiencies of methanogens."""



    def __init__(self, orgname='M_jannaschii', target={'GrowthRate':0.0003}, P_bar=False, *args, **kwargs):
        T = kwargs.pop('Temp',355.15)
        P = kwargs.pop('Pressure', 2e5)
        if P_bar==True:
            self.P_b = 101325
        else:
            self.P_b = P
        self.params = {'Tdef':kwargs.pop('Tdef','None'),
          'mol_CH4':kwargs.pop('mol_CH4',1e-9),
          'MP':kwargs.pop('MP',0),
          'radius':kwargs.pop('radius',0.8e-6),
          'CH4rate':kwargs.pop('CH4rate',0.4e-12/3600),
          'n_ATP':kwargs.pop('n_ATP',1.0),
          'Pconc':kwargs.pop('Pconc',0.1),
          'uptake_consts':kwargs.pop('uptake_consts',{'P':1e-10}),
          'pH':kwargs.pop('pH',6.8),
          'lifespan':kwargs.pop('lifespan',float('inf')),
          'inputs':kwargs.pop('imputs',{}),
          'Temp':T,
          'Pressure':P,
          'shape':kwargs.pop('shape', 'sphere'),
          'mol_CO2':VenusDrop.getgasconc('CO2(aq)', 0.2*P, T, P_bar=P, S=0),
          'mol_H2':VenusDrop.getgasconc('H2(aq)', 0.8*P, T, P_bar=P, S=0),
          'ESfrac': kwargs.pop('ESfrac', 1.0)}

        self.orgname=orgname
        self.target=target
        self.k_RTP = (self.params['CH4rate']/(self.params['mol_CO2']*(self.params['mol_H2']**4)))/((2**((self.params['Temp']-298)/10)))
        # self.k_RTP = (self.params['CH4rate']/(self.params['mol_CO2']*(self.params['mol_H2']**4)))/(3600*(2**((self.params['Temp']-298)/10)))

    @classmethod
    def get_eff(cls, orgname, Temp=None, paramchange={}):
        org_uniqueparams, Gr = None, None

        if orgname=='M_jannaschii':
            org_uniqueparams, GR = unique_efficiencies.M_jannaschii_efficiencies(paramchange=paramchange)
        elif orgname=='averageMethanogen':
            org_uniqueparams, GR = unique_efficiencies.avg_org_efficiencies(Temp=Temp, paramchange=paramchange)
        else:
            raise ValueError('orgname '+orgname+ ' is not recognised')

        return cls(orgname=orgname, target={'GrowthRate':GR}, **org_uniqueparams)



    def setup_methanogenesis(self, R):

        CO2 = reaction.reagent('CO2(aq)', R.env, phase='aq')
        H2aq = reaction.reagent('H2(aq)', R.env, phase='aq')
        CH4aq = reaction.reagent('CH4(g)', R.env, phase='g')
        H2O = reaction.reagent('H2O(l)', R.env, phase='l')

        thermalMG = reaction.reaction({CO2:1, H2aq:4}, {CH4aq:1, H2O:2},
          R.env)


        thermalMG.rate_constant_RTP = self.k_RTP
        # default is the one we calculated base on methanogens (LB3 pg 11)

        return thermalMG


    def initial_conditions(self, R, scale=1.0, NH3conc=0.1, H2Sconc=0.1, setupcomp=True):
        #scale : scaler for concentrations

        if setupcomp:
            mol_CO2 = self.params['mol_CO2']*scale #LB3 pg11
            mol_CH4 = self.params['mol_CH4']
            mol_H2 =  self.params['mol_H2']*scale
            Pconc = self.params['Pconc']

            mol_H=10**(-R.pH)

            # reagents
            CO2 = reaction.reagent('CO2(aq)', R.env, phase='aq', conc=mol_CO2,
              activity=mol_CO2)
            H2aq = reaction.reagent('H2(aq)', R.env, phase='aq', conc=mol_H2,
              activity=mol_H2)
            CH4aq = reaction.reagent('CH4(g)', R.env, phase='g', conc=mol_CH4,
              activity=mol_CH4)
            H2O = reaction.reagent('H2O(l)', R.env, phase='l', conc=55.5,
              phase_ss=True)
            # el = reaction.reagent('e-', self.env, charge=-1)
            H = reaction.reagent('H+', R.env, charge=1, conc=mol_H,
              phase='aq', activity=mol_H)


            R.composition = {CO2.name:CO2, H2aq.name:H2aq,
              CH4aq.name:CH4aq, H2O.name:H2O, H.name:H}

            # add in CHNOPS, in elemental form for now
            # we already have C in the form of CO2 and CH4
            R.composition['NH3(aq)'] = reaction.reagent('NH3(aq)', R.env, phase='aq', conc=NH3conc,
              activity=NH3conc)
            # P and S we dont actually know, so let's arbitrarily say 1 micromole
            R.composition['P(aq)'] = reaction.reagent('P(aq)', R.env, phase='aq', conc=Pconc,
              activity=Pconc, thermo=False)
            R.composition['H2S(aq)'] = reaction.reagent('H2S(aq)', R.env, phase='aq', conc=H2Sconc,
              activity=H2Sconc, thermo=False)


        el = reaction.reagent('e-', R.env, charge=-1)
        # overall
        r = {R.composition['CO2(aq)']:1, R.composition['H2(aq)']:4}
        p = {R.composition['CH4(g)']:1, R.composition['H2O(l)']:2}
        thermaloa = reaction.reaction(r,p,R.env)

        R.add_reaction(thermaloa)

    def getvol(self):
        if self.params['shape'] == 'sphere':
            vol=(4/3)*math.pi*((self.params['radius'])**3)
        elif self.params['shape'] == 'rod':
            vol=(math.pi*(self.params['radius'][0]**2)*self.params['radius'][1])
        drym=vol*300
        mass=vol*1000
        return vol, drym, mass


    def get_efficiency(self, Rname='MM8020', startnum=500, ESfrac=1.0, scale=1.0, dbpath=nmp.std_dbpath):

        R = es.reactor(Rname, workoutID=False, pH=self.params['pH'], dbpath=dbpath)
        R.change_T(self.params['Temp'])
        R.change_P(self.P_b)

        self.initial_conditions(R, scale=scale) # sets up composition, reaction

        volume, dry_mass, mass = self.getvol()

        H = es.horde(self.orgname, R, self.setup_methanogenesis(R), startnum, Tdef=self.params['Tdef'], mass=mass, dry_mass=dry_mass, volume=volume, workoutID=False, n_ATP=self.params['n_ATP'], dbpath=dbpath)

        H.E_synth = H.E_synth*ESfrac

        # horde and reactor are set up, time to fix their IDs

        R.dbh.workoutID()
        H.dbh.workoutID()

        Cu = es.culture(hordes=[H])

        ES = es.ecosystem(R, Cu, dbpath=dbpath)

        NMM = NutMEmatcher(H, R)
        e, S, O = NMM.match(level='NutME', target=self.target)


        PSupply = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S, 'PowerSupply_'+O, dbpath=dbpath)[50][0]
        PThrottle = PSupply*e
        PGrowth = es.ecosystem_dbhelper.db_helper.extract_param_db_Sim(S, 'GrowthPower_'+O, dbpath=dbpath)[50][0]


        return e, PThrottle, PSupply, PGrowth, S, O



    def perform_sim(self, Rname='MM8020', startnum=500, scale=1.0, dbpath=nmp.std_dbpath, dt=None, stoppers={}):
        """ perform a standard simulation """

        R = es.reactor(Rname, workoutID=False, pH=self.params['pH'], composition_inputs=self.params['inputs'], dbpath=dbpath)
        R.change_T(self.params['Temp'])
        R.change_P(self.P_b)

        self.initial_conditions(R, scale=scale) # sets up composition, reaction

        volume, dry_mass, mass = self.getvol()

        H = es.horde(self.orgname, R, self.setup_methanogenesis(R), startnum, Tdef=self.params['Tdef'], mass=mass, dry_mass=dry_mass, volume=volume, Basal=self.params['MP'], n_ATP=self.params['n_ATP'], base_life_span=self.params['lifespan'], workoutID=False, dbpath=dbpath, uptake_consts=self.params['uptake_consts'])

        H.E_synth = H.E_synth*self.params['ESfrac']

        # horde and reactor are set up, time to fix their IDs
        H.dbh.workoutID()
        R.dbh.workoutID()

        Cu = es.culture(hordes=[H])

        ES = es.ecosystem(R, Cu, dbpath=dbpath)
        ES.stoppingdict.update(stoppers)

        ES.predict_growth(dt=dt, tmax=5e7)

        return R.LocID, H.OrgID, ES.dbh.SimID
