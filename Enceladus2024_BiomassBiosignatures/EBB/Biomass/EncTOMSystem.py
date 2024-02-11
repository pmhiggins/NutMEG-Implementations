import sys, os
sys.path.append(os.path.dirname(__file__)+'/../')
from copy import deepcopy
import numpy as np

from EncBmBs_utils import EnceladusFetcher
from EncBmBs_utils import TOMFetcher

class EncTOMSystem:


    def __init__(self, JH2, EncParams={}, TOMParams={}, fluxes_in=None):
        """

        Parameters
        ----------
        JH2 : float
            Total biological consumption of H2 to consider in mol/s
        EncParams : dict
            Dictionary of arguments for the EnceladusFetcher.fetch() method
        TOMParams : dict
            Dictionary of arguments for the TOMFetcher.fetch() method
        fluxes_in : Nonetype or dict
            pass as None for no inward fluxes to the Enc. Alternatively pass a
            dict of the molar input per second using chemical species as the keys.
        """

        _EncParams = deepcopy(EncParams)
        _Enc0Params = deepcopy(EncParams)

        _TOMParams = deepcopy(TOMParams)

        self.fluxes_in = fluxes_in

        self.netfluxes = {'H2(aq)':JH2,
          'CO2(aq)':JH2/4,
          'Methane(aq)':-JH2/4}

        if self.fluxes_in != None:
            self.fluxes_out = {k:self.fluxes_in[k]-v for k, v in self.netfluxes.items()}
            # later functions may add, for example, biomass flows to this
        else:
            # if there is no input, only output is CH4
            # (nonsense in reality, but included here for mass balance)
            # later functions may add, for example, biomass flows to this
            self.fluxes_out = {'Methane(aq)':-JH2/4}


        self.Enc = EnceladusFetcher().fetch(
          _EncParams.pop('T', 273.15),
          _EncParams.pop('pH_bo', 8.0),
          _EncParams.pop('saltlevel', 'nom'),
          P = _EncParams.pop('P', 1),
          model = _EncParams.pop('CarbonatesModel', 'pitzerPHREEQCnoGases'),
          spec_dr = _EncParams.pop('spec_dr','spec_T_273-473_pH_7-12'),
          fixGases=_EncParams.pop('fixGases', True),
          gases = _EncParams.pop('gases','nom'),
          GasScale = _EncParams.pop('GasScale',1.0))

        self.Enc0 = EnceladusFetcher().fetch(
          273.15,
          _Enc0Params.pop('pH_bo', 8.0),
          _Enc0Params.pop('saltlevel', 'nom'),
          P = _Enc0Params.pop('P', 1),
          model = _Enc0Params.pop('CarbonatesModel', 'pitzerPHREEQCnoGases'),
          spec_dr = _Enc0Params.pop('spec_dr','spec_T_273-473_pH_7-12'),
          fixGases=_Enc0Params.pop('fixGases', True),
          gases = _Enc0Params.pop('gases','default'))

        # for a time-integrated simulation, this may be needed.
        # self.Enc.composition_inputs = self.netfluxes


        self.TOM = TOMFetcher().fetch(
          self.Enc,
          _TOMParams.pop('n_ATP','best'),
          _TOMParams.pop('k_corr',0.0),
          _TOMParams.pop('Tdef', 'TOM'),
          _TOMParams.pop('cap_mr', True))



    @classmethod
    def fromEncTOM(Enc, TOM):
        """
        Create an EncTOMSystem with predefined Enceladus-like and TOM-like
        objects.
        """
        self.Enc = Enc
        self.TOM = TOM

    def find_maxhab(self):
        """
        Find the maximum maintenance power, from a suite of estimates
        (TOM, Tijuis+1993, Lever+2015 2pc, Lever+2015 10pc) that is survivable
        by the TOM and assign its identifier to self. habitability_level
        """
        TOM_Tdef = deepcopy(self.TOM.maintenance.Tdef)
        _PS = self.TOM.get_supplied_power()
        hab = 'TOM'
        for Tdef in ['TOM', 'TijhuisAnaerobe', 'Lever2pc', 'Lever10pc']:
            self.TOM.maintenance.Tdef = Tdef
            self.TOM.maintenance.get_P_T()
            _PT = self.TOM.maintenance.net_dict['T']
            if _PS > _PT:
                hab = Tdef
                break
        # reset TOM's Tdef to its original
        self.TOM.maintenance.Tdef = TOM_Tdef
        self.TOM.maintenance.get_P_T()
        # set the habitability level
        self.habitability_level = hab

    def setup_isotopes(self, R_CO2_plume):
        """ To be added in future """
        ## build relevant SSMBdf objects
        return

        

    def update_steadystate_BMDRTO(self):
        """
        Using the current state of the Enceladus and TOM system, calculate the
        steady state biomass, turnover, and cell-specific death rate needed to
        maintain the system. Assumes no nutrient limitation and that biomass
        assembly is not the rate-limiting step of growth (i.e., instead,
         energy uptake is).
        """
        if self.TOM.respiration.rate ==0.:
            self.ssBM = 1
            self.ssDR = 0.
            self.ssTO = 0.
        else:
            self.ssBM = int((0.25*self.netfluxes['H2(aq)'])/self.TOM.respiration.rate)
            self.TOM.num = self.ssBM
            if self.ssBM > 1:
                __TOM = deepcopy(self.TOM)

                # total cell-specific power supply
                PS = __TOM.get_supplied_power(update_energetics=True)

                # update fractional cell-specific power loss owing to maintenance
                __TOM.maintenance.update_P_loss(PS)
                PM_frac = __TOM.maintenance.P_loss

                # fraction of new cells that can be made with the power supply
                # assuming that biomass synthesis is slower than energy uptake
                # and no nutrient limitation
                newcellspercell = max([0 , (PS*(1-PM_frac)) /__TOM.E_synth])

                # steady state turnover is this amound multiplied by the total
                # steady stats population
                self.ssTO = newcellspercell * self.ssBM

                # associated cell-specific 'death rate' is equivalent to growth rate
                # in  steady state
                self.ssDR = ((self.ssBM*newcellspercell) / self.ssBM)

            else:
                self.ssBM = 1
                self.ssDR = 0.
                self.ssTO=0.
