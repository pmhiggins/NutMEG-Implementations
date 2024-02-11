from EncSpecSaver import EncSpecSaver as ESSaver
import numpy as np
import math


""" use these to save speciation output across T, pH, and salt level """


""" Complete speciation across parameter space with 2 models """
def broad_pH():
    Cl_low = np.logspace(math.log10(0.05),math.log10(0.2), num=21)
    Cl_low = np.round(Cl_low, decimals=5)
    Cl_high = np.logspace(math.log10(0.2),math.log10(0.4), num=11)
    Cl_high = np.round(Cl_high, decimals=5)

    pHfloats = [7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0]
    Tfloats = np.linspace(273.15, 473.15, num=41)
    print(Tfloats)

    for Cl in [Cl_low, Cl_high]:
        for s in Cl:
            for P in [1., 100.]:
                for m in ['pitzerPHREEQCnoGases', 'SUPCRTnoGases']:

                    ESS = ESSaver(
                      pHfloats,
                      Tfloats,
                      s,
                      'spec_T_273-473_pH_7-12',
                      P=P,
                      model=m,
                      fixGases=False)

                    ESS.save_spec_params_1salt(s)

                for m in ['pitzerPHREEQC', 'SUPCRT']:
                    # for investigating saturation index

                    ESS = ESSaver(
                      pHfloats,
                      Tfloats,
                      s,
                      'spec_T_273-473_pH_7-12',
                      P=P,
                      model=m,
                      fixGases=True)

                    ESS.save_spec_params_1salt(s)




""" Speciation with a more narrow pH window """

def narrow_pH():
    Cl_low = np.logspace(math.log10(0.05),math.log10(0.2), num=21)
    Cl_low = np.round(Cl_low, decimals=5)
    Cl_high = np.logspace(math.log10(0.2),math.log10(0.4), num=11)
    Cl_high = np.round(Cl_high, decimals=5)

    pHfloats = np.linspace(8.,9.,num=6)
    pHfloats = np.round(pHfloats, decimals=1)
    Tfloats = np.linspace(273.15, 473.15, num=41)

    for Cl in [Cl_low, Cl_high]:
        for s in Cl:
            for P in [1., 100.]:

                ESS = ESSaver(
                  pHfloats,
                  Tfloats,
                  s,
                  'spec_T_273-473_pH_8-9',
                  P=P,
                  model='pitzerPHREEQCnoGases',
                  fixGases=False)

                ESS.save_spec_params_1salt(s)


# broad_pH()
# narrow_pH()

""" For a very fine grid. This was not used for results in H+2024"""

# def dense():
#     Cl = np.linspace(0.05,0.4, num=71)
#     Cl = np.round(Cl,decimals=5)
#
#     pHfloats = np.linspace(7.,12.,num=51)
#     pHfloats = np.round(pHfloats, decimals=1)
#     Tfloats = np.linspace(273.15, 403.15, num=53)
#
#     for s in Cl:
#         for P in [1., 100.]:
#
#             ESS = ESSaver(
#               pHfloats,
#               Tfloats,
#               s,
#               'spec_T_273-403_pH7-12_highdensity',
#               P=P,
#               model='pitzerPHREEQC',
#               fixGases=True)
#
#             ESS.save_spec_params_1salt(s)
