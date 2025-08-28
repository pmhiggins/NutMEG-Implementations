import math

class CarbonateFracFactorDeines1974:

    def __init__(self, T=298):

        self.T = T
        self.a_CO3 = self.CO3_wrt_CO2(T)
        self.a_HCO3 = self.HCO3_wrt_CO2(T)


    def CO3_wrt_CO2(self, T):

        KH2CO2exp = (-0.91 + (0.0063e6/(T*T)))
        CO3exp = (-3.4 + (0.87e6/(T*T)))

        return math.exp(0.001*(CO3exp-KH2CO2exp))

    def HCO3_wrt_CO2(self, T):

        KH2CO2exp = (-0.91 + (0.0063e6/(T*T)))
        HCO3exp = (-4.54 + (1.009e6/(T*T)))

        return math.exp(0.001*(HCO3exp-KH2CO2exp))

    def update(self, T):
        self.__init__(T)
