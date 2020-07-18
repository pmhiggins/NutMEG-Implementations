from methanogen_implementer import efficiencies
import matplotlib.pyplot as plt



def concs_vs_t(plot=True, Trange=range(275,375),CH4=3e-8):
    """ plot the concentration of CO2, H2, and CH4 with temperature"""
    CO2lst=[]
    H2lst=[]
    for T in Trange:
        peff = efficiencies.get_eff('averageMethanogen', Temp=T)
        CO2lst.append(peff.params['mol_CO2'])
        H2lst.append(peff.params['mol_H2'])
    if plot:
        plt.plot(Trange, CO2lst, label=r'[$\mathregular{CO}_2$]')
        plt.plot(Trange, H2lst, label=r'[$\mathregular{H}_2$]')
        plt.plot(Trange, [CH4]*len(Trange), label=r'[$\mathregular{CH}_4$]')
        plt.ylabel('Concentration [M]')
        plt.xlabel('Temperature [K]')
        # plt.yscale('log')
        plt.legend()
        plt.show()
    return CO2lst, H2lst

# concs_vs_t()
