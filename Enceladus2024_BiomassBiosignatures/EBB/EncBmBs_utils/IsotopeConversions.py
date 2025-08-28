# Some handy isotope conversions

class IsotopeConversions:

    def RtodC(R):
        return ((R/0.0112372)-1)*1e3

    def dCtoR(dC):
        return ((dC/1e3)+1)*0.0112372

    def atoe(a):
        return a -1.

    def etoa(e):
        return e + 1.

    def FtoR(F):
        return F / (1-F)

    def RtoF(R):
        return R / (1+R)
