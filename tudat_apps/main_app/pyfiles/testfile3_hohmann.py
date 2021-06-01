import numpy as np
from astropy import units as u

pi = np.pi





def getTOF(a, mu, unitsOut=u.s):
    T = 2*pi*np.sqrt( a**3 / mu )

    T = T.to(unitsOut)

    return 0.5*T

def getV(mu, r, a, unitsOut=u.m/u.s):
    V = np.sqrt( mu* (2/r - 1/a) )
    V = V.to(unitsOut)
    return V

def getDVHohmann(a1, a2, mu, unitsOut=u.m/u.s):
    hohmanna = 0.5*(a1 + a2)
    V1circ = getV(mu, a1, a1, unitsOut=unitsOut)
    V2circ = getV(mu, a2, a2, unitsOut=unitsOut)
    V1 = getV(mu, a1, hohmanna, unitsOut=unitsOut)
    V2 = getV(mu, a2, hohmanna, unitsOut=unitsOut)

    DV1 = abs(V1circ - V1)
    DV2 = abs(V2circ - V2)

    return [DV1, DV2]


def getOrbitalParameters(V, r, mu, rPe=True):

    a = (2/r - V**2/mu)**(-1)

    if rPe:
        e = 1 - r/a
    else:
        e = r/a - 1

    Ap = a*(1+e)
    Pe = a*(1-e)


    return (a, e, Ap, Pe)




R_E = 1*u.AU
muE = 3.986E14 * u.m**3*u.s**-2

R_J = 5.2*u.AU
mu_J = 1.267E17 * u.m**3*u.s**-2

R_Sat = 9.54 * u.AU
mu_Sat = 3.793E16 * u.m**3*u.s**-2

R_M = 1.524 * u.AU
mu_M = 0.107 * muE

mu_S = 1.327E20 * u.m**3*u.s**-2

a_EJ = (R_E + R_J)*0.5
a_ES = (R_E + R_Sat) * 0.5
a_EM = (R_E + R_M) * 0.5
TOFJ = getTOF(a_EJ, mu_S, unitsOut=u.year)
TOFS = getTOF(a_ES, mu_S, unitsOut=u.year)
TOFM = getTOF(a_EM, mu_S, unitsOut=u.year)
DV1J = getDVHohmann(R_E, R_J, mu_S, unitsOut=u.km/u.s)
DV1S = getDVHohmann(R_E, R_Sat, mu_S, unitsOut=u.km / u.s)
DV1M = getDVHohmann(R_E, R_M, mu_S, unitsOut=u.km/u.s)

print("Hohmann TOF to Jupiter: ", TOFJ)
print("Hohmann TOF to Saturn: ", TOFS)
print("Hohmann TOF to Mars: ", TOFM)
print("Hohmann DV to Jupiter: ", DV1J[0])
print("Hohmann DV to Saturn: ", DV1S[0])
print("Hohmann DV to Mars: ", DV1M[0])



DVMax = 10*u.km/u.s

V_E = getV(mu_S, R_E, R_E, unitsOut=u.km/u.s)
Vmax = V_E + DVMax
Vmin = V_E - DVMax

aMax, eMax, ApMax, PeMax = getOrbitalParameters(Vmax, R_E, mu_S, rPe=True)
aMin, eMin, ApMin, PeMin = getOrbitalParameters(Vmin, R_E, mu_S, rPe=False)


print("\nWith +DVMax = %s:" %DVMax)
print("Earth V: %s\n"
      "Vmax: %s\n"
      "aMax: %s\n"
      "eMax: %s\n"
      "ApMax: %s\n"
      "PeMax: %s" %(V_E, Vmax, aMax, eMax, ApMax, PeMax))

print("\nWith -DVMax = %s:" %DVMax)
print("Earth V: %s\n"
      "Vmin: %s\n"
      "aMin: %s\n"
      "eMin: %s\n"
      "ApMin: %s\n"
      "PeMin: %s" %(V_E, Vmin, aMin, eMin, ApMin, PeMin))