from json_to_dict import constants
from astropy import units as u
import numpy as np
from numpy import pi, sin, cos, tan
from scipy import special
u.imperial.enable()

def get_frmax_J2(mu, J2, R, r, units=None):

    frmax = 3*mu*J2*(R**2/r**4)

    if units is not None:
        frmax = frmax.to(units)

    return frmax

g = 9.81*u.m/u.s**2

r1_earth = constants["RE"] + 250*u.km


frmax1_J2_earth = get_frmax_J2(constants["muEarth"],
             constants["J2"],
             constants["RE"],
             r1_earth, units=u.cm/u.s**2)
print("J2 acceleration at 250km Earth orbit = %s (%s g)" %(frmax1_J2_earth, (frmax1_J2_earth/g).decompose()))

r_AU = 1

frmax_1au_J2_sun = get_frmax_J2(constants["muSun"],
                                constants["J2_sun"],
                                constants["RSun"],
                               r_AU*constants["AU"],
                                units=u.m/u.s**2)
print("J2 acceleration at %sAU solar orbit = %s (%s g)" %(r_AU, frmax_1au_J2_sun, (frmax_1au_J2_sun/g).decompose()))

######################################################################################################################

def get_fd__fs(md, ms, ri, rd, units=False):
    """
    gets ratio of perturbing accel (fd) to solar accel (fs)
    :param md: perturbing mass
    :param ms: solar mass
    :param ri: distance to main body (sun)
    :param rd: distance to perturbing body
    :return:
    """

    fdfs = 2*(md/ms)*(ri/rd)**3
    if units is True:
        fdfs = fdfs.decompose()

    return fdfs

def get_a(mu, r, units=None):
    a = mu/r**2

    if units is not None:
        a = a.to(units)

    return a

print(get_fd__fs(constants["muJupiter"],
                 constants["muEarth"],
                 42164*u.km,
                 constants["rJupiter"] - 1*u.AU,
                 units=True))

print("\n========================\n")

aE_1_d = 1.3*u.AU
aE_1 = get_a(constants["muEarth"], aE_1_d, units=u.m/u.s**2)

aE_2 = get_a(constants["muEarth"], 1*u.AU, units=u.m/u.s**2)

aS_1_d = 1*u.AU
aS_1 = get_a(constants["muSun"], aS_1_d, units=u.m/u.s**2)

aJ_1_d = 23.5*u.AU
aJ_1 = get_a(constants["muJupiter"], aJ_1_d, units=u.m/u.s**2)

aSat_1_d = 13*u.AU
aSat_1 = get_a(constants["muSaturn"], aSat_1_d, units=u.m/u.s**2)

aU_1_d = 5*u.AU
aU_1 = get_a(constants["muUranus"], aU_1_d, units=u.m/u.s**2)

aN_1_d = 5.5*u.AU
aN_1 = get_a(constants["muNeptune"], aN_1_d, units=u.m/u.s**2)



print("Earth at %s: %s" %(aE_1_d, aE_1 ))
# print(aE_2)
# print(aS_1)
print("Jupiter at %s: %s" %(aJ_1_d, aJ_1 ))
# print(constants["rSaturn"])
print("Saturn at %s: %s" %(aSat_1_d, aSat_1 ))
print("Uranus at %s: %s" %(aU_1_d, aU_1 ))
print("Neptune at %s: %s" %(aN_1_d, aN_1 ))


print("\n========================\n")

def get_fbarrad(CR, r, A, M, W0=constants["SC"], r0=1*u.AU, c=constants["c"], units=u.m/u.s**2):
    """
    get radiation accel (see wakker)
    :param CR: Radiation coefficient
    :param r: Distance from sun
    :param A: Sat area
    :param M: Sat mass
    :param W0: Solar flux at r0
    :param r0: Dikstance of known W0
    :param c: Speed of light
    :param units:
    :return:
    """
    W = W0 * (r0/r)**2

    fbar = -CR*(W/c)*(A/M)

    if units is not None:
        fbar = fbar.to(units)

    return fbar

CR = 1.5
A = (2*(27*14*u.imperial.inch**2) + 3*u.mm*4*u.km).to(u.m**2)    # First term for sat, second term for tether
M = (127*u.imperial.lb).to(u.kg)

r=12*u.AU

fbar = get_fbarrad(CR, r, A, M)
# fbar = get_fbarrad(1.9, r, 12*u.m**2, 1*u.kg, units=u.mm/u.s**2)
print("fbar is")
print(fbar)


#####################################################################################################################
# Plotting planetary distances

def get_dbar(r1, r2, units=None):
    """
    Find average distance between 2 orbit radii r1 and r2
    :param r1: radius 1
    :param r2: radius 2
    :param units: units to output
    :return:
    """

    if units is not None:
        r1 = (r1.to(u.AU)).value
        r2 = (r2.to(u.AU)).value

    dbar = (2/pi)*(r1 + r2) * special.ellipe((2*np.sqrt(r1*r2))/(r1+r2))

    if units is not None:
        dbar = (dbar*u.AU).to(units)

    return dbar

print(get_dbar(1*u.AU, constants["rJupiter"], units=u.AU))

