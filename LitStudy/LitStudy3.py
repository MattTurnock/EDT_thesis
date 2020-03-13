# In this one we graph the planets effect, assuming closest approach

import numpy as np
from json_to_dict import constants
from astropy import units as u
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
from scipy import special
from numpy import pi

from LitStudy1 import getB_2, B0_min, B0_max, phi0_min, phi0_max, R0, fig1savespace
u.imperial.enable()

def list_to_dimensionless(old_list):
    newlist = []
    for i in old_list:
        newlist.append(i.value)

    return newlist

def get_a(mu, r, units=None):
    a = mu/r**2

    if units is not None:
        a = a.to(units)

    return a


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

def get_frmax_J2(mu, J2, R, r, units=None):

    frmax = 3*mu*J2*(R**2/r**4)

    if units is not None:
        frmax = frmax.to(units)

    return frmax

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


ranges = []
ranges_base = np.logspace(-1,3, 1000)

for i in ranges_base:
    ranges.append(i*u.AU)



# Data nested list contains grav constant, radius from sun, empty list for appending accels per radius
data = [[ constants["muSun"]    , 0*u.AU               , [], "Sun",    [], [], [] ],
        [ constants["muMercury"], constants["rMercury"], [], "Mercury",[], [], [] ],
        [ constants["muVenus"], constants["rVenus"]    , [], "Venus",  [], [], [] ],
        [ constants["muMars"], constants["rMars"]      , [], "Mars",   [], [], [] ],
        [ constants["muEarth"], constants["rEarth"]    , [], "Earth",  [], [], [] ],
        [ constants["muJupiter"], constants["rJupiter"], [], "Jupiter",[], [], [] ],
        [ constants["muSaturn"], constants["rSaturn"]  , [], "Saturn", [], [], [] ],
        [ constants["muUranus"], constants["rUranus"]  , [], "Uranus", [], [], [] ],
        [ constants["muNeptune"], constants["rNeptune"], [], "Neptune",[], [], [] ]]

muSun = constants["muSun"]

fr_ratios = []
fbarrad_ratios = []
EDT_ratios = []
EDT_raw = []

for r in ranges:
    a_sun = get_a(muSun, r, units=u.m / u.s ** 2)
    #3rd body changes:
    for dat in data:
        mu = dat[0]
        r_body = dat[1]
        accelList_ca = dat[2]
        name = dat[3]
        ratioList_ca = dat[4]
        accelList_av = dat[5]
        ratioList_av = dat[6]

        d_body_ca = abs(r - r_body) # calculate closest approach distance between s/c and body
        d_body_av = get_dbar(r, r_body, units=u.AU)

        a_ca = get_a(mu, d_body_ca, units=u.m/u.s**2)
        a_av = get_a(mu, d_body_av, units=u.m/u.s**2)


        a_ca_ratio = a_ca/a_sun
        a_av_ratio = a_av / a_sun

        accelList_ca.append(a_ca)
        accelList_av.append(a_av)
        ratioList_ca.append(a_ca_ratio)
        ratioList_av.append(a_av_ratio)

    #SH Perts:
    frmax = get_frmax_J2(muSun, constants["J2_sun"], constants["RSun"], r, units=u.m / u.s ** 2)

    fr_ratios.append(abs(frmax / a_sun))

    # Radiation Perturbations:
    CR = 1.5
    A = (2 * (27 * 14 * u.imperial.inch ** 2) + 3 * u.mm * 4 * u.km).to(u.m ** 2)  # First term for sat, second term for tether
    M = (127 * u.imperial.lb).to(u.kg)
    fbarrad = get_fbarrad(CR, r, A, M)
    fbarrad_ratios.append( abs(fbarrad / a_sun))

    # EDT accels:
    a_EDT_0 = 0.325 * u.mm/u.s**2
    B_EDT_0 = 0.5*(25+65)*u.uT
    B = getB_2( 0.5*(B0_max + B0_min),
                0.5*(phi0_max + phi0_min),
                R0,
                r,
                u.nT )
    a_EDT = a_EDT_0 * (B / B_EDT_0)
    a_EDT = a_EDT.to(u.m/u.s**2)
    EDT_ratios.append(a_EDT / a_sun)
    EDT_raw.append(a_EDT)


# plt.figure()
# labels = []
# for dat in data:
#     plt.plot(list_to_dimensionless(ranges), list_to_dimensionless(dat[4]))
#     labels.append(dat[3])
# plt.legend(labels)
# ax=plt.gca()
# ax.set_xscale('log')
# ax.set_yscale('log')
# plt.ylabel("ratio of accel to solar accel [-]")
# plt.xlabel("distance from sun [AU]")
# plt.grid(which='both')
#
# ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
# ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
# plt.show()








planetsToInclude = ["Sun", "Earth", "Jupiter", "Saturn", "Uranus", "Neptune"]

# Plot average acceleration vs distance from sun
plt.figure(figsize=(8, 6))
labels = []
for dat in data:
    name = dat[3]
    if name in planetsToInclude:
        if name == "Sun":
            linestyle = "-"
        else:
            linestyle = "--"
        plt.plot(list_to_dimensionless(ranges), list_to_dimensionless(dat[6]), linestyle=linestyle)
        labels.append(dat[3])
labels.append("SH")
labels.append("Radiation")
labels.append("EDT")
plt.plot(list_to_dimensionless(ranges), list_to_dimensionless(fr_ratios))
plt.plot(list_to_dimensionless(ranges), list_to_dimensionless(fbarrad_ratios))
plt.plot(list_to_dimensionless(ranges), list_to_dimensionless(EDT_ratios))
plt.legend(labels)
ax=plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.ylabel("Normalised perturbation acceleration [-]")
plt.xlabel("Distance from sun [AU]")
plt.ylim([1e-10, 1e1])
plt.grid(which='both')

ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
plt.savefig(fig1savespace %"OD_perturbations.pdf" )
# plt.show()



plt.figure()
plt.plot(list_to_dimensionless(ranges), list_to_dimensionless(EDT_raw))
ax=plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.ylabel("accel [m/s**2]")
plt.xlabel("distance from sun [AU]")
plt.grid(which='both')

ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
plt.show()


