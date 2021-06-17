import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
from astropy import units as u
import itertools

showing = True
onlyDefineFunctions = False

# def getB_p10(B0, R0, R, units=None):
#     """
#     Get magnetic field strength above 10AU
#     :param B0:
#     :param R0:
#     :param R:
#     :param units:
#     :return:
#     """
#     B = B0 * (R0/R)
#
#     if units is not None:
#         B = B.to(units)
#
#     return B
#
# def getB_b10(B0, R0, R, units=None):
#     """
#     Magnetic field strength below 10AU
#     :param B0:
#     :param R0:
#     :param R:
#     :param units:
#     :return:
#     """
#     B = B0 * (R0/R)**2
#
#     if units is not None:
#         B = B.to(units)
#
#     return B

def getB_2(B0, phi0, R0, R, units=None):
    """
    Magnetic field strength
    :param B0:
    :param R0:
    :param R:
    :param units:
    :return:
    """
    BR_0 = B0*np.cos(phi0)
    Bphi_0 = B0*np.sin(phi0)

    BR = BR_0 * (R0/R)**2
    Bphi = Bphi_0 * (R0/R)

    B = np.sqrt(BR**2 + Bphi**2)

    if units is not None:
        B = B.to(units)

    return B

def list_to_dimensionless(old_list):
    newlist = []
    for i in old_list:
        newlist.append(i.value)

    return newlist

if onlyDefineFunctions is False:

    B0_min = 4 * u.nT
    B0_max = 12 * u.nT
    B0s = [B0_min, B0_max]
    phi0_min = 35*u.deg
    phi0_max = 55*u.deg
    phi0s = [phi0_min, phi0_max]

    permutations = itertools.product(B0s, phi0s)

    R0 = 1 * u.AU

    radii = []
    radii_tmp = np.linspace(0.1, 80, 100)

    for R_tmp in radii_tmp:
        radii.append(R_tmp*u.AU)

    B11 = []
    B12 = []
    B21 = []
    B22 = []

    for R in radii:
        # if R <= 10*u.AU:
        #     Bmax10 = getB_b10(B0_max, R0, R, u.nT)
        #     Bmin10 = getB_b10(B0_min, R0, R, u.nT)
        #     Bmax.append(Bmax10)
        #     Bmin.append(Bmin10)
        #     R10 = 10*u.AU
        #
        # elif R > 10*u.AU:
        #     Bmax.append(getB_p10(Bmax10, R10, R, u.nT))
        #     Bmin.append(getB_p10(Bmin10, R10, R, u.nT))
        # else:
        #     print("failed")
        #
        # print(Bmin10)
        # print(Bmax10)

        B11.append(getB_2(B0_min, phi0_min, R0, R, u.nT))
        B12.append(getB_2(B0_min, phi0_max, R0, R, u.nT))
        B21.append(getB_2(B0_max, phi0_min, R0, R, u.nT))
        B22.append(getB_2(B0_max, phi0_max, R0, R, u.nT))


    # Source surface is a

    # fig, ax = plt.subplots()
    # # ax.plot(list_to_dimensionless(radii), list_to_dimensionless(Bmax))
    # # ax.plot(list_to_dimensionless(radii), list_to_dimensionless(Bmin))
    # # ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    # # ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
    # # plt.grid(which='both')
    # # plt.yscale("log")
    # # plt.xscale("log")
    # # plt.xlabel("Distance from Sun, R [AU]")
    # # plt.ylabel("Magnetic Field Magnitude, B [nT]")
    # # plt.legend(["Solar Maximum", "Solar Minimum"])
    # # plt.show()


    fig1savespace = "E:\\Documents\\MEGA\\NEW\\delft\\schoolwork\\Master\\Thesis\\Literature Study\\images\\%s"
    plt.figure()
    plt.plot(list_to_dimensionless(radii), list_to_dimensionless(B11))
    plt.plot(list_to_dimensionless(radii), list_to_dimensionless(B12))
    plt.plot(list_to_dimensionless(radii), list_to_dimensionless(B21))
    plt.plot(list_to_dimensionless(radii), list_to_dimensionless(B22))
    plt.axvline(x=0.39, linestyle=':', color='black')
    plt.axvline(x=5.20, linestyle='--', color='black')
    plt.axvline(x=80, linestyle='-', color='black')
    ax=plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.xlabel("Distance from Sun, R [AU]")
    plt.ylabel("Magnetic Field Strength, B [nT]")
    tmp_string = "$B_0$ = %s, $\phi_0$ = %s"
    plt.legend([tmp_string %(B0_min, phi0_min),
                tmp_string %(B0_min, phi0_max),
                tmp_string %(B0_max, phi0_min),
                tmp_string %(B0_max, phi0_max),
                "Mercury Orbit Radius",
                "Jupiter Orbit Radius",
                "(Aprroximate) Termination Shock"])
    plt.grid(which='both')

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
    # ax.xaxis.get_minor_formatter().set_scientific(False)
    # ax.xaxis.get_minor_formatter().set_useOffset(False)
    plt.savefig(fig1savespace %"HMF_char.pdf" )
    if showing: plt.show()