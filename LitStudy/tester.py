import numpy as np
from json_to_dict import constants
from astropy import units as u
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
from scipy import special
from numpy import pi

def get_a(mu, r, units=None):
    a = mu/r**2

    if units is not None:
        a = a.to(units)

    return a

a_1au = get_a(constants["muSun"], 1*u.AU, units=u.m/u.s**2)

print("solar accel at 1AU: ", a_1au)