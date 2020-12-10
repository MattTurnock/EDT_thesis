import numpy as np

pi = np.pi
cos=np.cos
sin=np.sin

phi0_default = 55*(pi/180)

def getMagfield(R_, B0_=12E-9, phi0_=phi0_default, R0_=1, printing=True):
    Bx = B0_ * cos(phi0_) * pow((R0_/R_), 2)
    By = B0_ * sin(phi0_) * pow((R0_/R_), 1)
    Bz = 0
    Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)

    result = [Bx, By, Bz, Bmag]
    if printing == True:
        print("Magfield for alt %s AU: %s" %(R_, result))

    return result

getMagfield(1)
getMagfield(0.9947)
getMagfield(1.028)
getMagfield(2.7)