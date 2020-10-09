# import numpy as np
# from matplotlib import pyplot as plt
#
# x = np.arange(1, 10)
# y = x.reshape(-1, 1)
# h = x * y
#
# cs = plt.contourf(h, levels=[10, 30, 50],
#                   colors=['#808080', '#A0A0A0', '#C0C0C0'], extend='both')
# cs.cmap.set_over('red')
# cs.cmap.set_under('blue')
# cs.changed()
# print(h)
# plt.show()

import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import os
import scipy.interpolate as si

TOFS = [1,2,3]
LaunchYears = [2020, 2021, 2022]
DeltaVs = [10, 12, 13]



# delta = 0.025
# x = np.arange(-3.0, 3.0, delta)
# y = np.arange(-2.0, 2.0, delta)
# X, Y = np.meshgrid(x, y)
# Z1 = np.exp(-X**2 - Y**2)
# Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
# Z = (Z1 - Z2) * 2

exmaplePath = "/home/matt/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/SimulationOutput"
porkchopEarthMars_x_data = np.genfromtxt(os.path.join(exmaplePath, "porkchopEarthMars_x_data.dat"))
porkchopEarthMars_y_data = np.genfromtxt(os.path.join(exmaplePath, "porkchopEarthMars_y_data.dat"))
porkchopEarthMars = np.genfromtxt(os.path.join(exmaplePath, "porkchopEarthMars.dat"))
porkchopEarthMars = np.transpose(porkchopEarthMars)

porkchopExample_x_data = np.genfromtxt(os.path.join(exmaplePath, "porkchopEXAMPLE_x_data.dat"))
porkchopExample_y_data = np.genfromtxt(os.path.join(exmaplePath, "porkchopEXAMPLE_y_data.dat"))
porkchopExample = np.genfromtxt(os.path.join(exmaplePath, "porkchopEXAMPLE.dat"))
porkchopExample = np.transpose(porkchopExample)


plt.figure(1)
plt.contour((porkchopEarthMars_x_data-2451545)/365,porkchopEarthMars_y_data,porkchopEarthMars, 50,
            linewidths=0.1, colors="black")
plt.contourf((porkchopEarthMars_x_data-2451545)/365,porkchopEarthMars_y_data,porkchopEarthMars, 50, cmap=plt.cm.viridis)
plt.xlabel('Departure date [years since J2000]')
plt.ylabel('Travel time [days]')
plt.title('Porkchop plot, two-shot impulsive Earth-Mars transfer [m/s]')
plt.colorbar()

plt.figure(2)
plt.contour(porkchopExample_x_data, porkchopExample_y_data, porkchopExample,
            linewidths=0.5, colors="black")
plt.contourf(porkchopExample_x_data, porkchopExample_y_data, porkchopExample)
plt.xlabel('Departure date [Years since J2000]')
plt.ylabel('Travel time [Years]')
plt.title('Porkchop plot, impulsive Earth-Jupiter transfer [km/s]')
plt.colorbar()

Xflat, Yflat, Zflat = porkchopExample_x_data.flatten(), porkchopExample_y_data.flatten(), porkchopExample.flatten()
# def fmt(x, y):
#     # get closest point with known data
#     dist = np.linalg.norm(np.vstack([Xflat - x, Yflat - y]), axis=0)
#     idx = np.argmin(dist)
#     z = Zflat[idx]
#     return 'x={x:.5f}  y={y:.5f}  z={z:.5f}'.format(x=x, y=y, z=z)

def fmt(x, y):
    z = np.take(si.interp2d(porkchopExample_x_data, porkchopExample_y_data, porkchopExample)(x, y), 0)
    return 'x={x:.5f}  y={y:.5f}  z={z:.5f}'.format(x=x, y=y, z=z)

plt.gca().format_coord = fmt

minDV = np.amin(porkchopExample)
minDVIndex = np.where(porkchopExample == minDV)
minDVYear = porkchopExample_x_data[ minDVIndex[1][0] ]
minDVTOF = porkchopExample_y_data[minDVIndex[0][0] ]
print(minDV)
print(minDVIndex)
print(minDVYear)
print(minDVTOF)

plt.show()