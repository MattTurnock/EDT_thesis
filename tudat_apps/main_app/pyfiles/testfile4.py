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
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import sys



# utils.porkchopPlot("/home/matt/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/SimulationOutput",
#              "porkchopEarthMars", EarthMarsCorrection=True, ylabel="Travel Time [days]")
# utils.porkchopPlot("/home/matt/LinkToEDT_thesis/tudat_apps/main_app/SimulationOutput/GAJupiter/GACalculatorNominal_Jupiter_2020.00-2030.00",
#              "porkchop_GA_EJ")
#
# plt.show()

# minDV = np.amin(porkchopExample)
# minDVIndex = np.where(porkchopExample == minDV)
# minDVYear = porkchopExample_x_data[ minDVIndex[1][0] ]
# minDVTOF = porkchopExample_y_data[minDVIndex[0][0] ]
# print(minDV)
# print(minDVIndex)
# print(minDVYear)
# print(minDVTOF)


# Very slow for many datapoints.  Fastest for many costs, most readable
def getParetoArray(costs, returnParetoArray=True):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
    """

    is_efficient = np.ones(costs.shape[0], dtype = bool)
    for i, c in enumerate(costs):
        is_efficient[i] = np.all(np.any(costs[:i]>c, axis=1)) and np.all(np.any(costs[i+1:]>c, axis=1))

    if returnParetoArray:
        paretoArray = costs[is_efficient]
        return paretoArray
    else:
        return is_efficient

def getParetoLists(XInput, YInput, sortOutput=True):
    if len(XInput) != len(YInput):
        print("ERROR: Input lengths of lists for pareto are not equal")
        sys.exit()

    inputArray = np.zeros((len(XInput), 2))
    inputArray[:,0] = XInput
    inputArray[:,1] = YInput

    paretoArray = getParetoArray(inputArray, returnParetoArray=True)
    if sortOutput:
        paretoArray = paretoArray[paretoArray[:,0].argsort()]
    paretoXs = paretoArray[:,0]
    paretoYs = paretoArray[:,1]

    return paretoXs, paretoYs


xs = [1,0.9,2,2]
ys = [1,2,0.9,2,]
# dvTofs = np.zeros((len(xs), 2))
# dvTofs[:,0] = xs
# dvTofs[:,1] = ys

newXs, newYs = getParetoLists(xs, ys)
print(newXs, newYs)

# dvTofs = np.array([[1,1], [0.9,2], [2,0.9], [2,2]])
# dvTofsPar = is_pareto_efficient_dumb(dvTofs, returnParetoArray=True)
# # newdvTofs = dvTofs[dvTofsPar]
# print(dvTofsPar)
# print(dvTofs)
# print(newdvTofs)

plt.figure()
plt.scatter(xs, ys)
plt.scatter(newXs, newYs)
plt.plot(newXs, newYs)
plt.show()