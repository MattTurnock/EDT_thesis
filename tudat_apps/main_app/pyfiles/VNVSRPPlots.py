import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
from tudatApplications.EDT_thesis.LitStudy import LitStudy3
from astropy import units as u

matplotlib.rcParams.update({'font.size': 20})

########################################################################################################################
#####################  #############################################################
########################################################################################################################

VnVSRPNominalSubDir = "SRPVnV"
VnVSRPZeroAreaSubDir = "SRPZeroAreaVnV"
VnVSaveFolder = os.path.join("pyplots", "VnV")
utils.checkFolderExist(VnVSaveFolder)
showing=True
defaultFigsize = [16,9]

dataToLoad = ["propData", "depVarData"]
plotLegend = ["Tudat Simulation", "Reference Calculations"]

################## Do zero Area data #############################

allSimDataSRPZeroArea = utils.getAllSimDataFromFolder(VnVSRPZeroAreaSubDir, todoList=dataToLoad)
propDataArrayZeroArea = allSimDataSRPZeroArea[5]
depVarDataArrayZeroArea = allSimDataSRPZeroArea[2]

accelerationArrayZeroArea = depVarDataArrayZeroArea[:, 8:11]
coordsArrayZeroArea = propDataArrayZeroArea[:, 1:4]

magnitudeAccelerationArrayZeroArea = np.linalg.norm(accelerationArrayZeroArea, axis=1)
timesZeroArea = propDataArrayZeroArea[:, 0]
radiiZeroArea = np.linalg.norm(coordsArrayZeroArea, axis=1)

plt.figure(1, figsize=defaultFigsize)
# plt.plot(radiiZeroArea, magnitudeAccelerationArrayZeroArea)


####################### Do nominal data #####################################

allSimDataSRPNominal = utils.getAllSimDataFromFolder(VnVSRPNominalSubDir, todoList=dataToLoad)
propDataArrayNominal = allSimDataSRPNominal[5]
depVarDataArrayNominal = allSimDataSRPNominal[2]

accelerationArrayNominal = depVarDataArrayNominal[:, 8:11]
coordsArrayNominal = propDataArrayNominal[:, 1:4]

magnitudeAccelerationArrayNominal = np.linalg.norm(accelerationArrayNominal, axis=1)
timesNominal = propDataArrayNominal[:, 0]
radiiNominal = np.linalg.norm(coordsArrayNominal, axis=1) / utils.AU

plt.figure(1, figsize=defaultFigsize)
# plt.plot(radiiNominal, magnitudeAccelerationArrayNominal)


indices = utils.getLogarithmicResizeIndices(radiiNominal, 30)
plt.scatter(radiiNominal[indices], magnitudeAccelerationArrayNominal[indices], marker='x')





################################## Do reference plotting from "hand" calculations ############################

# Values taken from output of cpp files
vehicleRadiationCoefficient = 0.5
vehicleArea = 2.4
vehicleMass = 100

radiiReference = radiiNominal
magnitudeAccelerationArrayReference = np.zeros(np.shape(magnitudeAccelerationArrayNominal))
for i in range(len(magnitudeAccelerationArrayReference)):
    # print(radiiReference[i])
    magnitudeAccelerationArrayReference[i] = abs(LitStudy3.get_fbarrad(vehicleRadiationCoefficient,
                                                                   radiiReference[i]*utils.AU,
                                                                   vehicleArea,
                                                                   vehicleMass,
                                                                   W0=utils.W0,
                                                                   r0=utils.AU,
                                                                   c=utils.c,
                                                                   units=None))

differenceArray = magnitudeAccelerationArrayNominal - magnitudeAccelerationArrayReference
print("Maximum difference in SRP acceleration: %s m/s^2" %abs(np.max(differenceArray)))

plt.figure(1)
plt.plot(radiiReference, magnitudeAccelerationArrayReference, c='C1')
# plt.scatter(radiiReference, magnitudeAccelerationArrayReference, c='C1', marker='x')

## Add graph markings etc ##

plt.xlabel("Radius [AU]")
plt.ylabel("SRP Acceleration Magnitude [m/s$^2$]")
plt.grid(which="both")
plt.yscale("log")
plt.xscale("log")
plt.xlim([1E-1, 1E3])
plt.ylim([1E-14, 1E-6])
plt.legend(plotLegend)
plt.savefig(os.path.join(VnVSaveFolder, "SRP-Nominal.png"))



if showing: plt.show()

