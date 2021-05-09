# File to create compression module for saved data (ie if like 6gb should help shrink for doing plots etc)

import numpy as np
import os
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils

matplotlib.rcParams.update({'font.size': 20})

year = 365*24*60*60
AU = 1.496E11

def compressSimData(allSimData, interpolatingTimestep = 1E-5, savingSubDirectory=None, filenamePrefix="TRAJ"):

    interpolatedDataArraysList = []

    baseSaveFilename = filenamePrefix + "-%s-Compressed.dat"

    for i in range(len(allSimData)):
        ithDataArray = allSimData[i]

        ithDataTimes = ithDataArray[:,0]
        newtimes = np.arange(ithDataTimes[0], ithDataTimes[-1], interpolatingTimestep)
        ithDataInterpolated = utils.interpolateArrays(ithDataArray, newtimes)

        interpolatedDataArraysList.append(ithDataInterpolated)

        if savingSubDirectory is not None:
            if i == 0:
                dataType = "bodyData"
            elif i == 1:
                dataType = "currentData"
            elif i == 2:
                dataType = "depVarData"
            elif i == 3:
                dataType = "ionoData"
            elif i == 4:
                dataType = "magData"
            elif i == 5:
                dataType = "propData"
            elif i == 6:
                dataType = "thrustData"
            elif i == 7:
                dataType = "currentVNVData"
            elif i == 8:
                dataType = "configInfo"
            else:
                print("Too many indices for compression, ignoring.")

            if dataType is not "configInfo":
                saveFilename =  baseSaveFilename %dataType
                savePath = os.path.join(utils.simulation_output_dir, savingSubDirectory, saveFilename)
                np.savetxt(savePath, ithDataInterpolated, delimiter=", ")

    interpolatedDataArraysTuple = tuple(interpolatedDataArraysList)

    return interpolatedDataArraysTuple

##################### GA Plots ########################
DEFAULTSIZE = None
SaveFolder = "pyplots/finalSimsTemp"
interpTimestep = 1E5

# fullLoadTodoList = ["propData"]
fullLoadTodoList = ["bodyData", "currentData", "currentVNVData", "depVarData", "ionoData", "magData", "propData", "thrustData"]

dataSubDir_InO = "InO-Stage2"
allSimData_InO = utils.getAllSimDataFromFolder(dataSubDir_InO, todoList=fullLoadTodoList, useCompressed=False)




propData_InO = allSimData_InO[5]
# propDataTimes_InO = propData_InO[:,0]
# newtimes_InO = np.arange(propDataTimes_InO[0], propDataTimes_InO[-1], interpTimestep)
# propData_InO_Interpolated = utils.interpolateArrays(propData_InO, newtimes_InO)


# np.savetxt(os.path.join(utils.simulation_output_dir, dataSubDir_InO, "InO-Stage2-proData-Compressed.dat"), propData_InO_Interpolated)


allDataInterpolated = compressSimData(allSimData_InO, interpTimestep, savingSubDirectory=dataSubDir_InO, filenamePrefix="InO-Stage2")
propData_InO_Interpolated = allDataInterpolated[5]



utils.plotTrajectoryData(propData_InO, sameScale=True, planetsToPlot=["Earth"], plotSun=True, fignumber=1, plotOnlyTrajectory=False, trajectoryLabel="SSO Spacecraft")

utils.plotTrajectoryData(propData_InO_Interpolated, sameScale=True, planetsToPlot=["Earth"], plotSun=True, fignumber=2, plotOnlyTrajectory=False, trajectoryLabel="SSO Spacecraft")


plt.show()












