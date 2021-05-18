import numpy as np
import os
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils

from matplotlib import pyplot as plt
import matplotlib

year = 365*24*60*60
AU = 1.496E11

##################################### Set some common parameters ########################################################################
matplotlib.rcParams.update({'font.size': 20})
plotFolder = "pyplots/testCasePlots"

testCaseJsonPath = os.path.join(utils.jsonInputs_dir, "testVariables.json")
testCaseAllData = utils.getAllSimDataFromJson(testCaseJsonPath)

testCaseBodyData       = testCaseAllData[0]
testCaseCurrentData    = testCaseAllData[1]
testCaseDepVarData     = testCaseAllData[2]
testCaseIonoData       = testCaseAllData[3]
testCaseMagData        = testCaseAllData[4]
testCasePropData       = testCaseAllData[5]
testCaseThrustData     = testCaseAllData[6]
testCaseCurrentVNVData = testCaseAllData[7]
testCaseConfigInfoData = testCaseAllData[8]
testCaseJsonData       = testCaseAllData[9]

plotExtraKepler = False
plotExtraTrajectory = False
plotSeparations = True

utils.plotTrajectoryData(testCasePropData, dataArrayType="propData", plotType="x-y", plotSun=True, sameScale=True, plotTitle="Trajectory", planetsToPlot=["Jupiter"])
if plotExtraTrajectory:
    utils.plotTrajectoryData(testCasePropData, dataArrayType="propData", plotType="time-altitude", plotTitle="time-altitude")
    utils.plotTrajectoryData(testCasePropData, dataArrayType="propData", plotType="time-speed", plotTitle="time-speed")
    utils.plotTrajectoryData(testCasePropData, dataArrayType="propData", plotType="altitude-speed", plotTitle="altitude-speed")


utils.plotTrajectoryData(testCaseDepVarData, dataArrayType="depVarData", plotType="time-SMA", plotTitle="SMA")
utils.plotTrajectoryData(testCaseDepVarData, dataArrayType="depVarData", plotType="time-ECC", plotTitle="ECC")
if plotExtraKepler:
    utils.plotTrajectoryData(testCaseDepVarData, dataArrayType="depVarData", plotType="time-INC", plotTitle="INC")
    utils.plotTrajectoryData(testCaseDepVarData, dataArrayType="depVarData", plotType="time-AOP", plotTitle="AOP")
    utils.plotTrajectoryData(testCaseDepVarData, dataArrayType="depVarData", plotType="time-RAAN", plotTitle="RAAN")
    utils.plotTrajectoryData(testCaseDepVarData, dataArrayType="depVarData", plotType="time-TA", plotTitle="TA")

if plotSeparations:
    jupiterSeparationArray = np.empty([len(testCasePropData[:,0]), 2])

    JupiterCoordinateDifferenceArray = testCasePropData[:, 1:4] - testCaseDepVarData[:, 11:14]
    print(JupiterCoordinateDifferenceArray)

    jupiterSeparationArray = np.empty([len(testCasePropData[:,0]), 2])
    jupiterSeparationArray[:, 0] = testCasePropData[:, 0]
    jupiterSeparationArray[:, 1] = np.linalg.norm(JupiterCoordinateDifferenceArray, axis=1)

    utils.plotTrajectoryData(jupiterSeparationArray, dataArrayType="separation", plotType="time-separation", plotTitle="Vehicle Jupiter Separation")
    utils.plotTrajectoryData(testCaseDepVarData, dataArrayType="depVarData", plotType="x-y-jupiter", plotTitle="Jupiter trajectory plot", sameScale=True, plotSun=True)
    # print(jupiterSeparationArray)




plt.show()
