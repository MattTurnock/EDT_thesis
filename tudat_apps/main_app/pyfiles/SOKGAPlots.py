import numpy as np
import os
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import glob
import natsort

from matplotlib import pyplot as plt
import matplotlib

year = 365.25*24*60*60
AU = 1.496E11





##################################### Set some common parameters ########################################################################
normalising = False
matplotlib.rcParams.update({'font.size': 20})
plotFolder = "pyplots/SOKGA"
utils.checkFolderExist(plotFolder)
maxAltRequirement = utils.SOKGAMaxAltRequirement
# finalAltCutoff = 150

# doingJupiter=False
# doingSaturn=False
doProcessing = False
showing = False

processedSaveDir = "numpyBinaries/SOKGA/"
processedSavePath = os.path.join(processedSaveDir, "SOKGA_Processed.npy")
processedSavePathNoEDT = os.path.join(processedSaveDir, "SOKGA_Processed_NoEDT.npy")
loadingTodoList = ["bodyData", "currentData", "currentVNVData", "depVarData", "ionoData", "magData", "propData", "thrustData", "configInfo"]

if doProcessing:
    utils.checkFolderExist(processedSaveDir, emptyDirectory=True)

    SOKGAStage1JsonDirPath = os.path.join(utils.jsonInputs_dir, utils.SOKGAStage1JsonSubDir)
    SOKGAStage2JsonDirPath = os.path.join(utils.jsonInputs_dir, utils.SOKGAStage2JsonSubDir)
    SOKGAStage2NoEDTJsonDirPath = os.path.join(utils.jsonInputs_dir, utils.SOKGAStage2JsonSubDir_NoEDT)

    allStage1Jsons = utils.natsort.natsorted(utils.glob.glob(SOKGAStage1JsonDirPath + "/*"))
    allStage2Jsons = utils.natsort.natsorted(utils.glob.glob(SOKGAStage2JsonDirPath + "/*"))
    allStage2NoEDTJsons = utils.natsort.natsorted(utils.glob.glob(SOKGAStage2NoEDTJsonDirPath + "/*"))

    bar = utils.IncrementalBar("Post processing", max=len(allStage1Jsons), suffix='[%(percent).1f%%. ETA: %(eta)ds]   ')
    listOfAllNewSimDataArrays = []
    listOfAllNewSimDataArraysNoEDT = []
    # For loop to combine stage 1 and stage 2 data
    for i in range(len(allStage1Jsons)):
        finishLoop = True

        stage1JsonPath = allStage1Jsons[i]
        stage2JsonPath = allStage2Jsons[i]
        stage2NoEDTJsonPath = allStage2NoEDTJsons[i]

        saveIndex = int(stage1JsonPath.split(".")[0].split("_")[-1]) # Uses the same save index, found by string manipulation

        if "Jupiter" in stage1JsonPath:
            planet="Jupiter"
        elif "Saturn" in stage1JsonPath:
            planet="Saturn"
            # #### TEMPORARILY IGNORE SATURN ONES UNTIL FINISHED THIS CODE  ####
            # continue

        allSimDataStage2 = utils.getAllSimDataFromJson(stage2JsonPath, useCompressed=False, printInfo=True, todoList=loadingTodoList)
        allSimDataStage2NoEDT = utils.getAllSimDataFromJson(stage2NoEDTJsonPath, useCompressed=False, printInfo=True, todoList=loadingTodoList)
        # print(allSimDataStage2)
        propDataStage2 = allSimDataStage2[5]
        depVarDataStage2 = allSimDataStage2[2]
        maxAlt = np.amax(depVarDataStage2[:, 7]) / AU
        finalYear = 2000 + propDataStage2[-1, 0] / utils.year

        #Dont bother if final year does not reach the required one
        if finalYear < utils.SOKGASimulationEndYear:
            print("Removing (year): ", stage1JsonPath.split("/")[-1])
            finishLoop = False

        if maxAlt < maxAltRequirement:
            print("Removing (maxAlt = %s): " %maxAlt, stage1JsonPath.split("/")[-1])
            finishLoop = False

        if finishLoop:
            allSimDataStage1 = utils.getAllSimDataFromJson(stage1JsonPath, useCompressed=False, printInfo=True, todoList=loadingTodoList)

            newAllSimDataArray = utils.concatenateAllSimDatas(allSimDataStage1, allSimDataStage2, jsonToKeep=1)
            listOfAllNewSimDataArrays.append(newAllSimDataArray)

            newAllSimDataArrayNoEDT = utils.concatenateAllSimDatas(allSimDataStage1, allSimDataStage2NoEDT, jsonToKeep=1)
            listOfAllNewSimDataArraysNoEDT.append(newAllSimDataArrayNoEDT)


        bar.next()
    bar.finish()
    np.save(processedSavePath, np.array(listOfAllNewSimDataArrays))
    np.save(processedSavePathNoEDT, np.array(listOfAllNewSimDataArraysNoEDT))

print("Loading simdata from numpy binary")
listOfAllNewSimDataArrays = np.load(processedSavePath, allow_pickle=True)
listOfAllNewSimDataArraysNoEDT = np.load(processedSavePathNoEDT, allow_pickle=True)
# print(listOfAllNewSimDataArrays)

# for i in range(len(listOfAllNewSimDataArrays)):
#     allSimData = listOfAllNewSimDataArrays[i]
#     depVarData = allSimData[2]
#     finalAlt = depVarData[-1, 7]
#     print(finalAlt / AU)
JupiterListOfSimDataArrays = []
JupiterVelAtAltsList = []
JupiterTOFToAltsList = []




SaturnListOfSimDataArrays = []
SaturnVelAtAltsList = []
SaturnTOFToAltsList = []

bothVelAtAltsList = []
bothTOFToAltsList = []
for i in range(len(listOfAllNewSimDataArrays)):

    # Get sim data and associated parts
    allSimData = listOfAllNewSimDataArrays[i]
    jsonDataDict = allSimData[9]

    TOFAtTargetAlt, velAtTargetAlt = utils.getVelTOFTarget(allSimData, maxAltRequirement)

    # Append to Jupiter or Saturn based on json
    simFilenameBase = jsonDataDict["saveDataConfigs"]["baseFilename"]

    if "Jupiter" in simFilenameBase:
        JupiterListOfSimDataArrays.append(allSimData)
        JupiterVelAtAltsList.append(velAtTargetAlt)
        JupiterTOFToAltsList.append(TOFAtTargetAlt)
    elif "Saturn" in simFilenameBase:
        SaturnListOfSimDataArrays.append(allSimData)
        SaturnVelAtAltsList.append(velAtTargetAlt)
        SaturnTOFToAltsList.append(TOFAtTargetAlt)

    bothVelAtAltsList.append(velAtTargetAlt)
    bothTOFToAltsList.append(TOFAtTargetAlt)

# print(len(VelAtAltsList))
# print(len(TOFToAltsList))
# print(TOFAtMaxAlt)
# print(TOFToAltsList)

## Find pareto efficient trajectories
TOF_Vel_Raw = utils.arrayCoordsConvert(inputParameter1=bothTOFToAltsList, inputParameter2=bothVelAtAltsList)
efficiencyMatrix = (0, 1)
pareto_is_efficient, TOF_Vel_Pareto = utils.getParetoArray(TOF_Vel_Raw, returnBoth=True, sortOutput=True, efficiencyMatrix=efficiencyMatrix)
simDataArraysPareto = listOfAllNewSimDataArrays[pareto_is_efficient]
np.save(os.path.join(processedSaveDir, "SOKGA_simDataArrayParetoJupiter.npy"), simDataArraysPareto)
simDataArraysParetoNoEDT = listOfAllNewSimDataArraysNoEDT[pareto_is_efficient]

## Find Saturn pareto efficient trajectories
TOF_Vel_Raw_Saturn = utils.arrayCoordsConvert(inputParameter1=SaturnTOFToAltsList, inputParameter2=SaturnVelAtAltsList)
pareto_is_efficient_Saturn, TOF_Vel_Pareto_Saturn = utils.getParetoArray(TOF_Vel_Raw_Saturn, returnBoth=True, sortOutput=True, efficiencyMatrix=efficiencyMatrix)
simDataArraysPareto_Saturn = np.array(SaturnListOfSimDataArrays)[pareto_is_efficient_Saturn]

plt.figure(1, figsize=utils.figSizeDefault)
plt.scatter(np.array(JupiterTOFToAltsList) / utils.year, np.array(JupiterVelAtAltsList) / 1000, 2)
plt.scatter(np.array(SaturnTOFToAltsList) / utils.year, np.array(SaturnVelAtAltsList) / 1000, 2)
plt.scatter(TOF_Vel_Pareto[:,0] / utils.year, TOF_Vel_Pareto[:, 1] / 1000, 5, c='r')
plt.xlabel("Termination TOF [years]")
plt.ylabel("Termination Velocity [km/s]")
plt.grid()
plt.legend(["Jupiter Trajectories", "Saturn Trajectories", "Optimal Trajectories"])
plt.savefig(os.path.join(plotFolder, "SOKGA_TOF_V.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(plotFolder, "SOKGA_TOF_V.png"), bbox_inches="tight")


# print(simDataArraysPareto)


# plt.figure(1, figsize=utils.figSizeDefault)
# # plt.scatter(np.array(JupiterTOFToAltsList) / utils.year, np.array(JupiterVelAtAltsList) / 1000, 2)
#
# plt.xlabel("TOF At 100AU [years]")
# plt.ylabel("Velocity At 100AU [km/s]")
# plt.grid()




for i in range(len(listOfAllNewSimDataArrays)):
    allSimData = listOfAllNewSimDataArrays[i]
    propData = allSimData[5]
    TOFAtTargetAlt = bothTOFToAltsList[i]

    # print(propData[0,0]/year)

    if (TOFAtTargetAlt / year) < 1000:
        utils.plotTrajectoryData(propData, plotOnlyTrajectory=True, fignumber=2)

utils.plotTrajectoryData(propData, planetsToPlot=["Earth", "Jupiter", "Saturn"], plotOnlyTrajectory=False, fignumber=2, doNotPlot=True, sameScale=True, plotSun=True)


# matplotlib.rcParams.update({'font.size': 15})

#plot general efficient trajectories
for i in range(len(simDataArraysPareto)):
    allSimData = simDataArraysPareto[i]
    propData = allSimData[5]
    TOF, Vel = utils.getVelTOFTarget(allSimData, maxAltRequirement)
    print("FOR JUPITER ")
    print("Launch Year: ", propData[0,0]/year + 2000)
    print("TOF: ", TOF/year)
    print("Vel: ", Vel/1000)
    print("Initial Cartesian: ", simDataArraysPareto[0][-1]['GuidanceConfigs']['vehicleInitialCartesian'])
    print("")
    utils.plotTrajectoryData(propData, plotOnlyTrajectory=True, fignumber=3)

    allSimDataNoEDT = simDataArraysParetoNoEDT[i]
    propDataNoEDT = allSimDataNoEDT[5]
    TOFNoEDT, VelNoEDT = utils.getVelTOFTarget(allSimDataNoEDT, maxAltRequirement)
    print("FOR JUPITER (No EDT) ")
    print("Launch Year: ", propDataNoEDT[0,0]/year + 2000)
    print("TOF: ", TOFNoEDT/year)
    print("Vel: ", VelNoEDT/1000)
    print("")
    utils.plotTrajectoryData(propDataNoEDT, plotOnlyTrajectory=True, fignumber=3)

#plot saturn efficient ones
for i in range(len(simDataArraysPareto_Saturn)):
    allSimData = simDataArraysPareto_Saturn[i]
    propData = allSimData[5]
    TOF, Vel = utils.getVelTOFTarget(allSimData, maxAltRequirement)
    print("FOR Saturn ")
    print("Launch Year: ", propData[0,0]/year + 2000)
    print("TOF: ", TOF/year)
    print("Vel: ", Vel/1000)
    print("")
    utils.plotTrajectoryData(propData, plotOnlyTrajectory=True, fignumber=3)

utils.plotTrajectoryData(propData, planetsToPlot=["Earth", "Jupiter", "Saturn"], plotOnlyTrajectory=False, fignumber=3,
                         doNotPlot=True, sameScale=True, plotSun=True, legendLabelsCustom=["Jupiter GA", "Jupiter GA (No EDT)", "Saturn GA", "Earth", "Jupiter", "Saturn", "Sun", "6"],
                         saveFolder=plotFolder, savename="SOKGA_Pareto_Trajectories.pdf", savePngAndPdf=True, xlims=[-95,50], ylims=[-50, 100])

# print("Optimal trajectory launch date: %s" %(propData[0,0]/year))


# firstAllSimData = listOfAllNewSimDataArrays[10]
# utils.plotTrajectoryData(firstAllSimData[5], planetsToPlot=["Earth", "Jupiter"], plotSun=True, sameScale=True)

if showing:
    plt.show()


