import numpy as np
import os
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils

from matplotlib import pyplot as plt
import matplotlib
import natsort

matplotlib.rcParams.update({'font.size': 30})

year = 365*24*60*60
AU = 1.496E11

doingSSO = False
doingSOKGA = False
doingInO = False
doingCompare = False
doingSSOSensitivity = True

fullLoadTodoList = ["bodyData", "propData", "thrustData", "currentVNVData"]

DEFAULTSIZE = None
SaveFolder = "pyplots/finalSimsTempv2"

##################### SSO ########################

if doingSSO:

    # Plot the base case
    dataSubDir_SSO = "SSO"
    SSOJsonFilepath = "/home/matt/LinkToEDT_thesis/tudat_apps/main_app/JsonInputs/finalSims/SSO.json"
    # allSimData_SSO = utils.getAllSimDataFromFolder(dataSubDir_SSO, todoList=fullLoadTodoList, printInfo=False, jsonFilePath=SSOJsonFilepath)
    allSimData_SSO = utils.getAllSimDataFromJson(SSOJsonFilepath, todoList=fullLoadTodoList, printInfo=False)
    propData_SSO = allSimData_SSO[5]
    bodyData_SSO = allSimData_SSO[0]

    thrustData_SSO = allSimData_SSO[6]
    thrustLocal = thrustData_SSO[:,5:8]
    timesThrust = thrustData_SSO[:,0]/year
    altitudes = bodyData_SSO[:, 1]/AU

    # print("JSON: ", allSimData_SSO[9])
    # plt.figure()
    # plt.plot(timesThrust, thrustLocal[:,0])
    # plt.plot(timesThrust, thrustLocal[:,1])
    #
    # plt.figure()
    # plt.scatter(altitudes, thrustLocal[:,0])
    # plt.scatter(altitudes, thrustLocal[:,1])

    utils.plotTrajectoryData(propData_SSO, sameScale=True, planetsToPlot=["Earth"], plotSun=True, fignumber=1, plotOnlyTrajectory=False, trajectoryLabel="SSO Spacecraft", saveFolder=SaveFolder, savename="SSO-trajectory.pdf", xlims=[-10,10], ylims=[-10,10])

    utils.plotTrajectoryData(propData_SSO, fignumber=2, trajectoryLabel="BOOB", dataArrayType="propData", plotType="time-altitude", logScaleY=True, saveFolder=SaveFolder, savename="SSO-time-altitude.pdf")

    utils.plotTrajectoryData(propData_SSO, fignumber=3, trajectoryLabel="BOOB2", dataArrayType="propData", plotType="time-speed", logScaleY=True, saveFolder=SaveFolder, savename="SSO-time-speed.pdf")



######################### SOKGA ###########################################

if doingSOKGA:
    # Plot the SOKGA case
    dataSubDir_SOKGA = "SOKGA"
    allSimData_SOKGA = utils.getAllSimDataFromFolder(dataSubDir_SOKGA, todoList=fullLoadTodoList)
    propData_SOKGA = allSimData_SOKGA[5]
    bodyData_SOKGA = allSimData_SOKGA[0]

    thrustData_SOKGA = allSimData_SOKGA[6]
    # thrustLocal = thrustData_SSO[:,5:8]
    # timesThrust = thrustData_SSO[:,0]/year
    # altitudes = bodyData_SSO[:, 1]/AU

    # utils.plotTrajectoryData(propData_SOKGA, sameScale=True, planetsToPlot=["Earth", "Jupiter"], plotSun=True, fignumber=2, plotOnlyTrajectory=False, trajectoryLabel="SOKGA Spacecraft", saveFolder=SaveFolder, savename="SOKGA-test.png")



    dataSubDir_SOKGA_stage2 = "SOKGA-Stage2"
    allSimData_SOKGA_stage2 = utils.getAllSimDataFromFolder(dataSubDir_SOKGA_stage2, todoList=fullLoadTodoList)
    propData_SOKGA_stage2 = allSimData_SOKGA_stage2[5]
    bodyData_SOKGA_stage2 = allSimData_SOKGA_stage2[0]

    thrustData_SOKGA_stage2 = allSimData_SOKGA_stage2[6]
    # thrustLocal = thrustData_SSO[:,5:8]
    # timesThrust = thrustData_SSO[:,0]/year
    # altitudes = bodyData_SSO[:, 1]/AU

    propData_SOKGA_combined = np.concatenate((propData_SOKGA, propData_SOKGA_stage2))





    dataSubDir_SOKGA_reference = "SOKGA-Reference"
    allSimData_SOKGA_reference = utils.getAllSimDataFromFolder(dataSubDir_SOKGA_reference, todoList=fullLoadTodoList)
    propData_SOKGA_reference = allSimData_SOKGA_reference[5]
    bodyData_SOKGA_reference = allSimData_SOKGA_reference[0]

    thrustData_SOKGA_reference = allSimData_SOKGA_reference[6]
    # thrustLocal = thrustData_SSO[:,5:8]
    # timesThrust = thrustData_SSO[:,0]/year
    # altitudes = bodyData_SSO[:, 1]/AU





    V_reference = np.linalg.norm(propData_SOKGA_reference[-1, 4:])
    V_SOKGA = np.linalg.norm(propData_SOKGA_combined[-1, 4:])

    print("Reference V: ", V_reference)
    print("SOKGA V: ", V_SOKGA)


    customLegendTrajectory = ["SOKGA Nominal", "GA Reference", "Earth", "Jupiter", "Sun"]
    utils.plotTrajectoryData(propData_SOKGA_combined, figsize=[14,14], sameScale=True, planetsToPlot=["Earth", "Jupiter"], plotSun=True, fignumber=10, plotOnlyTrajectory=True)#, trajectoryLabel="SOKGA Spacecraft", saveFolder=SaveFolder, savename="SOKGA-test2.png")
    utils.plotTrajectoryData(propData_SOKGA_reference,  sameScale=True, planetsToPlot=["Earth", "Jupiter"], plotSun=True, fignumber=10, plotOnlyTrajectory=False, trajectoryLabel="SOKGA Spacecraft Reference", saveFolder=SaveFolder, savename="SOKGA-trajectory.pdf", legendLabelsCustom=customLegendTrajectory)


    customLegendTimeAltitude = ["SOKGA Nominal", "GA Reference"]
    utils.plotTrajectoryData(propData_SOKGA_combined, fignumber=11, trajectoryLabel="BOOB", dataArrayType="propData", plotType="time-altitude", logScaleY=True, plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_reference, fignumber=11, trajectoryLabel="PEEPEE", dataArrayType="propData", plotType="time-altitude", logScaleY=True, plotOnlyTrajectory=False, legendLabelsCustom=customLegendTimeAltitude, saveFolder=SaveFolder, savename="SOKGA-TimeAltitude.pdf")

    customLegendTimeSpeed = ["SOKGA Nominal", "GA Reference"]
    utils.plotTrajectoryData(propData_SOKGA_stage2, fignumber=12, trajectoryLabel="BOOB2", dataArrayType="propData", plotType="time-speed", logScaleY=True, plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_reference, fignumber=12, trajectoryLabel="PEEPEE2", dataArrayType="propData", plotType="time-speed", logScaleY=True, plotOnlyTrajectory=False, legendLabelsCustom=customLegendTimeSpeed, saveFolder=SaveFolder, savename="SOKGA-TimeSpeed.pdf")



##################################### InO #################################################

if doingInO:

    # Plot the InO
    dataSubDir_InO = "InO"
    allSimData_InO = utils.getAllSimDataFromFolder(dataSubDir_InO, todoList=fullLoadTodoList, useCompressed=False)
    propData_InO = allSimData_InO[5]
    bodyData_InO = allSimData_InO[0]

    thrustData_InO = allSimData_InO[6]
    thrustLocal = thrustData_InO[:,5:8]
    timesThrust = thrustData_InO[:,0]/year
    altitudes = bodyData_InO[:, 1]/AU



    print("Time of termination: ", timesThrust[-1] + 2000)



    dataSubDir_InO_stage2 = "InO-Stage2"
    allSimData_InO_stage2 = utils.getAllSimDataFromFolder(dataSubDir_InO_stage2, todoList=fullLoadTodoList, useCompressed=True)
    propData_InO_stage2 = allSimData_InO_stage2[5]
    bodyData_InO_stage2 = allSimData_InO_stage2[0]

    thrustData_InO_stage2 = allSimData_InO_stage2[6]
    thrustLocal = thrustData_InO_stage2[:,5:8]
    timesThrust = thrustData_InO_stage2[:,0]/year
    altitudes = bodyData_InO_stage2[:, 1]/AU


    propData_InO_Combined = np.concatenate((propData_InO, propData_InO_stage2))

    InOLegendLabelsCustom = ["InO Spacecraft Stage 1", "InO Spacecraft Stage 2", "Earth", "Mercury", "Sun"]
    utils.plotTrajectoryData(propData_InO, sameScale=True, planetsToPlot=["Earth", "Mercury"], plotSun=True, fignumber=30, plotOnlyTrajectory=True)#, trajectoryLabel="InO Spacecraft", saveFolder=SaveFolder, savename="InO-test.pdf")
    utils.plotTrajectoryData(propData_InO_stage2, sameScale=True, planetsToPlot=["Earth", "Mercury"], plotSun=True, fignumber=30, plotOnlyTrajectory=False, trajectoryLabel="InO Spacecraft", legendLabelsCustom=InOLegendLabelsCustom, saveFolder=SaveFolder, savename="InO-trajectory.pdf")

    utils.rescaleAndSavePlot(fignumber=30, xlims=[-10,1.2], ylims=[-4,1.5], saveFolder=SaveFolder, savename="InO-trajectory-zoomed.pdf")


    utils.plotTrajectoryData(propData_InO_Combined, fignumber=31, dataArrayType="propData", plotType="time-altitude", logScaleY=True, saveFolder=SaveFolder, savename="InO-time-altitude.pdf")

    utils.plotTrajectoryData(propData_InO_Combined, fignumber=32, dataArrayType="propData", plotType="time-speed", logScaleY=True, saveFolder=SaveFolder, savename="InO-time-speed.pdf")

if doingCompare:


    customLegendTrajectoriesCombined = ["SSO", "Nominal SOKGA", "GA Reference", "InO", "Earth", "Jupiter", "Sun"]
    utils.plotTrajectoryData(propData_SSO, fignumber=40, plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_combined, fignumber=40, plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_reference, fignumber=40, plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_InO_Combined, fignumber=40, plotOnlyTrajectory=False, sameScale=True, planetsToPlot=["Earth", "Jupiter"], plotSun=True, legendLabelsCustom=customLegendTrajectoriesCombined, saveFolder=SaveFolder, savename="Compared-trajectories.pdf")

    utils.rescaleAndSavePlot(fignumber=40, xlims=[-10,10], ylims=[-5.5,5.5], saveFolder=SaveFolder, savename="Compared-trajectories-zoomed.pdf")

    customLegendTimeAltitudeCombined = ["SSO", "Nominal SOKGA", "GA Reference", "InO"]
    utils.plotTrajectoryData(propData_SSO, fignumber=41, plotType="time-altitude", plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_combined, fignumber=41, plotType="time-altitude", plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_reference, fignumber=41, plotType="time-altitude", plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_InO_Combined, fignumber=41, plotType="time-altitude", plotOnlyTrajectory=False, logScaleY=True, legendLabelsCustom=customLegendTimeAltitudeCombined, saveFolder=SaveFolder, savename="Compared-time-altitudes.pdf")


    customLegendTimeSpeedCombined = ["SSO", "Nominal SOKGA", "GA Reference", "InO"]
    utils.plotTrajectoryData(propData_SSO, fignumber=42, plotType="time-speed", plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_combined, fignumber=42, plotType="time-speed", plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_reference, fignumber=42, plotType="time-speed", plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_InO_Combined, fignumber=42, plotType="time-speed", plotOnlyTrajectory=False, logScaleY=True, legendLabelsCustom=customLegendTimeSpeedCombined, saveFolder=SaveFolder, savename="Compared-time-speeds.pdf")

if doingSSOSensitivity:

    configPlotFolder = "pyplots/configSensitivity/"

    # typesTodo = ["currents", "diameters", "endmassMasses", "lengthRatios", "lengths", "lineSeparationCoefficients", "noLines", "occultationCoefficients", "slackCoefficients", "areaRatios"]
    baseTypeValues = [10*1E3, 1*1E-3, 100*1E-3, 0.5, 10, 0.1, 1.005, 5E-6, 0.7, 10]
    typesTodo = ["lengths", "diameters", "currents", "areaRatios", "noLines", "lengthRatios", "slackCoefficients", "lineSeparationCoefficients", "occultationCoefficients", "endmassMasses"]
    typesUnits = ["m", "m", "A", "", "", "", "", "", "", "kg"]
    listOfJsons = []
    listOfDatas = []
    rejectedJsons = []
    rejectionValue = 1

    SSOSensitivityMasterSubdir = "SSO-Configs-Sensitivity/"
    SSOSensitivityJsonDirPath = os.path.join(utils.jsonInputs_dir, "finalSims", "SSO_Config_Sensitivity")
    # allSSOSensitivitySubdirs = natsort.natsorted(os.listdir(os.path.join(utils.simulation_output_dir, SSOSensitivityMasterSubdir)))
    allSSOSensitivityJsonFiles = natsort.natsorted(os.listdir(SSOSensitivityJsonDirPath))

    # Load base / nominal case manually
    dataSetBase = utils.getAllSimDataFromJson(os.path.join(SSOSensitivityJsonDirPath, "SSO_Base.json"))
    thrustDataBase = dataSetBase[6]
    configDataBase = dataSetBase[8]
    spacecraftMassBase = configDataBase[4]
    propDataBase = dataSetBase[5]

    meanThrustBase = np.mean(thrustDataBase[:, 1])
    meanAccelBase = np.mean(thrustDataBase[:, 1])/spacecraftMassBase


    for i in range(len(typesTodo)):
        thisType = typesTodo[i]
        thisTypeJsonList = []
        thisTypeDataList = []
        # print(thisType)
        for j in range(len(allSSOSensitivityJsonFiles)):
            thisPath = os.path.join(SSOSensitivityJsonDirPath, allSSOSensitivityJsonFiles[j])

            if thisType in thisPath:
                thisTypeJsonList.append(allSSOSensitivityJsonFiles[j])
                try:
                    thisDataSet = utils.getAllSimDataFromJson(thisPath)
                except FileNotFoundError:
                    print("ERROR loading data from: %s" %thisPath)
                thisTypeDataList.append(thisDataSet)

        listOfJsons.append(thisTypeJsonList)
        listOfDatas.append(thisTypeDataList)


    # print(listOfJsons)
    bestThrustIndices = []
    bestThrustValues = []
    bestAccelIndices = []
    bestAccelValues = []
    worstThrustIndices = []
    worstThrustValues = []
    worstAccelIndices = []
    worstAccelValues = []
    for i in range(len(typesTodo)):
        thisType = typesTodo[i]
        thisJsonList = listOfJsons[i]
        thisDataList = listOfDatas[i]

        bestThrust = -999
        bestAccel = -999
        worstThrust = 1E10
        worstAccel = 1E10
        bestThrustIndex = -1
        bestAccelIndex = -1
        worstThrustIndex = -1
        worstAccelIndex = -1

        for j in range(len(thisDataList)):
            thisData = thisDataList[j]
            thisJson = thisJsonList[j]

            bodyData = thisData[0]
            propData = thisData[5]
            thrustData = thisData[6]
            currentVNVData = thisData[7]
            configData = thisData[8]

            iavgData = currentVNVData[:, 5]
            spacecraftMass = configData[4]

            meaniavg = np.mean(iavgData)
            maxiavg = np.nanmax(iavgData)
            # print(maxiavg)
            if abs(maxiavg) > rejectionValue:
                rejectedJsons.append(thisJson)
            else:
                meanThrust = np.mean(thrustData[:, 1])
                meanAccel = np.mean(thrustData[:, 1])/spacecraftMass
                if meanThrust > bestThrust:
                    bestThrust = meanThrust
                    bestThrustIndex = j
                if meanAccel > bestAccel:
                    bestAccel = meanAccel
                    bestAccelIndex = j
                if meanThrust < worstThrust:
                    worstThrust = meanThrust
                    worstThrustIndex = j
                if meanAccel < worstAccel:
                    worstAccel = meanAccel
                    worstAccelIndex = j

                # print("MASS: ", spacecraftMass)
                # print("THRUST: ", meanThrust)
                # print("ACCEL: ", meanAccel, "\n")


        bestThrustIndices.append(bestThrustIndex)
        bestThrustValues.append(bestThrust)

        bestAccelIndices.append(bestAccelIndex)
        bestAccelValues.append(bestAccel)

        worstThrustIndices.append(worstThrustIndex)
        worstThrustValues.append(worstThrust)

        worstAccelIndices.append(worstAccelIndex)
        worstAccelValues.append(worstAccel)

    # thrustAccelIndices = [bestThrustIndices, bestAccelIndices, worstThrustIndices, worstAccelIndices]
    # print(bestThrustIndices)
    print("Rejected: ", rejectedJsons)



    for i in range(len(bestThrustIndices)):

        thisParameterType = typesTodo[i]
        thisParameterValueBase = baseTypeValues[i]
        thisRunnerValuesList = utils.configSensitivityRunnerValues[i]

        thisBestThrustIndex = bestThrustIndices[i]
        thisBestThrustValue = bestThrustValues[i]
        thisBestThrustJson = listOfJsons[i][thisBestThrustIndex]
        thisBestThrustAllData = listOfDatas[i][thisBestThrustIndex]
        thisBestThrust_ThrustData = thisBestThrustAllData[6]
        thisBestThrust_CurrentVNVData = thisBestThrustAllData[7]
        thisBestThrust_CurrentData = thisBestThrustAllData[1]
        thisBestThrust_MagfieldData = thisBestThrustAllData[4]
        thisBestThrustParameterValue = thisRunnerValuesList[thisBestThrustIndex]

        thisWorstThrustIndex = worstThrustIndices[i]
        thisWorstThrustValue = worstThrustValues[i]
        thisWorstThrustJson = listOfJsons[i][thisWorstThrustIndex]
        thisWorstThrustAllData = listOfDatas[i][thisWorstThrustIndex]
        thisWorstThrust_ThrustData = thisWorstThrustAllData[6]
        thisWorstThrustParameterValue = thisRunnerValuesList[thisWorstThrustIndex]





        # Set plot legend and general parameters
        thisThrustLegend = ["Best Thrust: %0.4g nN, %0.4g" %(1E9* thisBestThrustValue, thisBestThrustParameterValue),
                            "Worst Thrust: %0.4g nN, %0.4g" %(1E9* thisWorstThrustValue, thisWorstThrustParameterValue),
                            "Base Thrust: %0.4g nN, %0.4g" %(1E9* meanThrustBase, thisParameterValueBase)]
        scatterDotSize = 5
        plt.figure(figsize=2* utils.figSizeDefault)
        ax = plt.gca()
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        plt.xlabel("Year")
        plt.ylabel("Thrust [nN]")
        plt.grid()
        # plt.title(thisParameterType)


        # plt.plot(thisBestThrust_CurrentVNVData[:,0]/year + 2000, thisBestThrust_CurrentVNVData[:, -1])                               # Shows whether thrust is on or off over time
        # plt.plot(thisBestThrust_CurrentData[:,0]/year + 2000, thisBestThrust_CurrentData[:, 1])                                      # Shows overall current magnitude over time
        # plt.plot(thisBestThrust_CurrentVNVData[:,0]/year + 2000, thisBestThrust_CurrentVNVData[:, 13]*1E-9)                            # Shows direction of spacecraft over time
        # plt.plot(thisBestThrust_MagfieldData[:,0]/year + 2000, thisBestThrust_MagfieldData[:, 1])                                    # Shows overall magnetic field strenght over time

        # Do scatter plots - best worst and base
        plt.scatter(thisBestThrust_ThrustData[:,0] /year + 2000, 1E9* thisBestThrust_ThrustData[:,1], c="C0", s=scatterDotSize )                                          # Shows thrust data over time
        plt.scatter(thisWorstThrust_ThrustData[:,0] /year + 2000, 1E9* thisWorstThrust_ThrustData[:,1], c="C1", s=scatterDotSize)                                   # Shows worst thrust data over time
        plt.scatter(thrustDataBase[:, 0] / year + 2000, 1E9* thrustDataBase[:,1], c="C2", s=scatterDotSize)                                                         # Shows nominal thrust data over time

        # Do hlines for mean values
        plt.axhline(1E9* thisBestThrustValue, c="C0")
        plt.axhline(1E9* thisWorstThrustValue, c="C1")
        plt.axhline(1E9* meanThrustBase, c="C2")

        # Finalise the plot with legend and saving
        plt.legend(thisThrustLegend)

        plt.savefig(os.path.join(configPlotFolder, "%s-thrustPlot.png" %thisParameterType))



        ################# ACCELERATION PLOT ##################################################



        thisBestAccelIndex = bestAccelIndices[i]
        thisBestAccelValue = bestAccelValues[i]
        thisBestAccelJson = listOfJsons[i][thisBestAccelIndex]
        thisBestAccelAllData = listOfDatas[i][thisBestAccelIndex]
        thisBestAccel_ThrustData = thisBestAccelAllData[6]
        thisBestAccel_CurrentVNVData = thisBestAccelAllData[7]
        thisBestAccel_CurrentData = thisBestAccelAllData[1]
        thisBestAccel_MagfieldData = thisBestAccelAllData[4]
        thisBestAccel_ConfigInfo = thisBestAccelAllData[8]
        thisBestAccel_SpacecraftMass = thisBestAccel_ConfigInfo[4]
        thisBestAccelParameterValue = thisRunnerValuesList[thisBestAccelIndex]

        thisWorstAccelIndex = worstAccelIndices[i]
        thisWorstAccelValue = worstAccelValues[i]
        thisWorstAccelJson = listOfJsons[i][thisWorstAccelIndex]
        thisWorstAccelAllData = listOfDatas[i][thisWorstAccelIndex]
        thisWorstAccel_ThrustData = thisWorstAccelAllData[6]
        thisWorstAccel_ConfigInfo = thisWorstAccelAllData[8]
        thisWorstAccel_SpacecraftMass = thisWorstAccel_ConfigInfo[4]
        thisWorstAccelParameterValue = thisRunnerValuesList[thisWorstAccelIndex]





        # Set plot legend and general parameters
        thisThrustLegend = ["Best Accel: %0.4g m/s^2, %0.4g" %(thisBestAccelValue, thisBestAccelParameterValue),
                            "Worst Accel: %0.4g m/s^2, %0.4g" %(thisWorstAccelValue, thisWorstAccelParameterValue),
                            "Base Accel: %0.4g m/s^2, %0.4g" %(meanAccelBase, thisParameterValueBase)]
        scatterDotSize = 5
        plt.figure(figsize=2* utils.figSizeDefault)
        ax = plt.gca()
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        plt.xlabel("Year")
        plt.ylabel("Acceleration [m/s^2]")
        plt.grid()
        # plt.title(thisParameterType)


        # plt.plot(thisBestThrust_CurrentVNVData[:,0]/year + 2000, thisBestThrust_CurrentVNVData[:, -1])                               # Shows whether thrust is on or off over time
        # plt.plot(thisBestThrust_CurrentData[:,0]/year + 2000, thisBestThrust_CurrentData[:, 1])                                      # Shows overall current magnitude over time
        # plt.plot(thisBestThrust_CurrentVNVData[:,0]/year + 2000, thisBestThrust_CurrentVNVData[:, 13]*1E-9)                            # Shows direction of spacecraft over time
        # plt.plot(thisBestThrust_MagfieldData[:,0]/year + 2000, thisBestThrust_MagfieldData[:, 1])                                    # Shows overall magnetic field strenght over time

        # Do scatter plots - best worst and base
        plt.scatter(thisBestAccel_ThrustData[:,0] /year + 2000, (thisBestAccel_ThrustData[:,1])/thisBestAccel_SpacecraftMass, c="C0", s=scatterDotSize )                                          # Shows thrust data over time
        plt.scatter(thisWorstAccel_ThrustData[:,0] /year + 2000, (thisWorstAccel_ThrustData[:,1])/thisWorstAccel_SpacecraftMass, c="C1", s=scatterDotSize)                                   # Shows worst thrust data over time
        plt.scatter(thrustDataBase[:, 0] / year + 2000, thrustDataBase[:,1]/spacecraftMassBase, c="C2", s=scatterDotSize)                                                         # Shows nominal thrust data over time

        # Do hlines for mean values
        plt.axhline(thisBestAccelValue, c="C0")
        plt.axhline(thisWorstAccelValue, c="C1")
        plt.axhline(meanAccelBase, c="C2")

        # Finalise the plot with legend and saving
        plt.legend(thisThrustLegend)

        plt.savefig(os.path.join(configPlotFolder, "%s-accelPlot.png" %thisParameterType))

    
    #plot trajectory of base case
    utils.plotTrajectoryData(propDataBase, plotType="x-y")


plt.show()