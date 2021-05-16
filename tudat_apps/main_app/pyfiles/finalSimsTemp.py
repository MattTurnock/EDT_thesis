import numpy as np
import os
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils

from matplotlib import pyplot as plt
import matplotlib
import natsort


matplotlib.rcParams.update({'font.size': 50,
                            "lines.linewidth": 3})

thisFigureSize = 3*np.array([12, 10])

year = 365*24*60*60
AU = 1.496E11

doingSSO = False
doingSOKGA = False
doingInO = False
doingCompare = False
doingSSOSensitivity = True

showing = False

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

    doExample = False
    plottingVsTime = False
    plottingVsParameterValue = True


    configPlotFolder = "pyplots/configSensitivity/"
    sensitivityGeneralFolderName = "SSO-Configs-Sensitivity"
    if doExample:
        SSOSensitivityMasterSubdir = "SSO-Configs-Sensitivity_EXAMPLE/"
        SSOSensitivityJsonDirPath = os.path.join(utils.jsonInputs_dir, "finalSims", "SSO_Config_Sensitivity_EXAMPLE")
    else:
        SSOSensitivityMasterSubdir = "SSO-Configs-Sensitivity/"
        SSOSensitivityJsonDirPath = os.path.join(utils.jsonInputs_dir, "finalSims", "SSO_Config_Sensitivity")

    # typesTodo = ["currents", "diameters", "endmassMasses", "lengthRatios", "lengths", "lineSeparationCoefficients", "noLines", "occultationCoefficients", "slackCoefficients", "areaRatios"]
    baseTypeValues =      [10*1E3,    10*1E-3,      100,    0.5,           10,                        0.5,                            10,              0.7]
    typesTodo =           ["lengths", "diameters", "currents", "areaRatios",  "noLines",                 "lengthRatios",                 "endmassMasses", "rotationCoefficients"]
    typesTodoPlotTitles = ["Length",  "Diameter",  "Current",  "Area Ratio",  "Number of Primary Lines", "Primary Segment Length Ratio", "Endmass Mass",  "Rotation Coefficient"]
    typesUnits =          ["m",       "m",         "mA",        "-",           "-",                       "-",                            "kg",            "-"]
    logPlotTypesThrust =  [                        "currents",                "noLines"                                                                                         ]
    logPlotTypesAccel  =  ["lengths", "diameters", "currents",                "noLines"                                                                                         ]
    listOfJsons = []
    listOfDatas = []
    rejectedJsons = []
    rejectedParameterValues = []
    rejectediValues = []
    rejectionValue = 0.1


    # allSSOSensitivitySubdirs = natsort.natsorted(os.listdir(os.path.join(utils.simulation_output_dir, SSOSensitivityMasterSubdir)))
    allSSOSensitivityJsonFiles = natsort.natsorted(os.listdir(SSOSensitivityJsonDirPath))

    # Load base / nominal case manually
    dataSetBase = utils.getAllSimDataFromJson(os.path.join(SSOSensitivityJsonDirPath, "SSO_Bare_Base.json"))
    thrustDataBase = dataSetBase[6]
    configDataBase = dataSetBase[8]
    spacecraftMassBase = configDataBase[4]
    propDataBase = dataSetBase[5]

    meanThrustBase = np.mean(thrustDataBase[:, 1])
    meanAccelBase = np.mean(thrustDataBase[:, 1])/spacecraftMassBase

    bar = utils.IncrementalBar("Loading data ", max=len(typesTodo), suffix='[%(percent).1f%%. ETA: %(eta)ds]   ')
    for i in range(len(typesTodo)):
        thisType = typesTodo[i]
        thisTypeJsonList = []
        thisTypeDataList = []
        print("Loading type: %s." %(thisType))

        for j in range(len(allSSOSensitivityJsonFiles)):
            thisPath = os.path.join(SSOSensitivityJsonDirPath, allSSOSensitivityJsonFiles[j])

            if thisType in thisPath:
                thisTypeJsonList.append(allSSOSensitivityJsonFiles[j])
                try:
                    thisDataSet = utils.getAllSimDataFromJson(thisPath)
                    thisTypeDataList.append(thisDataSet)
                except FileNotFoundError:
                    print("ERROR loading data from: %s" %thisPath)


        listOfJsons.append(thisTypeJsonList)
        listOfDatas.append(thisTypeDataList)
        bar.next()
    bar.finish()

    configTypesToTest = ["Bare", "Trans"]

    foundDataBare = {"bestThrustIndices": [],
                    "bestThrustValues": [],
                    "bestAccelIndices": [],
                    "bestAccelValues": [],
                    "worstThrustIndices": [],
                    "worstThrustValues": [],
                    "worstAccelIndices": [],
                    "worstAccelValues": []}

    foundDataTrans = {"bestThrustIndices": [],
                     "bestThrustValues": [],
                     "bestAccelIndices": [],
                     "bestAccelValues": [],
                     "worstThrustIndices": [],
                     "worstThrustValues": [],
                     "worstAccelIndices": [],
                     "worstAccelValues": []}

    listOfMeanThrustsBare = []
    listOfMeanAccelsBare = []

    listOfMeanThrustsTrans = []
    listOfMeanAccelsTrans = []

    configSensitivityRunnerValuesToPlotBare = []
    configSensitivityRunnerValuesToPlotTrans = []

    for h in range(len(configTypesToTest)):
        thisConfigType = configTypesToTest[h]
        if thisConfigType == "Bare":
            libraryToWorkWith = foundDataBare

            listOfMeanThrusts = listOfMeanThrustsBare
            listOfMeanAccels = listOfMeanAccelsBare

            configSensitivityRunnerValuesToPlot = configSensitivityRunnerValuesToPlotBare

        elif thisConfigType == "Trans":
            libraryToWorkWith = foundDataTrans

            listOfMeanThrusts = listOfMeanThrustsTrans
            listOfMeanAccels = listOfMeanAccelsTrans

            configSensitivityRunnerValuesToPlot = configSensitivityRunnerValuesToPlotTrans

        # print(listOfJsons)
        # bestThrustIndices = []
        # bestThrustValues = []
        # bestAccelIndices = []
        # bestAccelValues = []
        # worstThrustIndices = []
        # worstThrustValues = []
        # worstAccelIndices = []
        # worstAccelValues = []


        for i in range(len(typesTodo)):
            thisType = typesTodo[i]
            thisJsonList = listOfJsons[i]
            thisDataList = listOfDatas[i]
            thisRunnerValuesList = utils.configSensitivityRunnerValues[i]
            # print("DATA: ", len(thisDataList))
            # print("RUNNER: ", len(thisRunnerValuesList))

            bestThrust = -999
            bestAccel = -999
            worstThrust = 1E10
            worstAccel = 1E10
            bestThrustIndex = -1
            bestAccelIndex = -1
            worstThrustIndex = -1
            worstAccelIndex = -1

            thisMeanThrusts = []
            thisMeanAccels = []

            thisRunnerValuesListToPlot = []

            for j in range(len(thisDataList)):
                thisData = thisDataList[j]
                thisJson = thisJsonList[j]
                if "Bare" in thisJson:
                    thisRunnerValue = thisRunnerValuesList[j]
                elif "Trans" in thisJson:
                    thisRunnerValue = thisRunnerValuesList[j - len(thisRunnerValuesList)]


                if thisConfigType in thisJson:

                    bodyData = thisData[0]
                    propData = thisData[5]
                    thrustData = thisData[6]
                    currentVNVData = thisData[7]
                    configData = thisData[8]

                    iavgData = currentVNVData[:, 5]
                    spacecraftMass = configData[4]

                    meaniavg = np.mean(iavgData)
                    maxiavg = np.nanmax(iavgData)

                    # meanThrust = np.mean(thrustData[:, 1])
                    # meanAccel = np.mean(thrustData[:, 1])/spacecraftMass
                    # thisMeanThrusts.append(meanThrust)
                    # thisMeanAccels.append(meanAccel)
                    # print(maxiavg)
                    if abs(maxiavg) > rejectionValue:
                        # Directly append json and i values to lists
                        rejectedJsons.append(thisJson)
                        rejectediValues.append(maxiavg)
                        # rejectedParameterValues.append()

                        with open(os.path.join(SSOSensitivityJsonDirPath, thisJson), 'r+') as f:
                            allJsonVariables = utils.json.load(f)

                            if "Bare_diameters" in thisJson:
                                parameterValueToAppend = allJsonVariables["EDTConfigs"]["tetherDiameter"]
                            elif "Bare_currents" in thisJson:
                                parameterValueToAppend = allJsonVariables["EDTConfigs"]["emitterCurrentmA"]
                            elif "Bare_noLines" in thisJson:
                                parameterValueToAppend = allJsonVariables["EDTConfigs"]["hoytether"]["noPrimaryLines"]
                            else:
                                parameterValueToAppend = float("nan")
                                print("WARNING: Rejected values unable to be appended, replacing with nan")

                            rejectedParameterValues.append(parameterValueToAppend)

                    else:

                        meanThrust = np.mean(thrustData[:, 1])
                        meanAccel = np.mean(thrustData[:, 1])/spacecraftMass
                        thisMeanThrusts.append(meanThrust)
                        thisMeanAccels.append(meanAccel)

                        thisRunnerValuesListToPlot.append(thisRunnerValue)

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
                            # print(j)

            configSensitivityRunnerValuesToPlot.append(thisRunnerValuesListToPlot)

                    # print("MASS: ", spacecraftMass)
                    # print("THRUST: ", meanThrust)
                    # print("ACCEL: ", meanAccel, "\n")


            libraryToWorkWith["bestThrustIndices"].append(bestThrustIndex)
            libraryToWorkWith["bestThrustValues"].append(bestThrust)

            libraryToWorkWith["bestAccelIndices"].append(bestAccelIndex)
            libraryToWorkWith["bestAccelValues"].append(bestAccel)

            libraryToWorkWith["worstThrustIndices"].append(worstThrustIndex)
            libraryToWorkWith["worstThrustValues"].append(worstThrust)

            libraryToWorkWith["worstAccelIndices"].append(worstAccelIndex)
            libraryToWorkWith["worstAccelValues"].append(worstAccel)

            listOfMeanThrusts.append(thisMeanThrusts)
            listOfMeanAccels.append(thisMeanAccels)




    listOfMeanThrustsBare = np.array(listOfMeanThrustsBare)
    listOfMeanAccelsBare = np.array(listOfMeanAccelsBare)

    listOfMeanThrustsTrans = np.array(listOfMeanThrustsTrans)
    listOfMeanAccelsTrans = np.array(listOfMeanAccelsTrans)




    # thrustAccelIndices = [bestThrustIndices, bestAccelIndices, worstThrustIndices, worstAccelIndices]
    # print(bestThrustIndices)
    print("Rejected: ", rejectedJsons)
    print("Rejected Values: ", rejectedParameterValues)
    print("Rejected i values: ", rejectediValues)
    rejectedDataToSave = np.empty([len(rejectedJsons), 3], dtype='S100')
    rejectedDataToSave[:, 0] = rejectedJsons
    rejectedDataToSave[:, 1] = rejectedParameterValues
    rejectedDataToSave[:, 2] = rejectediValues
    # print(rejectedDataToSave)
    np.savetxt(os.path.join(configPlotFolder, "rejectedData.txt"), rejectedDataToSave, fmt="%s", delimiter=",\t")

    for i in range(len(typesTodo)):

        thisParameterType = typesTodo[i]
        thisParameterUnits = typesUnits[i]
        if doExample:
            thisRunnerValuesList = [utils.configSensitivityRunnerValues[i][0], utils.configSensitivityRunnerValues[i][-1]]
        else:
            thisRunnerValuesList = utils.configSensitivityRunnerValues[i]


        if plottingVsTime:

            ##############################################################################################################
            ############################### DO PLOTTING FOR BEST / WORST OVER TIME ######################################
            ##############################################################################################################

            # Do bare tether ones
            # print(foundDataBare["bestThrustIndices"])
            thisBestThrustIndex_Bare = foundDataBare["bestThrustIndices"][i]
            thisBestThrustValue_Bare = foundDataBare["bestThrustValues"][i]
            thisBestThrustJson_Bare = listOfJsons[i][thisBestThrustIndex_Bare]
            thisBestThrustAllData_Bare = listOfDatas[i][thisBestThrustIndex_Bare]
            thisBestThrust_ThrustData_Bare = thisBestThrustAllData_Bare[6]
            thisBestThrustParameterValue_Bare = thisRunnerValuesList[thisBestThrustIndex_Bare]

            thisWorstThrustIndex_Bare = foundDataBare["worstThrustIndices"][i]
            thisWorstThrustValue_Bare = foundDataBare["worstThrustValues"][i]
            thisWorstThrustJson_Bare = listOfJsons[i][thisWorstThrustIndex_Bare]
            thisWorstThrustAllData_Bare = listOfDatas[i][thisWorstThrustIndex_Bare]
            thisWorstThrust_ThrustData_Bare = thisWorstThrustAllData_Bare[6]
            thisWorstThrustParameterValue_Bare = thisRunnerValuesList[thisWorstThrustIndex_Bare]

            # Do Trans ones
            thisBestThrustIndex_Trans = foundDataTrans["bestThrustIndices"][i]
            thisBestThrustValue_Trans = foundDataTrans["bestThrustValues"][i]
            thisBestThrustJson_Trans = listOfJsons[i][thisBestThrustIndex_Trans]
            thisBestThrustAllData_Trans = listOfDatas[i][thisBestThrustIndex_Trans]
            thisBestThrust_ThrustData_Trans = thisBestThrustAllData_Trans[6]
            thisBestThrustParameterValue_Trans = thisRunnerValuesList[thisBestThrustIndex_Trans - len(thisRunnerValuesList)] # Minus accounts for that all trans values are 10 larger

            thisWorstThrustIndex_Trans = foundDataTrans["worstThrustIndices"][i]
            thisWorstThrustValue_Trans = foundDataTrans["worstThrustValues"][i]
            thisWorstThrustJson_Trans = listOfJsons[i][thisWorstThrustIndex_Trans]
            thisWorstThrustAllData_Trans = listOfDatas[i][thisWorstThrustIndex_Trans]
            thisWorstThrust_ThrustData_Trans = thisWorstThrustAllData_Trans[6]
            thisWorstThrustParameterValue_Trans = thisRunnerValuesList[thisWorstThrustIndex_Trans - len(thisRunnerValuesList)] # Minus accounts for that all trans values are 10 larger





            # Set plot legend and general parameters
            thisThrustLegend = ["Best Thrust Bare: %0.4g nN, %0.4g" %(1E9* thisBestThrustValue_Bare, thisBestThrustParameterValue_Bare),
                                "Worst Thrust Bare: %0.4g nN, %0.4g" %(1E9* thisWorstThrustValue_Bare, thisWorstThrustParameterValue_Bare),
                                "Base Thrust: %0.4g nN, %0.4g" %(1E9* meanThrustBase, baseTypeValues[i]),
                                "Best Thrust Trans: %0.4g nN, %0.4g" %(1E9* thisBestThrustValue_Trans, thisBestThrustParameterValue_Trans),
                                "Worst Thrust Trans: %0.4g nN, %0.4g" %(1E9* thisWorstThrustValue_Trans, thisWorstThrustParameterValue_Trans)]
            scatterDotSize = 5
            transLinestyle = "dashed"
            plt.figure(figsize=thisFigureSize)
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

            # Do scatter plots - best worst and base, then again for trans values
            plt.plot(thisBestThrust_ThrustData_Bare[1:, 0] /year + 2000, 1E9* thisBestThrust_ThrustData_Bare[1:, 1], c="C0")#, s=scatterDotSize )                                          # Shows best thrust data over time
            plt.plot(thisWorstThrust_ThrustData_Bare[1:, 0] /year + 2000, 1E9* thisWorstThrust_ThrustData_Bare[1:, 1], c="C1")#, s=scatterDotSize)                                   # Shows worst thrust data over time
            plt.plot(thrustDataBase[1:, 0] / year + 2000, 1E9* thrustDataBase[1:,1], c="C2")#, s=scatterDotSize)                                                         # Shows nominal thrust data over time

            plt.plot(thisBestThrust_ThrustData_Trans[1:, 0] /year + 2000, 1E9* thisBestThrust_ThrustData_Trans[1:, 1], c="C3", linestyle=transLinestyle)#, s=scatterDotSize )                                          # Shows best thrust data over time
            plt.plot(thisWorstThrust_ThrustData_Trans[1:, 0] /year + 2000, 1E9* thisWorstThrust_ThrustData_Trans[1:, 1], c="C4", linestyle=transLinestyle)#, s=scatterDotSize)                                   # Shows worst thrust data over time

            # Do hlines for mean values
            plt.axhline(1E9* thisBestThrustValue_Bare, c="C0")
            plt.axhline(1E9* thisWorstThrustValue_Bare, c="C1")
            plt.axhline(1E9* meanThrustBase, c="C2")

            plt.axhline(1E9* thisBestThrustValue_Trans, c="C3", linestyle=transLinestyle)
            plt.axhline(1E9* thisWorstThrustValue_Trans, c="C4", linestyle=transLinestyle)

            # Finalise the plot with legend and saving
            plt.legend(thisThrustLegend)

            plt.savefig(os.path.join(configPlotFolder, "time-thrust-%s.png" %thisParameterType))



            ################# ACCELERATION PLOT ##################################################


            # Data for Bare
            thisBestAccelIndex_Bare = foundDataBare["bestAccelIndices"][i]
            thisBestAccelValue_Bare = foundDataBare["bestAccelValues"][i]
            thisBestAccelJson_Bare = listOfJsons[i][thisBestAccelIndex_Bare]
            thisBestAccelAllData_Bare = listOfDatas[i][thisBestAccelIndex_Bare]
            thisBestAccel_ThrustData_Bare = thisBestAccelAllData_Bare[6]
            thisBestAccel_ConfigInfo_Bare = thisBestAccelAllData_Bare[8]
            thisBestAccel_SpacecraftMass_Bare = thisBestAccel_ConfigInfo_Bare[4]
            thisBestAccelParameterValue_Bare = thisRunnerValuesList[thisBestAccelIndex_Bare]

            thisWorstAccelIndex_Bare = foundDataBare["worstAccelIndices"][i]
            thisWorstAccelValue_Bare = foundDataBare["worstAccelValues"][i]
            thisWorstAccelJson_Bare = listOfJsons[i][thisWorstAccelIndex_Bare]
            thisWorstAccelAllData_Bare = listOfDatas[i][thisWorstAccelIndex_Bare]
            thisWorstAccel_ThrustData_Bare = thisWorstAccelAllData_Bare[6]
            thisWorstAccel_ConfigInfo_Bare = thisWorstAccelAllData_Bare[8]
            thisWorstAccel_SpacecraftMass_Bare = thisWorstAccel_ConfigInfo_Bare[4]
            thisWorstAccelParameterValue_Bare = thisRunnerValuesList[thisWorstAccelIndex_Bare]


            # Data for Trans
            thisBestAccelIndex_Trans = foundDataTrans["bestAccelIndices"][i]
            thisBestAccelValue_Trans = foundDataTrans["bestAccelValues"][i]
            thisBestAccelJson_Trans = listOfJsons[i][thisBestAccelIndex_Trans]
            thisBestAccelAllData_Trans = listOfDatas[i][thisBestAccelIndex_Trans]
            thisBestAccel_ThrustData_Trans = thisBestAccelAllData_Trans[6]
            thisBestAccel_ConfigInfo_Trans = thisBestAccelAllData_Trans[8]
            thisBestAccel_SpacecraftMass_Trans = thisBestAccel_ConfigInfo_Trans[4]
            thisBestAccelParameterValue_Trans = thisRunnerValuesList[thisBestAccelIndex_Trans - len(thisRunnerValuesList)] # Minus accounts for that all trans values are 10 larger

            thisWorstAccelIndex_Trans = foundDataTrans["worstAccelIndices"][i]
            thisWorstAccelValue_Trans = foundDataTrans["worstAccelValues"][i]
            thisWorstAccelJson_Trans = listOfJsons[i][thisWorstAccelIndex_Trans]
            thisWorstAccelAllData_Trans = listOfDatas[i][thisWorstAccelIndex_Trans]
            thisWorstAccel_ThrustData_Trans = thisWorstAccelAllData_Trans[6]
            thisWorstAccel_ConfigInfo_Trans = thisWorstAccelAllData_Trans[8]
            thisWorstAccel_SpacecraftMass_Trans = thisWorstAccel_ConfigInfo_Trans[4]
            thisWorstAccelParameterValue_Trans = thisRunnerValuesList[thisWorstAccelIndex_Trans - len(thisRunnerValuesList)] # Minus accounts for that all trans values are 10 larger






            # Set plot legend and general parameters
            thisAccelLegend = ["Best Accel Bare: %0.4g m/s^2, %0.4g" %(thisBestAccelValue_Bare, thisBestAccelParameterValue_Bare),
                                "Worst Accel Bare: %0.4g m/s^2, %0.4g" %(thisWorstAccelValue_Bare, thisWorstAccelParameterValue_Bare),
                                "Base Accel: %0.4g m/s^2, %0.4g" %(meanAccelBase, baseTypeValues[i]),
                               "Best Accel Trans: %0.4g m/s^2, %0.4g" %(thisBestAccelValue_Trans, thisBestAccelParameterValue_Trans),
                               "Worst Accel Trans: %0.4g m/s^2, %0.4g" %(thisWorstAccelValue_Trans, thisWorstAccelParameterValue_Trans)]
            scatterDotSize = 5
            plt.figure(figsize=thisFigureSize)
            ax = plt.gca()
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            plt.xlabel("Year")
            plt.ylabel("Acceleration [m/s$^2$]")
            plt.grid(which="both")
            # plt.title(thisParameterType)


            # plt.plot(thisBestThrust_CurrentVNVData[:,0]/year + 2000, thisBestThrust_CurrentVNVData[:, -1])                               # Shows whether thrust is on or off over time
            # plt.plot(thisBestThrust_CurrentData[:,0]/year + 2000, thisBestThrust_CurrentData[:, 1])                                      # Shows overall current magnitude over time
            # plt.plot(thisBestThrust_CurrentVNVData[:,0]/year + 2000, thisBestThrust_CurrentVNVData[:, 13]*1E-9)                            # Shows direction of spacecraft over time
            # plt.plot(thisBestThrust_MagfieldData[:,0]/year + 2000, thisBestThrust_MagfieldData[:, 1])                                    # Shows overall magnetic field strenght over time

            # Do scatter plots - best worst and base
            plt.plot(thisBestAccel_ThrustData_Bare[1:,0] /year + 2000, (thisBestAccel_ThrustData_Bare[1:,1])/thisBestAccel_SpacecraftMass_Bare, c="C0")#, s=scatterDotSize )                                          # Shows thrust data over time
            plt.plot(thisWorstAccel_ThrustData_Bare[1:,0] /year + 2000, (thisWorstAccel_ThrustData_Bare[1:,1])/thisWorstAccel_SpacecraftMass_Bare, c="C1")#, s=scatterDotSize)                                   # Shows worst thrust data over time
            plt.plot(thrustDataBase[1:, 0] / year + 2000, thrustDataBase[1:,1]/spacecraftMassBase, c="C2")#, s=scatterDotSize)                                                         # Shows nominal thrust data over time

            plt.plot(thisBestAccel_ThrustData_Trans[1:,0] /year + 2000, (thisBestAccel_ThrustData_Trans[1:,1])/thisBestAccel_SpacecraftMass_Trans, c="C3", linestyle=transLinestyle)#, s=scatterDotSize )                                          # Shows thrust data over time
            plt.plot(thisWorstAccel_ThrustData_Trans[1:,0] /year + 2000, (thisWorstAccel_ThrustData_Trans[1:,1])/thisWorstAccel_SpacecraftMass_Trans, c="C4", linestyle=transLinestyle)#, s=scatterDotSize)                                   # Shows worst thrust data over time


        # Do hlines for mean values
            plt.axhline(thisBestAccelValue_Bare, c="C0")
            plt.axhline(thisWorstAccelValue_Bare, c="C1")
            plt.axhline(meanAccelBase, c="C2")

            plt.axhline(thisBestAccelValue_Trans, c="C3", linestyle=transLinestyle)
            plt.axhline(thisWorstAccelValue_Trans, c="C4", linestyle=transLinestyle)

            # Finalise the plot with legend and saving
            plt.legend(thisAccelLegend)

            plt.savefig(os.path.join(configPlotFolder, "time-accel-%s.png" %thisParameterType))


        ##############################################################################################################
        ############################# DO PLOTTING FOR PARAMETER AGAINST THRUST / ACCEL ###############################
        ##############################################################################################################

        if plottingVsParameterValue:


            # Set plotting units for some itchy parameters
            parameterScaler = 1
            parameterUnits = typesUnits[i]
            startIndex = 0
            if thisParameterType is "diameters":
                parameterScaler = 1E3
                parameterUnits = "mm"
                startIndex = 27

            elif thisParameterType is "lengths":
                parameterScaler = 1E-3
                parameterUnits = "km"

            configSensitivityRunnerValuesToPlotBare_THISTIME = parameterScaler * np.array(configSensitivityRunnerValuesToPlotBare[i])
            configSensitivityRunnerValuesToPlotTrans_THISTIME = parameterScaler * np.array(configSensitivityRunnerValuesToPlotTrans[i])

            ############################### Do Thrust Parameter Plots ####################################################

            plt.figure(figsize=thisFigureSize)
            ax = plt.gca()
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            ax.ticklabel_format(useOffset=False, style="plain")
            plt.xlabel("%s [%s]" %(typesTodoPlotTitles[i], parameterUnits))
            plt.ylabel("Thrust [nN]")
            if thisParameterType in logPlotTypesThrust:
                ax.set_yscale("log")
                ax.set_xscale("log")

            plt.grid("both")




            plt.plot(configSensitivityRunnerValuesToPlotBare_THISTIME[startIndex:],  1E9* np.array(listOfMeanThrustsBare[i][startIndex:]), c="C0")
            plt.plot(configSensitivityRunnerValuesToPlotTrans_THISTIME, 1E9* np.array(listOfMeanThrustsTrans[i]), c="C1")

            idx = utils.findNearestInArray(configSensitivityRunnerValuesToPlotBare[i], baseTypeValues[i])[1]
            plt.plot(configSensitivityRunnerValuesToPlotBare_THISTIME[idx], 1E9*listOfMeanThrustsBare[i][idx], 'o', markersize=20, c="C2")
            # plt.axvline(parameterScaler * baseTypeValues[i], c="C2")

            plt.legend(["Bare", "Transient", "Base"])



            plt.savefig(os.path.join(configPlotFolder, "par-thrust-%s.png" %thisParameterType))



            ############################### Do Accel Parameter Plots ####################################################

            plt.figure(figsize=thisFigureSize)
            ax2 = plt.gca()
            ax2.get_yaxis().get_major_formatter().set_useOffset(False)
            ax2.ticklabel_format(useOffset=False, style="plain")
            plt.xlabel("%s [%s]" %(typesTodoPlotTitles[i], parameterUnits))
            plt.ylabel("Acceleration [pm/s$^2$]")
            if thisParameterType in logPlotTypesAccel:
                plt.xscale("log")
                if thisParameterType is not "lengths": plt.yscale("log")

            plt.grid(which="both")

            # plt.title()



            plt.plot(configSensitivityRunnerValuesToPlotBare_THISTIME[startIndex:],  np.array(1E12*np.array(listOfMeanAccelsBare[i][startIndex:])), c="C0")
            plt.plot(configSensitivityRunnerValuesToPlotTrans_THISTIME, np.array(1E12*np.array(listOfMeanAccelsTrans[i])), c="C1")

            plt.plot(configSensitivityRunnerValuesToPlotBare_THISTIME[idx], 1E12*listOfMeanAccelsBare[i][idx], 'o', markersize=20, c="C2")
            # plt.axvline(parameterScaler * baseTypeValues[i], c="C2")

            plt.legend(["Bare", "Transient", "Base"])



            plt.savefig(os.path.join(configPlotFolder, "par-accel-%s.png" %thisParameterType))





    #plot trajectory of base case
    utils.plotTrajectoryData(propDataBase, plotType="x-y", saveFolder=configPlotFolder, savename="Base trajectory", sameScale=True)


if showing: plt.show()