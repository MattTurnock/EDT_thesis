import os
import numpy as np
import json
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import copy
import glob

########################### General info #############################################
# Simulation excecution settings
runSims = True
runJupiter = False
runSaturn = False
runMars = False

runSecondStage = True
runStage2Jupiter = True
runStage2Saturn = True
printSetting=0

# Run paths and base filenames
GACalculatorRunPath = os.path.join(utils.cppApplications_dir, "application_GA_calculator")
algorithmConfigsSynodic = ["moead", 1000, 100, False, False, 1000, True, False]
algorithmConfigsGlobal = ["moead", 1000, 100, False, True, 1000, True, False]
jsonFilenameBase = "GAConfigs_%s_%s-%s.json"
templateJsonPath = os.path.join(utils.jsonInputs_dir, "GAConfigsNominal.json")

################### Create the json files for Jupiter  #####################################

print("======= Creating first stage json files ==========")

#Values to put into function
outputSubFolderBaseJupiterSynodic = "GAJupiterSynodic/GACalculatorNominal"
outputSubFolderBaseJupiterGlobal = "GAJupiterGlobal/GACalculatorNominal"
jsonSaveSubDirJupiter = "GAConfigs_Jupiter"
inputStartYearsRangeJupiter = [2020, 2050]

# Delete all old jsons:
fileList = glob.glob(os.path.join(utils.jsonInputs_dir, jsonSaveSubDirJupiter) + "/*")
for filePath in fileList:
    os.remove(filePath)

# Create jsons for the case of running synodic period optimisations
utils.createGARunnerJsons(utils.quickConfigsJupiter, outputSubFolderBaseJupiterSynodic, jsonSaveSubDirJupiter, jsonFilenameBase, inputStartYearsRangeJupiter,
                          utils.JupiterInfoList, templateJsonPath=templateJsonPath, createSynodicJsons=True,
                          algorithmConfigs=algorithmConfigsSynodic)
# Create jsons for the case of running glboal optimisations and grid search
utils.createGARunnerJsons(utils.quickConfigsJupiter, outputSubFolderBaseJupiterGlobal, jsonSaveSubDirJupiter, jsonFilenameBase, inputStartYearsRangeJupiter,
                          utils.JupiterInfoList, templateJsonPath=templateJsonPath, createSynodicJsons=False,
                          algorithmConfigs=algorithmConfigsGlobal)



################### Create the json files for Saturn #####################################

#Values to put into function
outputSubFolderBaseSaturnSynodic = "GASaturnSynodic/GACalculatorNominal"
outputSubFolderBaseSaturnGlobal = "GASaturnGlobal/GACalculatorNominal"
jsonSaveSubDirSaturn = "GAConfigs_Saturn"
inputStartYearsRangeSaturn = inputStartYearsRangeJupiter

# Delete all old jsons:
fileList = glob.glob(os.path.join(utils.jsonInputs_dir, jsonSaveSubDirSaturn) + "/*")
for filePath in fileList:
    os.remove(filePath)

# Create jsons for the case of running synodic period optimisations
utils.createGARunnerJsons(utils.quickConfigsSaturn, outputSubFolderBaseSaturnSynodic, jsonSaveSubDirSaturn, jsonFilenameBase, inputStartYearsRangeSaturn,
                          utils.SaturnInfoList, templateJsonPath=templateJsonPath, createSynodicJsons=True,
                          algorithmConfigs=algorithmConfigsSynodic)
# Create jsons for the case of running glboal optimisations and grid search
utils.createGARunnerJsons(utils.quickConfigsSaturn, outputSubFolderBaseSaturnGlobal, jsonSaveSubDirSaturn, jsonFilenameBase, inputStartYearsRangeSaturn,
                          utils.SaturnInfoList, templateJsonPath=templateJsonPath, createSynodicJsons=False,
                          algorithmConfigs=algorithmConfigsGlobal)



################### Create the json files for Mars #####################################

#Values to put into function
outputSubFolderBaseMarsSynodic = "GAMarsSynodic/GACalculatorNominal"
outputSubFolderBaseMarsGlobal = "GAMarsGlobal/GACalculatorNominal"
jsonSaveSubDirMars = "GAConfigs_Mars"
inputStartYearsRangeMars = [2020, 2025]
algorithmConfigsSynodicMars = ["moead", 1024, 100, False, False, 1000, True, True]
algorithmConfigsGlobalMars = ["moead", 1024, 100, False, True, 1000, True, True]

# Delete all old jsons:
fileList = glob.glob(os.path.join(utils.jsonInputs_dir, jsonSaveSubDirMars) + "/*")
for filePath in fileList:
    os.remove(filePath)


# Create jsons for the case of running synodic period optimisations
utils.createGARunnerJsons(utils.quickConfigsMars, outputSubFolderBaseMarsSynodic, jsonSaveSubDirMars, jsonFilenameBase, inputStartYearsRangeMars,
                          utils.MarsInfoList, templateJsonPath=templateJsonPath, createSynodicJsons=True,
                          algorithmConfigs=algorithmConfigsSynodicMars)
# Create jsons for the case of running glboal optimisations and grid search
utils.createGARunnerJsons(utils.quickConfigsMars, outputSubFolderBaseMarsGlobal, jsonSaveSubDirMars, jsonFilenameBase, inputStartYearsRangeMars,
                          utils.MarsInfoList, templateJsonPath=templateJsonPath, createSynodicJsons=False,
                          algorithmConfigs=algorithmConfigsGlobalMars)




############################################ Run First Stage simulations ########################################################

if runSims:
    print("======= Running First Stage Simulations ==========")

    if runJupiter: utils.runAllSimulations(jsonSaveSubDirJupiter, printSetting=printSetting, printProgress=True)

    if runSaturn: utils.runAllSimulations(jsonSaveSubDirSaturn, printSetting=printSetting, printProgress=True)

    if runMars: utils.runAllSimulations(jsonSaveSubDirMars, printSetting=printSetting, printProgress=True)



if runSecondStage:
    ################### Create the json files for Second stage runs #####################################
    print("======= Creating second stage json files ==========")

    simulationRunPath = os.path.join(utils.cppApplications_dir, "application_simulation_SSO-CHB")

    nominalTestVariablesJsonPath = os.path.join(utils.jsonInputs_dir, "GATestVariablesNominal.json")

    stage2JsonChangeKeys = [ ["Spice", "bodiesToInclude", "Jupiter"],
                             ["Spice", "bodiesToInclude", "Saturn"],
                             ["GuidanceConfigs", "initialEphemerisYear"],
                             ["GuidanceConfigs", "terminationSettings", "timeTerminationYears"],
                             ["GuidanceConfigs", "vehicleInitialCartesian", "x1_m"],
                             ["GuidanceConfigs", "vehicleInitialCartesian", "x2_m"],
                             ["GuidanceConfigs", "vehicleInitialCartesian", "x3_m"],
                             ["GuidanceConfigs", "vehicleInitialCartesian", "v1_ms"],
                             ["GuidanceConfigs", "vehicleInitialCartesian", "v2_ms"],
                             ["GuidanceConfigs", "vehicleInitialCartesian", "v3_ms"],
                             ["saveDataConfigs", "outputSubFolder"],
                             ["saveDataConfigs", "baseFilename"] ]


    ## Jupiter ones ##


    terminationYearsJupiter = 10
    outputSubFolderBaseJupiterStage2 = "GAJupiterStage2/GA_Jupiter_Stage2_syn_%s/GA_Stage2_EJ_%s" #ok
    baseFilenameJupiterStage2 = "GA_Stage2_EJ_syn_%s_%s-" #ok

    jsonSaveSubDirJupiterStage2Base = "GATestVariables_Jupiter_Stage2/GA_Stage2_EJ_Syn_%s" #'ok
    jsonSavenameJupiterStage2Base = "GATestVariables_Jupiter_Stage2_%s_%s.json" #ok

    # # Delete all old jsons:
    # fileList = glob.glob(os.path.join(utils.jsonInputs_dir, jsonSaveDirJupiterStage2) + "/*")
    # for filePath in fileList:
    #     os.remove(filePath)


    # GAJupiterGlobalResultsGenInitalStatesFilePath = os.path.join(utils.GAJupiterGlobalResultsDirPath, "initialState_GA_EJ_%s.dat" %utils.GAStage2GenerationNumberToUse)

    GAJupiterSynodicDirPath = os.path.join(utils.simulation_output_dir, "GAJupiterSynodic")
    GAJupiterSynodicSubDirList = utils.natsort.natsorted(glob.glob(GAJupiterSynodicDirPath + "/*"))

    GAJupiterJsonDirList = []
    bar = utils.IncrementalBar("Jupiter json progress", max=len(GAJupiterSynodicSubDirList), suffix='[%(percent).1f%%. ETA: %(eta)ds]   ' )
    for j in range(len(GAJupiterSynodicSubDirList)):

        # jsonSaveDirJupiterStage2 = os.path.join(utils.jsonInputs_dir, jsonSaveSubDirJupiterStage2)
        jsonSaveSubDir = jsonSaveSubDirJupiterStage2Base %j
        GAJupiterJsonDirList.append(jsonSaveSubDir)
        utils.checkFolderExist(os.path.join(utils.jsonInputs_dir, jsonSaveSubDir))

        # Delete all old jsons:
        fileList = glob.glob(os.path.join(utils.jsonInputs_dir, jsonSaveSubDir) + "/*")
        for filePath in fileList:
            os.remove(filePath)

        # GAJupiterInitialStates = np.genfromtxt(utils.GAJupiterGlobalResultsGenInitalStatesFilePath, delimiter=",")

        subDirPath = GAJupiterSynodicSubDirList[j]
        initialStates = np.genfromtxt(os.path.join(subDirPath, utils.GAJupiterInitialStateFilename), delimiter=",")


        for i in range(len(initialStates)):

            initialState = initialStates[i]
            initialEphemerisYear = initialState[1] + 2000
            col1Value = int(initialState[0])

            thisStateChangeValues = [ 1,
                                      0,
                                      initialEphemerisYear,
                                      terminationYearsJupiter,
                                      initialState[2],
                                      initialState[3],
                                      initialState[4],
                                      initialState[5],
                                      initialState[6],
                                      initialState[7],
                                      outputSubFolderBaseJupiterStage2 %(j, col1Value),
                                      baseFilenameJupiterStage2 %(j, col1Value)]


            jsonSavename = jsonSavenameJupiterStage2Base %(j, col1Value)
            utils.createModifiedJson(nominalTestVariablesJsonPath, os.path.join(utils.jsonInputs_dir, jsonSaveSubDir), jsonSavename, stage2JsonChangeKeys, thisStateChangeValues)

        bar.next()
    bar.finish()



    ## Saturn ones ##


    terminationYearsSaturn = 10
    outputSubFolderBaseSaturnStage2 = "GASaturnStage2/GA_Saturn_Stage2_syn_%s/GA_Stage2_ES_%s" #ok
    baseFilenameSaturnStage2 = "GA_Stage2_ES_syn_%s_%s-" #ok

    jsonSaveSubDirSaturnStage2Base = "GATestVariables_Saturn_Stage2/GA_Stage2_ES_Syn_%s"
    jsonSavenameSaturnStage2Base = "GATestVariables_Saturn_Stage2_%s_%s.json" #ok


    GASaturnSynodicDirPath = os.path.join(utils.simulation_output_dir, "GASaturnSynodic")
    GASaturnSynodicSubDirList = utils.natsort.natsorted(glob.glob(GASaturnSynodicDirPath + "/*"))

    GASaturnJsonDirList = []
    bar = utils.IncrementalBar("Saturn json progress", max=len(GASaturnSynodicSubDirList), suffix='[%(percent).1f%%. ETA: %(eta)ds]   ' )
    for j in range(len(GASaturnSynodicSubDirList)):

        jsonSaveSubDir = jsonSaveSubDirSaturnStage2Base %j
        GASaturnJsonDirList.append(jsonSaveSubDir)
        utils.checkFolderExist(os.path.join(utils.jsonInputs_dir, jsonSaveSubDir))

        # Delete all old jsons:
        fileList = glob.glob(os.path.join(utils.jsonInputs_dir, jsonSaveSubDir) + "/*")
        for filePath in fileList:
            os.remove(filePath)


        subDirPath = GASaturnSynodicSubDirList[j]
        initialStates = np.genfromtxt(os.path.join(subDirPath, utils.GASaturnInitialStateFilename), delimiter=",")


        for i in range(len(initialStates)):

            initialState = initialStates[i]
            initialEphemerisYear = initialState[1] + 2000
            col1Value = int(initialState[0])

            thisStateChangeValues = [ 1,
                                      0,
                                      initialEphemerisYear,
                                      terminationYearsSaturn,
                                      initialState[2],
                                      initialState[3],
                                      initialState[4],
                                      initialState[5],
                                      initialState[6],
                                      initialState[7],
                                      outputSubFolderBaseSaturnStage2 %(j, col1Value),
                                      baseFilenameSaturnStage2 %(j, col1Value)]


            jsonSavename = jsonSavenameSaturnStage2Base %(j, col1Value)
            utils.createModifiedJson(nominalTestVariablesJsonPath, os.path.join(utils.jsonInputs_dir, jsonSaveSubDir), jsonSavename, stage2JsonChangeKeys, thisStateChangeValues)

        bar.next()
    bar.finish()






    ############################################ Run Second Stage simulations ########################################################

    if runSims:
        print("======= Running Second Stage Simulations ==========")

        if runStage2Jupiter:
            for i in range(len(GAJupiterJsonDirList)):
                print("Run: ", i)
                subDirPath = GAJupiterJsonDirList[i]
                utils.runAllSimulations(subDirPath, printSetting=printSetting, printProgress=True, runPath=simulationRunPath)

        if runStage2Saturn:
            for i in range(len(GASaturnJsonDirList)):
                print("Run: ", i)
                subDirPath = GASaturnJsonDirList[i]
                utils.runAllSimulations(subDirPath, printSetting=printSetting, printProgress=True, runPath=simulationRunPath)