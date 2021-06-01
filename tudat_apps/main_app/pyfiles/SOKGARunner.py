import os
import numpy as np
import json
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import copy
import sys



# # Set logfile
logfileDir = utils.pythonRunnerLogfilesDir
utils.checkFolderExist(logfileDir)
logfilePath = os.path.join(logfileDir, utils.createLogfileName("SOKGARunner"))
utils.logger = utils.createLogger(logfilePath)

########################### General info #############################################

runStage1s = True
runStage2s = True
printSetting=1



simulationRunPath = os.path.join(utils.cppApplications_dir, "application_simulation_SSO-CHB")

SOKGAStage1JsonSaveDirPath = os.path.join(utils.jsonInputs_dir, utils.SOKGAStage1JsonSubDir)
SOKGAStage2JsonSaveDirPath = os.path.join(utils.jsonInputs_dir, utils.SOKGAStage2JsonSubDir)

# Delete all old jsons:
fileList = utils.glob.glob(SOKGAStage1JsonSaveDirPath + "/*")
for filePath in fileList:
    os.remove(filePath)

fileList = utils.glob.glob(SOKGAStage2JsonSaveDirPath + "/*")
for filePath in fileList:
    os.remove(filePath)




SOKGABaseJson = os.path.join(utils.jsonInputs_dir, "finalSims", "SOKGA_Base.json")


############################################## Create Stage 1 simulation jsons ######################################################

print("Creating Stage 1 Jsons")

processedDataJupiter = np.load(os.path.join(utils.numpyBinary_dir, "GAStage2NewArrayData_Jupiter_Processed.npy"), allow_pickle=True)
processedDataSaturn = np.load(os.path.join(utils.numpyBinary_dir, "GAStage2NewArrayData_Saturn_Processed.npy"), allow_pickle=True)

SOKGAStage1JsonSavenameBase = "SOKGA_Stage1_%s_%s.json"
SOKGAStage1OutputSubFolderBase = "SOKGA/SOKGA_Stage1_%s_%s/"
SOKGAStage1FilenameBase = "SOGKA_Stage1_%s_%s-"

SOKGAStage2JsonSavenameBase = "SOKGA_Stage2_%s_%s.json"
SOKGAStage2OutputSubFolderBase = "SOKGA/SOKGA_Stage2_%s_%s/"
SOKGAStage2FilenameBase = "SOGKA_Stage2_%s_%s-"

SOKGAChangeKeys = [ ["Spice", "bodiesToInclude", "Mercury"],
                    ["Spice", "bodiesToInclude", "Venus"],
                    ["Spice", "bodiesToInclude", "Earth"],
                    ["Spice", "bodiesToInclude", "Mars"],
                    ["Spice", "bodiesToInclude", "Jupiter"],
                    ["Spice", "bodiesToInclude", "Saturn"],
                    ["Spice", "bodiesToInclude", "Uranus"],
                    ["Spice", "bodiesToInclude", "Neptune"],
                    ["GuidanceConfigs", "thrustMagnitudeConfig"],
                    ["GuidanceConfigs", "initialEphemerisYear"],
                    ["GuidanceConfigs", "terminationSettings", "terminationType"],
                    ["GuidanceConfigs", "terminationSettings", "absoluteTimeTerminationYear"],
                    ["GuidanceConfigs", "terminationSettings", "proximityTerminationBody2"],
                    ["GuidanceConfigs", "terminationSettings", "proximityTerminationCutoffAU"],
                    ["GuidanceConfigs", "vehicleInitialCartesian", "x1_m"],
                    ["GuidanceConfigs", "vehicleInitialCartesian", "x2_m"],
                    ["GuidanceConfigs", "vehicleInitialCartesian", "x3_m"],
                    ["GuidanceConfigs", "vehicleInitialCartesian", "v1_ms"],
                    ["GuidanceConfigs", "vehicleInitialCartesian", "v2_ms"],
                    ["GuidanceConfigs", "vehicleInitialCartesian", "v3_ms"],
                    ["saveDataConfigs", "outputSubFolder"],
                    ["saveDataConfigs", "baseFilename"],
                    ["scConfigs", "useSRP"]]

for i in range(len(processedDataJupiter)):
    thisJsonSavename = SOKGAStage1JsonSavenameBase %("Jupiter", i)
    thisOutputSubFolder = SOKGAStage1OutputSubFolderBase %("Jupiter", i)
    thisFilename = SOKGAStage1FilenameBase %("Jupiter", i)

    initialEphemerisYear = 2000 + processedDataJupiter[i, 1]

    changeValues = [0,
                    0,
                    0,
                    0,
                    1,
                    0,
                    0,
                    0,
                    "disabled",
                    initialEphemerisYear,
                    "proximityTermination",
                    utils.SOKGASimulationEndYear,
                    "Jupiter",
                    utils.closeApproachCutoffAU,
                    processedDataJupiter[i, 2],
                    processedDataJupiter[i, 3],
                    processedDataJupiter[i, 4],
                    processedDataJupiter[i, 5],
                    processedDataJupiter[i, 6],
                    processedDataJupiter[i, 7],
                    thisOutputSubFolder,
                    thisFilename,
                    False]

    utils.createModifiedJson(SOKGABaseJson, SOKGAStage1JsonSaveDirPath, thisJsonSavename, SOKGAChangeKeys, changeValues)



for i in range(len(processedDataSaturn)):
    thisJsonSavename = SOKGAStage1JsonSavenameBase %("Saturn", i)
    thisOutputSubFolder = SOKGAStage1OutputSubFolderBase %("Saturn", i)
    thisFilename = SOKGAStage1FilenameBase %("Saturn", i)

    initialEphemerisYear = 2000 + processedDataSaturn[i, 1]

    changeValues = [0,
                    0,
                    0,
                    0,
                    0,
                    1,
                    0,
                    0,
                    "disabled",
                    initialEphemerisYear,
                    "proximityTermination",
                    utils.SOKGASimulationEndYear,
                    "Saturn",
                    utils.closeApproachCutoffAU,
                    processedDataSaturn[i, 2],
                    processedDataSaturn[i, 3],
                    processedDataSaturn[i, 4],
                    processedDataSaturn[i, 5],
                    processedDataSaturn[i, 6],
                    processedDataSaturn[i, 7],
                    thisOutputSubFolder,
                    thisFilename,
                    False]

    utils.createModifiedJson(SOKGABaseJson, SOKGAStage1JsonSaveDirPath, thisJsonSavename, SOKGAChangeKeys, changeValues)


############################################## Run Stage 1 simulations ######################################################

if runStage1s:
    print("Running Stage 1 Simulations")
    utils.runAllSimulations(utils.SOKGAStage1JsonSubDir, printSetting=printSetting, runPath=simulationRunPath, printProgress=True)


############################################## Create Stage 2 simulation jsons ######################################################
print("Creating Stage 2 Jsons")


allStage1Jsons = utils.natsort.natsorted(utils.glob.glob(SOKGAStage1JsonSaveDirPath + "/*"))

for i in range(len(allStage1Jsons)):
    thisStage1Json = allStage1Jsons[i]
    allSimData = utils.getAllSimDataFromJson(thisStage1Json, printInfo=True, todoList=["propData"])
    propDataArray = allSimData[5]

    if "Jupiter" in thisStage1Json:
        planet = "Jupiter"
        index = i
    elif "Saturn" in thisStage1Json:
        planet = "Saturn"
        index = i - len(processedDataJupiter)

    thisJsonSavename = SOKGAStage2JsonSavenameBase %(planet, index)
    thisOutputSubFolder = SOKGAStage2OutputSubFolderBase %(planet, index)
    thisFilename = SOKGAStage2FilenameBase %(planet, index)

    initialEphemerisYear = 2000 + propDataArray[-1, 0] / utils.year

    changeValues = [1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    "nominal",
                    initialEphemerisYear,
                    "nominalTimeTermination",
                    utils.SOKGASimulationEndYear,
                    planet,
                    utils.closeApproachCutoffAU,
                    propDataArray[-1, 1],
                    propDataArray[-1, 2],
                    propDataArray[-1, 3],
                    propDataArray[-1, 4],
                    propDataArray[-1, 5],
                    propDataArray[-1, 6],
                    thisOutputSubFolder,
                    thisFilename,
                    True]

    utils.createModifiedJson(SOKGABaseJson, SOKGAStage2JsonSaveDirPath, thisJsonSavename, SOKGAChangeKeys, changeValues)


############################################## Run Stage 1 simulations ######################################################

if runStage2s:
    print("Running Stage 2 Simulations")
    utils.runAllSimulations(utils.SOKGAStage2JsonSubDir, printSetting=printSetting, runPath=simulationRunPath, printProgress=True)




# #### Create second stage SOKGA json, after running first stage init #######
#
# # Get SOKGA data
# dataSubDir_SOKGA = "SOKGA"
# allSimData_SOKGA = utils.getAllSimDataFromFolder(dataSubDir_SOKGA)
# propData_SOKGA = allSimData_SOKGA[5]
# bodyData_SOKGA = allSimData_SOKGA[0]
#
# # Create new json
#
# SOKGAStage2JsonSavename = "SOKGA_stage2.json"
# SOKGAChangeKeys = [ ["GuidanceConfigs", "thrustMagnitudeConfig"],
#                     ["GuidanceConfigs", "thrustDirectionConfig"],
#                     ["GuidanceConfigs", "initialEphemerisYear"],
#                     ["GuidanceConfigs", "terminationSettings", "terminationType"],
#                     ["GuidanceConfigs", "terminationSettings", "timeTerminationYears"],
#                     ["GuidanceConfigs", "terminationSettings", "absoluteTimeTerminationYear"],
#                     ["GuidanceConfigs", "initialStateType"],
#                     ["GuidanceConfigs", "vehicleInitialCartesian", "x1_m"],
#                     ["GuidanceConfigs", "vehicleInitialCartesian", "x2_m"],
#                     ["GuidanceConfigs", "vehicleInitialCartesian", "x3_m"],
#                     ["GuidanceConfigs", "vehicleInitialCartesian", "v1_ms"],
#                     ["GuidanceConfigs", "vehicleInitialCartesian", "v2_ms"],
#                     ["GuidanceConfigs", "vehicleInitialCartesian", "v3_ms"],
#                     ["EDTConfigs", "emitterCurrentmA"],
#                     ["saveDataConfigs", "outputSubFolder"],
#                     ["saveDataConfigs", "baseFilename"]]
#
# SOKGAChangeValues = [ "nominal",
#                       "nominalPrograde",
#                       2000 + propData_SOKGA[-1, 0]/utils.year,
#                       "nominalTimeTermination",
#                       999999,
#                       2500,
#                       "Cartesian",
#                       propData_SOKGA[-1, 1],
#                       propData_SOKGA[-1, 2],
#                       propData_SOKGA[-1, 3],
#                       propData_SOKGA[-1, 4],
#                       propData_SOKGA[-1, 5],
#                       propData_SOKGA[-1, 6],
#                       1000,
#                       "SOKGA-Stage2/",
#                       "SOKGA-Stage2-"]
#
# utils.createModifiedJson(SOKGABaseJson, jsonSaveDir, SOKGAStage2JsonSavename, SOKGAChangeKeys, SOKGAChangeValues)
#
# #### Run second stage SOKGA powered trajectory #######
#
# utils.runAllSimulations(SOKGAJsonSubDir, printSetting=printSetting, runPath=simulationRunPath, runOnlyThisFile="SOKGA_stage2.json")
#
# sys.stdout.close()
#
#
# #### Create and run SOKGA reference trajectory, using initial trajectory but with nominalTimeTermination #####
#
# SOKGAReferenceJsonSavename = "SOKGA_reference.json"
# SOKGAReferenceChangeKeys = [ ["GuidanceConfigs", "terminationSettings", "terminationType"],
#                              ["GuidanceConfigs", "terminationSettings", "timeTerminationYears"],
#                              ["GuidanceConfigs", "terminationSettings", "absoluteTimeTerminationYear"],
#                              ["EDTConfigs", "emitterCurrentmA"],
#                              ["saveDataConfigs", "outputSubFolder"],
#                              ["saveDataConfigs", "baseFilename"]]
#
# SOKGAReferenceChangeValues = [ "nominalTimeTermination",
#                                999999,
#                                2500,
#                                0,
#                                "SOKGA-Reference/",
#                                "SOKGA-Reference-"]
#
# utils.createModifiedJson(SOKGABaseJson, jsonSaveDir, SOKGAReferenceJsonSavename, SOKGAReferenceChangeKeys, SOKGAReferenceChangeValues)
#
# #### Run reference SOKGA powered trajectory #######
#
# utils.runAllSimulations(SOKGAJsonSubDir, printSetting=printSetting, runPath=simulationRunPath, runOnlyThisFile="SOKGA_reference.json")
#
# sys.stdout.close()





