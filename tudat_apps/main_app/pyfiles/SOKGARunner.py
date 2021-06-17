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

runStage1s = False
runStage2s = False
runStage2s_noEDT = True
printSetting=0



simulationRunPath = os.path.join(utils.cppApplications_dir, "application_simulation_SSO-CHB")

SOKGAStage1JsonSaveDirPath = os.path.join(utils.jsonInputs_dir, utils.SOKGAStage1JsonSubDir)
SOKGAStage2JsonSaveDirPath = os.path.join(utils.jsonInputs_dir, utils.SOKGAStage2JsonSubDir)
SOKGAStage2JsonSaveDirPath_noEDT = os.path.join(utils.jsonInputs_dir, utils.SOKGAStage2JsonSubDir_NoEDT)

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

SOKGAStage2JsonSavenameBase_NoEDT = "SOKGA_Stage2_NoEDT_%s_%s.json"
SOKGAStage2OutputSubFolderBase_NoEDT = "SOKGA/SOKGA_Stage2_NoEDT_%s_%s/"
SOKGAStage2FilenameBase_NoEDT = "SOGKA_Stage2_NoEDT_%s_%s-"

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
    allSimData = utils.getAllSimDataFromJson(thisStage1Json, printInfo=False, todoList=["propData"])
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


############################################## Run Stage 2 simulations ######################################################

if runStage2s:
    print("Running Stage 2 Simulations")
    utils.runAllSimulations(utils.SOKGAStage2JsonSubDir, printSetting=printSetting, runPath=simulationRunPath, printProgress=True)


############################################## Create Stage 2 (no EDT) simulation jsons ######################################################
print("Creating Stage 2 NoEDT Jsons")


# allStage1Jsons = utils.natsort.natsorted(utils.glob.glob(SOKGAStage1JsonSaveDirPath + "/*"))

for i in range(len(allStage1Jsons)):
    thisStage1Json = allStage1Jsons[i]
    allSimData = utils.getAllSimDataFromJson(thisStage1Json, printInfo=False, todoList=["propData"])
    propDataArray = allSimData[5]

    if "Jupiter" in thisStage1Json:
        planet = "Jupiter"
        index = i
    elif "Saturn" in thisStage1Json:
        planet = "Saturn"
        index = i - len(processedDataJupiter)

    thisJsonSavename = SOKGAStage2JsonSavenameBase_NoEDT %(planet, index)
    thisOutputSubFolder = SOKGAStage2OutputSubFolderBase_NoEDT %(planet, index)
    thisFilename = SOKGAStage2FilenameBase_NoEDT %(planet, index)

    initialEphemerisYear = 2000 + propDataArray[-1, 0] / utils.year

    changeValues = [1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    "disabled",
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

    utils.createModifiedJson(SOKGABaseJson, SOKGAStage2JsonSaveDirPath_noEDT, thisJsonSavename, SOKGAChangeKeys, changeValues)


############################################## Run Stage 2 (no EDT) simulations ######################################################

if runStage2s_noEDT:
    print("Running Stage 2 NoEDT Simulations")
    utils.runAllSimulations(utils.SOKGAStage2JsonSubDir_NoEDT, printSetting=printSetting, runPath=simulationRunPath, printProgress=True)


