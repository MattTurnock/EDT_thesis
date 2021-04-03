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

simulationRunPath = os.path.join(utils.cppApplications_dir, "application_simulation_SSO-CHB")
SOKGAJsonSubDir = "finalSims"


initJson = os.path.join(utils.jsonInputs_dir, SOKGAJsonSubDir, "SOKGA_init.json")
jsonSaveDir = os.path.join(utils.jsonInputs_dir, "finalSims/")


#### Run intial SOKGA unpowered trajectory #######

printSetting=1

utils.runAllSimulations(SOKGAJsonSubDir, printSetting=printSetting, runPath=simulationRunPath, runOnlyThisFile="SOKGA_init.json")

sys.stdout.close()


#### Create second stage SOKGA json, after running first stage init #######

# Get SOKGA data
dataSubDir_SOKGA = "SOKGA"
allSimData_SOKGA = utils.getAllSimDataFromFolder(dataSubDir_SOKGA)
propData_SOKGA = allSimData_SOKGA[5]
bodyData_SOKGA = allSimData_SOKGA[0]

# Create new json

SOKGAStage2JsonSavename = "SOKGA_stage2.json"
SOKGAChangeKeys = [ ["GuidanceConfigs", "thrustMagnitudeConfig"],
                    ["GuidanceConfigs", "thrustDirectionConfig"],
                    ["GuidanceConfigs", "initialEphemerisYear"],
                    ["GuidanceConfigs", "terminationSettings", "terminationType"],
                    ["GuidanceConfigs", "terminationSettings", "timeTerminationYears"],
                    ["GuidanceConfigs", "terminationSettings", "absoluteTimeTerminationYear"],
                    ["GuidanceConfigs", "initialStateType"],
                    ["GuidanceConfigs", "vehicleInitialCartesian", "x1_m"],
                    ["GuidanceConfigs", "vehicleInitialCartesian", "x2_m"],
                    ["GuidanceConfigs", "vehicleInitialCartesian", "x3_m"],
                    ["GuidanceConfigs", "vehicleInitialCartesian", "v1_ms"],
                    ["GuidanceConfigs", "vehicleInitialCartesian", "v2_ms"],
                    ["GuidanceConfigs", "vehicleInitialCartesian", "v3_ms"],
                    ["EDTConfigs", "emitterCurrentmA"],
                    ["saveDataConfigs", "outputSubFolder"],
                    ["saveDataConfigs", "baseFilename"]]

SOKGAChangeValues = [ "nominal",
                      "nominalPrograde",
                      2000 + propData_SOKGA[-1, 0]/utils.year,
                      "nominalTimeTermination",
                      999999,
                      2500,
                      "Cartesian",
                      propData_SOKGA[-1, 1],
                      propData_SOKGA[-1, 2],
                      propData_SOKGA[-1, 3],
                      propData_SOKGA[-1, 4],
                      propData_SOKGA[-1, 5],
                      propData_SOKGA[-1, 6],
                      1000,
                      "SOKGA-Stage2/",
                      "SOKGA-Stage2-"]

utils.createModifiedJson(initJson, jsonSaveDir, SOKGAStage2JsonSavename, SOKGAChangeKeys, SOKGAChangeValues)

#### Run second stage SOKGA powered trajectory #######

utils.runAllSimulations(SOKGAJsonSubDir, printSetting=printSetting, runPath=simulationRunPath, runOnlyThisFile="SOKGA_stage2.json")

sys.stdout.close()


#### Create and run SOKGA reference trajectory, using initial trajectory but with nominalTimeTermination #####

SOKGAReferenceJsonSavename = "SOKGA_reference.json"
SOKGAReferenceChangeKeys = [ ["GuidanceConfigs", "terminationSettings", "terminationType"],
                             ["GuidanceConfigs", "terminationSettings", "timeTerminationYears"],
                             ["GuidanceConfigs", "terminationSettings", "absoluteTimeTerminationYear"],
                             ["EDTConfigs", "emitterCurrentmA"],
                             ["saveDataConfigs", "outputSubFolder"],
                             ["saveDataConfigs", "baseFilename"]]

SOKGAReferenceChangeValues = [ "nominalTimeTermination",
                               999999,
                               2500,
                               0,
                               "SOKGA-Reference/",
                               "SOKGA-Reference-"]

utils.createModifiedJson(initJson, jsonSaveDir, SOKGAReferenceJsonSavename, SOKGAReferenceChangeKeys, SOKGAReferenceChangeValues)

#### Run reference SOKGA powered trajectory #######

utils.runAllSimulations(SOKGAJsonSubDir, printSetting=printSetting, runPath=simulationRunPath, runOnlyThisFile="SOKGA_reference.json")

sys.stdout.close()





