import os
import numpy as np
import json
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import copy
import sys



# # Set logfile
logfileDir = utils.pythonRunnerLogfilesDir
utils.checkFolderExist(logfileDir)
logfilePath = os.path.join(logfileDir, utils.createLogfileName("InORunner"))
utils.logger = utils.createLogger(logfilePath)

########################### General info #############################################

simulationRunPath = os.path.join(utils.cppApplications_dir, "application_simulation_SSO-CHB")
InOJsonSubDir = "finalSims"


initJson = os.path.join(utils.jsonInputs_dir, InOJsonSubDir, "InO_init.json")
jsonSaveDir = os.path.join(utils.jsonInputs_dir, "finalSims/")


#### Run intial SOKGA unpowered trajectory #######

printSetting=1

utils.runAllSimulations(InOJsonSubDir, printSetting=printSetting, runPath=simulationRunPath, runOnlyThisFile="InO_init.json")

sys.stdout.close()


#### Create second stage SOKGA json, after running first stage init #######

# Get SOKGA data
dataSubDir_InO = "InO"
allSimData_InO = utils.getAllSimDataFromFolder(dataSubDir_InO)
propData_InO = allSimData_InO[5]
bodyData_SOKGA = allSimData_InO[0]

# Create new json

InOStage2JsonSavename = "InO_stage2.json"
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

InOChangeValues = [ "nominal",
                      "nominalPrograde",
                      2000 + propData_InO[-1, 0]/utils.year,
                      "nominalTimeTermination",
                      99999,
                      2200,
                      "Cartesian",
                      propData_InO[-1, 1],
                      propData_InO[-1, 2],
                      propData_InO[-1, 3],
                      propData_InO[-1, 4],
                      propData_InO[-1, 5],
                      propData_InO[-1, 6],
                      2000000000,
                      "InO-Stage2/",
                      "InO-Stage2-"]

utils.createModifiedJson(initJson, jsonSaveDir, InOStage2JsonSavename, SOKGAChangeKeys, InOChangeValues)

#### Run second stage SOKGA powered trajectory #######

utils.runAllSimulations(InOJsonSubDir, printSetting=printSetting, runPath=simulationRunPath, runOnlyThisFile="InO_stage2.json")

sys.stdout.close()


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
#                                2100,
#                                0,
#                                "SOKGA-Reference/",
#                                "SOKGA-Reference-"]
#
# utils.createModifiedJson(initJson, jsonSaveDir, SOKGAReferenceJsonSavename, SOKGAReferenceChangeKeys, SOKGAReferenceChangeValues)
#
# #### Run second stage SOKGA powered trajectory #######
#
# utils.runAllSimulations(SOKGAJsonSubDir, printSetting=printSetting, runPath=simulationRunPath, runOnlyThisFile="SOKGA_reference.json")
#
# sys.stdout.close()





