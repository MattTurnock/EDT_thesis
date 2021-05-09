import os
import numpy as np
import json
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import copy
import sys
import glob



# # Set logfile
logfileDir = utils.pythonRunnerLogfilesDir
utils.checkFolderExist(logfileDir)
logfilePath = os.path.join(logfileDir, utils.createLogfileName("InORunner"))
utils.logger = utils.createLogger(logfilePath)

########################### General info #############################################

runSimulations = True


simulationRunPath = os.path.join(utils.cppApplications_dir, "application_simulation_SSO-CHB")
SSOJsonSubDir = os.path.join("finalSims", "SSO_Config_Sensitivity")


initJson = os.path.join(utils.jsonInputs_dir, SSOJsonSubDir, "SSO_Base.json")
jsonSaveDir = os.path.join(utils.jsonInputs_dir, SSOJsonSubDir)

# Delete all old jsons:
fileList = glob.glob(jsonSaveDir + "/*")
for filePath in fileList:
    if "SSO_Base.json" not in filePath:
        os.remove(filePath)

###### Create data for new jsons to use #######

# SEE UTILS FOR THE CREATION DATA
lengths_m, diameters_m, currents_mA, areaRatios, noLinesRange, lengthRatios, slackCoefficients, lineSeparationRatios, occultationCoefficients, endmassMasses = utils.configSensitivityRunnerValues


changesList = [lengths_m, diameters_m, currents_mA, areaRatios, noLinesRange, lengthRatios, slackCoefficients, lineSeparationRatios, occultationCoefficients, endmassMasses]
changeNames = ["lengths", "diameters", "currents", "areaRatios", "noLines", "lengthRatios", "slackCoefficients", "lineSeparationCoefficients", "occultationCoefficients", "endmassMasses"]

SSOConfigsSavenameBase = "SSO_%s-%s.json"

for i in range(len(changesList)):



    changeName = changeNames[i]
    currentChangeList = changesList[i]

    for j in range(len(currentChangeList)):

        nameIndex = j+1

        if changeName == "lengths":
            SSOChangeKeysTemp = [ ["EDTConfigs", "tetherLength"],
                                  ["saveDataConfigs", "outputSubFolder"],
                                  ["saveDataConfigs", "baseFilename"]]

        elif changeName == "diameters":
            SSOChangeKeysTemp = [ ["EDTConfigs", "tetherDiameter"],
                                  ["EDTConfigs", "hoytether", "tetherDiameterSecondary"],
                                  ["saveDataConfigs", "outputSubFolder"],
                                  ["saveDataConfigs", "baseFilename"]]

        elif changeName == "currents":
            SSOChangeKeysTemp = [ ["EDTConfigs", "emitterCurrentmA"],
                                  ["saveDataConfigs", "outputSubFolder"],
                                  ["saveDataConfigs", "baseFilename"]]

        elif changeName == "areaRatios":
            SSOChangeKeysTemp = [ ["EDTConfigs", "tetherAreaRatio"],
                                  ["EDTConfigs", "hoytether", "tetherAreaRatioSecondary"],
                                  ["saveDataConfigs", "outputSubFolder"],
                                  ["saveDataConfigs", "baseFilename"]]

        elif changeName == "noLines":
            SSOChangeKeysTemp = [ ["EDTConfigs", "hoytether", "noPrimaryLines"],
                                  ["saveDataConfigs", "outputSubFolder"],
                                  ["saveDataConfigs", "baseFilename"]]

        elif changeName == "lengthRatios":
            SSOChangeKeysTemp = [ ["EDTConfigs", "hoytether", "primaryLineSegmentLengthRatio"],
                                  ["saveDataConfigs", "outputSubFolder"],
                                  ["saveDataConfigs", "baseFilename"]]

        elif changeName == "slackCoefficients":
            SSOChangeKeysTemp = [ ["EDTConfigs", "hoytether", "slackCoefficient"],
                                  ["saveDataConfigs", "outputSubFolder"],
                                  ["saveDataConfigs", "baseFilename"]]

        elif changeName == "lineSeparationCoefficients":
            SSOChangeKeysTemp = [ ["EDTConfigs", "hoytether", "primaryLineSeparationRatio"],
                                  ["saveDataConfigs", "outputSubFolder"],
                                  ["saveDataConfigs", "baseFilename"]]

        elif changeName == "occultationCoefficients":
            SSOChangeKeysTemp = [ ["EDTConfigs", "hoytether", "occultationCoefficient"],
                                  ["saveDataConfigs", "outputSubFolder"],
                                  ["saveDataConfigs", "baseFilename"]]

        elif changeName == "endmassMasses":
            SSOChangeKeysTemp = [ ["scConfigs", "SRP", "endmassMass1"],
                                  ["scConfigs", "SRP", "endmassMass2"],
                                  ["scConfigs", "endmassMass1"],
                                  ["scConfigs", "endmassMass2"],
                                  ["saveDataConfigs", "outputSubFolder"],
                                  ["saveDataConfigs", "baseFilename"]]

        else:
            utils.logger.info("Change name %s not recognised" %changeName)



        if (changeName == "diameters") or (changeName == "areaRatios"):
            SSOChangeValuesTemp = [currentChangeList[j],
                                   currentChangeList[j],
                                   "SSO-Configs-Sensitivity/SSO-%s-%s/" %(changeName, nameIndex),
                                   "SSO-%s-%s-" %(changeName, nameIndex)]

        elif changeName == "endmassMasses":
            SSOChangeValuesTemp = [currentChangeList[j],
                                   currentChangeList[j],
                                   currentChangeList[j],
                                   currentChangeList[j],
                                   "SSO-Configs-Sensitivity/SSO-%s-%s/" %(changeName, nameIndex),
                                   "SSO-%s-%s-" %(changeName, nameIndex)]

        else:
            SSOChangeValuesTemp = [currentChangeList[j],
                               "SSO-Configs-Sensitivity/SSO-%s-%s/" %(changeName, nameIndex),
                               "SSO-%s-%s-" %(changeName, nameIndex)]

        SSOJsonSavenameTemp = SSOConfigsSavenameBase %(changeName, nameIndex)

        utils.createModifiedJson(initJson, jsonSaveDir, SSOJsonSavenameTemp, SSOChangeKeysTemp, SSOChangeValuesTemp)


if runSimulations:
    utils.runAllSimulations(SSOJsonSubDir, printSetting=0, runPath=simulationRunPath, runOnlyThisFile=None)

sys.stdout.close()


# # Create new json
#
# InOStage2JsonSavename = "InO_stage2.json"
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
# InOChangeValues = [ "nominal",
#                       "nominalPrograde",
#                       2000 + propData_InO[-1, 0]/utils.year,
#                       "nominalTimeTermination",
#                       99999,
#                       2200,
#                       "Cartesian",
#                       propData_InO[-1, 1],
#                       propData_InO[-1, 2],
#                       propData_InO[-1, 3],
#                       propData_InO[-1, 4],
#                       propData_InO[-1, 5],
#                       propData_InO[-1, 6],
#                       2000000000,
#                       "InO-Stage2/",
#                       "InO-Stage2-"]
#
# utils.createModifiedJson(initJson, jsonSaveDir, InOStage2JsonSavename, SOKGAChangeKeys, InOChangeValues)
#
# #### Run second stage SOKGA powered trajectory #######
#
# utils.runAllSimulations(InOJsonSubDir, printSetting=printSetting, runPath=simulationRunPath, runOnlyThisFile="InO_stage2.json")
#
# sys.stdout.close()
#
#
# # #### Create and run SOKGA reference trajectory, using initial trajectory but with nominalTimeTermination #####
# #
# # SOKGAReferenceJsonSavename = "SOKGA_reference.json"
# # SOKGAReferenceChangeKeys = [ ["GuidanceConfigs", "terminationSettings", "terminationType"],
# #                              ["GuidanceConfigs", "terminationSettings", "timeTerminationYears"],
# #                              ["GuidanceConfigs", "terminationSettings", "absoluteTimeTerminationYear"],
# #                              ["EDTConfigs", "emitterCurrentmA"],
# #                              ["saveDataConfigs", "outputSubFolder"],
# #                              ["saveDataConfigs", "baseFilename"]]
# #
# # SOKGAReferenceChangeValues = [ "nominalTimeTermination",
# #                                999999,
# #                                2100,
# #                                0,
# #                                "SOKGA-Reference/",
# #                                "SOKGA-Reference-"]
# #
# # utils.createModifiedJson(initJson, jsonSaveDir, SOKGAReferenceJsonSavename, SOKGAReferenceChangeKeys, SOKGAReferenceChangeValues)
# #
# # #### Run second stage SOKGA powered trajectory #######
# #
# # utils.runAllSimulations(SOKGAJsonSubDir, printSetting=printSetting, runPath=simulationRunPath, runOnlyThisFile="SOKGA_reference.json")
# #
# # sys.stdout.close()
#




