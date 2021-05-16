import os
import numpy as np
import json
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import copy
import sys
import glob
import shutil


# # Set logfile
logfileDir = utils.pythonRunnerLogfilesDir
utils.checkFolderExist(logfileDir)
logfilePath = os.path.join(logfileDir, utils.createLogfileName("ConfigSensitivityRunner"))
utils.logger = utils.createLogger(logfilePath)

########################### General info #############################################

runSimulations = True
deleteOutputs = True


simulationRunPath = os.path.join(utils.cppApplications_dir, "application_simulation_SSO-CHB")
SSOJsonSubDir = os.path.join("finalSims", "SSO_Config_Sensitivity")

# Set json directory and create if not existing. Also set the base json
initJsonName = "SSO_Bare_Base.json"
initJson = os.path.join(utils.jsonInputs_dir, SSOJsonSubDir, initJsonName)
jsonSaveDir = os.path.join(utils.jsonInputs_dir, SSOJsonSubDir)
utils.checkFolderExist(jsonSaveDir)

# Delete all old jsons:
fileList = glob.glob(jsonSaveDir + "/*")
for filePath in fileList:
    if initJsonName not in filePath:
        os.remove(filePath)

# Create the base json for transition type configuration ie CHTr, using the bare one as a template
initJsonTransName = "SSO_Trans_Base.json"
initJsonTrans = os.path.join(utils.jsonInputs_dir, SSOJsonSubDir, initJsonTransName)
initJsonTransChangeKeys = [ ["EDTConfigs", "configType"],
                            ["saveDataConfigs", "outputSubFolder"],
                            ["saveDataConfigs", "baseFilename"]]
initJsonTransChangeValues = ["CHTr", "SSO-Configs-Sensitivity/SSO-Trans-Base/", "SSO-Trans-Base-"]
utils.createModifiedJson(initJson, jsonSaveDir, initJsonTransName, initJsonTransChangeKeys, initJsonTransChangeValues)

###### Create data for new jsons to use #######

configTypesToTest = ["CHB", "CHTr"]

# SEE UTILS FOR THE CREATION DATA
lengths_m, diameters_m, currents_mA, areaRatios, noLinesRange, lengthRatios, endmassMasses, rotationCoefficients = utils.configSensitivityRunnerValues

# Specifies changelist and changenames
changesList = [lengths_m, diameters_m, currents_mA, areaRatios,   noLinesRange, lengthRatios,   endmassMasses,   rotationCoefficients]
changeNames = ["lengths", "diameters", "currents",  "areaRatios", "noLines",    "lengthRatios", "endmassMasses", "rotationCoefficients"]

## Specifies base values for json and output filenames and directories
SSOConfigsSavenameBaseTop = "SSO_%s_%s-%s.json"

SSOSubdir = "SSO-Configs-Sensitivity/"
SSOSavenameBaseTop = "SSO_%s-%s-%s"



if deleteOutputs:
    # Deletes all folders and files in the output directory too
    SSOOutputPath = os.path.join(utils.simulation_output_dir, SSOSubdir)
    dirListOutput = glob.glob(SSOOutputPath + "/*")
    for filePath in dirListOutput:
        shutil.rmtree(filePath)


for k in range(len(configTypesToTest)):

    thisConfigType = configTypesToTest[k]
    if thisConfigType == "CHB":
        SSOConfigsSavenameBase = SSOConfigsSavenameBaseTop %("Bare", "%s", "%s")
        SSOSavenameBase = SSOSavenameBaseTop %("Bare", "%s", "%s")
        initJsonToUse = initJson

    elif thisConfigType == "CHTr":
        SSOConfigsSavenameBase = SSOConfigsSavenameBaseTop %("Trans", "%s", "%s")
        SSOSavenameBase = SSOSavenameBaseTop %("Trans", "%s", "%s")
        initJsonToUse = initJsonTrans




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

            elif changeName == "endmassMasses":
                SSOChangeKeysTemp = [ ["scConfigs", "SRP", "endmassMass1"],
                                      ["scConfigs", "SRP", "endmassMass2"],
                                      ["scConfigs", "endmassMass1"],
                                      ["scConfigs", "endmassMass2"],
                                      ["saveDataConfigs", "outputSubFolder"],
                                      ["saveDataConfigs", "baseFilename"]]

            elif changeName == "rotationCoefficients":
                SSOChangeKeysTemp = [ ["EDTConfigs", "generalRotationCoefficient"],
                                      ["saveDataConfigs", "outputSubFolder"],
                                      ["saveDataConfigs", "baseFilename"]]

            else:
                utils.logger.info("Change name %s not recognised" %changeName)

            # Specifies the output folder and filename for this run case, to use in changeValues, and also the json

            SSOJsonSavenameTemp = SSOConfigsSavenameBase %(changeName, nameIndex)

            thisOutputSubFolder = SSOSubdir + SSOSavenameBase %(changeName, nameIndex) + "/"
            thisBaseFilename = SSOSavenameBase %(changeName, nameIndex)



            # # Creates folder if not existing, and also empties it so that all values are new. Output only
            # utils.checkFolderExist(os.path.join(utils.jsonInputs_dir))

            if (changeName == "diameters") or (changeName == "areaRatios"):
                SSOChangeValuesTemp = [currentChangeList[j],
                                       currentChangeList[j],
                                       thisOutputSubFolder,
                                       thisBaseFilename]

            elif changeName == "endmassMasses":
                SSOChangeValuesTemp = [currentChangeList[j],
                                       currentChangeList[j],
                                       currentChangeList[j],
                                       currentChangeList[j],
                                       thisOutputSubFolder,
                                       thisBaseFilename]

            else:
                SSOChangeValuesTemp = [currentChangeList[j],
                                   thisOutputSubFolder,
                                   thisBaseFilename]




            utils.createModifiedJson(initJsonToUse, jsonSaveDir, SSOJsonSavenameTemp, SSOChangeKeysTemp, SSOChangeValuesTemp)


if runSimulations:
    utils.runAllSimulations(SSOJsonSubDir, printSetting=2, runPath=simulationRunPath, runOnlyThisFile=None, printProgress=True)

sys.stdout.close()





