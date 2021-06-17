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
logfilePath = os.path.join(logfileDir, utils.createLogfileName("ConfigTrueSensitivityRunner"))
utils.logger = utils.createLogger(logfilePath)

simulationRunPath = os.path.join(utils.cppApplications_dir, "application_simulation_SSO-CHB")
utils.checkFolderExist(utils.jsonSaveDir_ConfigSensitivity, emptyDirectory=False)

# Some common values
spacingValues = utils.configTrueSensitivitySpacingValues





############################## Create jsons #####################################################

# Primary line separation ratio (ka)
utils.createConfigSensitivityJsons(utils.changeKeys_PrimaryLineSeparationRatio, utils.sensitivityRange_PrimaryLineSeparationRatio, spacingValues, utils.outputSubFolder_PrimaryLineSeparationRatio, utils.baseFilenameBase_PrimaryLineSeparationRatio, utils.jsonNameBase_PrimaryLineSeparationRatio)

# Secondary tether diameter
utils.createConfigSensitivityJsons(utils.changeKeys_SecondaryTetherDiameter, utils.sensitivityRange_SecondaryTetherDiameter, spacingValues, utils.outputSubFolder_SecondaryTetherDiameter, utils.baseFilenameBase_SecondaryTetherDiameter, utils.jsonNameBase_SecondaryTetherDiameter)

# Secondary tether area ratio
utils.createConfigSensitivityJsons(utils.changeKeys_SecondaryTetherAreaRatio, utils.sensitivityRange_SecondaryTetherAreaRatio, spacingValues, utils.outputSubFolder_SecondaryTetherAreaRatio, utils.baseFilenameBase_SecondaryTetherAreaRatio, utils.jsonNameBase_SecondaryTetherAreaRatio)

# Slack coefficient
utils.createConfigSensitivityJsons(utils.changeKeys_SlackCoefficient, utils.sensitivityRange_SlackCoefficient, spacingValues, utils.outputSubFolder_SlackCoefficient, utils.baseFilenameBase_SlackCoefficient, utils.jsonNameBase_SlackCoefficient)

# # Occultation coefficient
# utils.createConfigSensitivityJsons(utils.changeKeys_OccultationCoefficient, utils.sensitivityRange_OccultationCoefficient, spacingValues, utils.outputSubFolder_OccultationCoefficient, utils.baseFilenameBase_OccultationCoefficient, utils.jsonNameBase_OccultationCoefficient)

# Al resistivity
utils.createConfigSensitivityJsons(utils.changeKeys_AlResistivity, utils.sensitivityRange_AlResistivity, spacingValues, utils.outputSubFolder_AlResistivity, utils.baseFilenameBase_AlResistivity, utils.jsonNameBase_AlResistivity)

# Cu Resistivity
utils.createConfigSensitivityJsons(utils.changeKeys_CuResistivity, utils.sensitivityRange_CuResistivity, spacingValues, utils.outputSubFolder_CuResistivity, utils.baseFilenameBase_CuResistivity, utils.jsonNameBase_CuResistivity)

# Al density
utils.createConfigSensitivityJsons(utils.changeKeys_AlDensity, utils.sensitivityRange_AlDensity, spacingValues, utils.outputSubFolder_AlDensity, utils.baseFilenameBase_AlDensity, utils.jsonNameBase_AlDensity)

# Cu density
utils.createConfigSensitivityJsons(utils.changeKeys_CuDensity, utils.sensitivityRange_CuDensity, spacingValues, utils.outputSubFolder_CuDensity, utils.baseFilenameBase_CuDensity, utils.jsonNameBase_CuDensity)



########################################## Run Jsons '###############################################################################################################

runSimulations = True

if runSimulations:
    utils.runAllSimulations(utils.configSensitivitySubDir, printSetting=0, runPath=simulationRunPath, runOnlyThisFile=None, printProgress=True)

sys.stdout.close()

################################################# Assess data and do plotting ##############################################################################



