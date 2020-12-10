import os
import numpy as np
import json
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import copy
import sys



# # Set logfile
logfileDir = utils.pythonRunnerLogfilesDir
utils.checkFolderExist(logfileDir)
logfilePath = os.path.join(logfileDir, utils.createLogfileName("VNVRunner"))
utils.logger = utils.createLogger(logfilePath)

########################### General info #############################################

simulationRunPath = os.path.join(utils.cppApplications_dir, "application_simulation_SSO-CHB")
VNVJsonSubDir = "VnV"

templateJson = os.path.join(utils.jsonInputs_dir, "testVariables.json")
jsonSaveDir = os.path.join(utils.jsonInputs_dir, "VnV/")
fileStringIgnores = ["Parker"]

#NOTE: jsons create manually for most of these, not by python function
# ############################################ Create json files ########################################################
#
# # Create json for the Voyager-based parker VnV
# voyagerJsonSavename = "testVariablesParkerVnV.json"
# voyagerChangeKeys = [["Spice", "bodiesToInclude", "Sun"],
#                      ["Spice", "bodiesToInclude", "Mercury"],
#                      ["Spice", "bodiesToInclude", "Venus"],
#                      ["Spice", "bodiesToInclude", "Earth"],
#                      ["Spice", "bodiesToInclude", "Mars"],
#                      ["Spice", "bodiesToInclude", "Jupiter"],
#                      ["Spice", "bodiesToInclude", "Saturn"],
#                      ["Spice", "bodiesToInclude", "Uranus"],
#                      ["Spice", "bodiesToInclude", "Neptune"],
#
#                      ["ParkerMagField", "twoSinePars", "a1"],
#                      ["ParkerMagField", "twoSinePars", "b1"],
#                      ["ParkerMagField", "twoSinePars", "c1"],
#                      ["ParkerMagField", "twoSinePars", "a2"],
#                      ["ParkerMagField", "twoSinePars", "b2"],
#                      ["ParkerMagField", "twoSinePars", "c2"],
#                      ["ParkerMagField", "twoSinePars", "d"],
#
#                      ["GuidanceConfigs", "initialEphemerisYear"],
#                      ["GuidanceConfigs", "terminationSettings", "timeTerminationYears"],
#                      ["GuidanceConfigs", "terminationSettings", "absoluteTimeTerminationYears"],]


############################################ Run simulations ########################################################

runSims = True
printSetting=1
if runSims:
    utils.runAllSimulations(VNVJsonSubDir, printSetting=printSetting, runPath=simulationRunPath, fileIgnores=fileStringIgnores)

sys.stdout.close()