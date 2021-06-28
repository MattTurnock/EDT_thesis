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
fileStringIgnores = ["Parker", "Current", "SRP", "BACKUP" ]
# fileStringIgnores = ["BACKUP" ]
# fileStringIgnores = ["SRP", "Current", "Integrator", "BACKUP"]

#NOTE: jsons create manually for most of these, not by python function
############################## Create relevant jsons to run, for the current vnv stuff ####################

allBaseVariables = [utils.baseVariables_1a, utils.baseVariables_1b, utils.baseVariables_2a]
utils.createCurrentVNVJsons(allBaseVariables)


############################################ Run simulations ########################################################

runSims = True
printSetting=1
if runSims:
    utils.runAllSimulations(VNVJsonSubDir, printSetting=printSetting, runPath=simulationRunPath, fileIgnores=fileStringIgnores)

sys.stdout.close()