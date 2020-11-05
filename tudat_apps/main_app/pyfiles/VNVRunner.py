import os
import numpy as np
import json
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import copy

########################### General info #############################################

simulationRunPath = os.path.join(utils.cppApplications_dir, "application_simulation_SSO-CHB")
VNVJsonSubDir = "VnV"


############################################ Run simulations ########################################################

runSims = True
printSetting=2
if runSims:
    utils.runAllSimulations(VNVJsonSubDir, printSetting=printSetting, runPath=simulationRunPath)