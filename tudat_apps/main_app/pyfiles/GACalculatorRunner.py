import os
import numpy as np
import json
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import copy

########################### General info #############################################

GACalculatorRunPath = os.path.join(utils.cppApplications_dir, "application_GA_calculator")
algorithmConfigsSynodic = ["moead", 1000, 100, False, False, 1000, True, False]
algorithmConfigsGlobal = ["moead", 1000, 100, False, True, 1000, True, False]
jsonFilenameBase = "GAConfigs_%s_%s-%s.json"
templateJsonPath = os.path.join(utils.jsonInputs_dir, "GAConfigsNominal.json")

################### Create the json files and run sims for Jupiter  #####################################

#Values to put into function
outputSubFolderBaseJupiterSynodic = "GAJupiterSynodic/GACalculatorNominal"
outputSubFolderBaseJupiterGlobal = "GAJupiterGlobal/GACalculatorNominal"
jsonSaveSubDirJupiter = "GAConfigs_Jupiter"
inputStartYearsRangeJupiter = [2020, 2050]

# Create jsons for the case of running synodic period optimisations
utils.createGARunnerJsons(utils.quickConfigsJupiter, outputSubFolderBaseJupiterSynodic, jsonSaveSubDirJupiter, jsonFilenameBase, inputStartYearsRangeJupiter,
                          utils.JupiterInfoList, templateJsonPath=templateJsonPath, createSynodicJsons=True,
                          algorithmConfigs=algorithmConfigsSynodic)
# Create jsons for the case of running glboal optimisations and grid search
utils.createGARunnerJsons(utils.quickConfigsJupiter, outputSubFolderBaseJupiterGlobal, jsonSaveSubDirJupiter, jsonFilenameBase, inputStartYearsRangeJupiter,
                          utils.JupiterInfoList, templateJsonPath=templateJsonPath, createSynodicJsons=False,
                          algorithmConfigs=algorithmConfigsGlobal)



################### Create the json files and run sims for Saturn #####################################

#Values to put into function
outputSubFolderBaseSaturnSynodic = "GASaturnSynodic/GACalculatorNominal"
outputSubFolderBaseSaturnGlobal = "GASaturnGlobal/GACalculatorNominal"
jsonSaveSubDirSaturn = "GAConfigs_Saturn"
inputStartYearsRangeSaturn = inputStartYearsRangeJupiter

# Create jsons for the case of running synodic period optimisations
utils.createGARunnerJsons(utils.quickConfigsSaturn, outputSubFolderBaseSaturnSynodic, jsonSaveSubDirSaturn, jsonFilenameBase, inputStartYearsRangeSaturn,
                          utils.SaturnInfoList, templateJsonPath=templateJsonPath, createSynodicJsons=True,
                          algorithmConfigs=algorithmConfigsSynodic)
# Create jsons for the case of running glboal optimisations and grid search
utils.createGARunnerJsons(utils.quickConfigsSaturn, outputSubFolderBaseSaturnGlobal, jsonSaveSubDirSaturn, jsonFilenameBase, inputStartYearsRangeSaturn,
                          utils.SaturnInfoList, templateJsonPath=templateJsonPath, createSynodicJsons=False,
                          algorithmConfigs=algorithmConfigsGlobal)



################### Create the json files and run sims for Mars #####################################

#Values to put into function
outputSubFolderBaseMarsSynodic = "GAMarsSynodic/GACalculatorNominal"
outputSubFolderBaseMarsGlobal = "GAMarsGlobal/GACalculatorNominal"
jsonSaveSubDirMars = "GAConfigs_Mars"
inputStartYearsRangeMars = [2020, 2025]
algorithmConfigsSynodicMars = ["moead", 1024, 100, False, False, 1000, True, True]
algorithmConfigsGlobalMars = ["moead", 1024, 100, False, True, 1000, True, True]

print("======= Creating relevant json files ==========")
# Create jsons for the case of running synodic period optimisations
utils.createGARunnerJsons(utils.quickConfigsMars, outputSubFolderBaseMarsSynodic, jsonSaveSubDirMars, jsonFilenameBase, inputStartYearsRangeMars,
                          utils.MarsInfoList, templateJsonPath=templateJsonPath, createSynodicJsons=True,
                          algorithmConfigs=algorithmConfigsSynodicMars)
# Create jsons for the case of running glboal optimisations and grid search
utils.createGARunnerJsons(utils.quickConfigsMars, outputSubFolderBaseMarsGlobal, jsonSaveSubDirMars, jsonFilenameBase, inputStartYearsRangeMars,
                          utils.MarsInfoList, templateJsonPath=templateJsonPath, createSynodicJsons=False,
                          algorithmConfigs=algorithmConfigsGlobalMars)


############################################ Run simulations ########################################################

runSims = True
runJupiter = False
runSaturn = False
runMars = True
printSetting=2
if runSims:
    if runJupiter: utils.runAllSimulations(jsonSaveSubDirJupiter, printSetting=printSetting)

    if runSaturn: utils.runAllSimulations(jsonSaveSubDirSaturn, printSetting=printSetting)

    if runMars: utils.runAllSimulations(jsonSaveSubDirMars, printSetting=printSetting)