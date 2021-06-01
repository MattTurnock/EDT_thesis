import numpy as np
import os
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import glob
import natsort

from matplotlib import pyplot as plt
import matplotlib

year = 365*24*60*60
AU = 1.496E11

##################################### Set some common parameters ########################################################################
loadingOnlyDepVars = False



GAJupiterSynodicDirPath = os.path.join(utils.simulation_output_dir, "GAJupiterSynodic")
GAJupiterSynodicJsonDirPath = os.path.join(utils.jsonInputs_dir, "GATestVariables_Jupiter_Stage2")
utils.GAStage2DataToNewArrayFormat(GAJupiterSynodicDirPath, GAJupiterSynodicJsonDirPath, utils.GAJupiterInitialStateFilename,
                                   loadDepVarsOnly=loadingOnlyDepVars, saveDirectory=utils.numpyBinary_dir, saveNameBase= "GAStage2NewArrayData_Jupiter_%s.npy",
                                   DVLimit=utils.launchDVLimit, launchDateRange=utils.launchDateRange,
                                 planet="Jupiter")

GASaturnSynodicDirPath = os.path.join(utils.simulation_output_dir, "GASaturnSynodic")
GASaturnSynodicJsonDirPath = os.path.join(utils.jsonInputs_dir, "GATestVariables_Saturn_Stage2")
utils.GAStage2DataToNewArrayFormat(GASaturnSynodicDirPath, GASaturnSynodicJsonDirPath, utils.GASaturnInitialStateFilename,
                                   loadDepVarsOnly=loadingOnlyDepVars, saveDirectory=utils.numpyBinary_dir, saveNameBase= "GAStage2NewArrayData_Saturn_%s.npy",
                                   DVLimit=utils.launchDVLimit, launchDateRange=utils.launchDateRange,
                                   planet="Saturn")



print("Loading numpy arrays: ")

completeJupiter = np.load("/home/matt/LinkToEDT_thesis/tudat_apps/main_app/pyfiles/numpyBinaries/GAStage2NewArrayData_Jupiter_Complete.npy", allow_pickle=True)
processedJupiter = np.load("/home/matt/LinkToEDT_thesis/tudat_apps/main_app/pyfiles/numpyBinaries/GAStage2NewArrayData_Jupiter_Processed.npy", allow_pickle=True)
completeSaturn = np.load("/home/matt/LinkToEDT_thesis/tudat_apps/main_app/pyfiles/numpyBinaries/GAStage2NewArrayData_Saturn_Complete.npy", allow_pickle=True)
processedSaturn = np.load("/home/matt/LinkToEDT_thesis/tudat_apps/main_app/pyfiles/numpyBinaries/GAStage2NewArrayData_Saturn_Processed.npy", allow_pickle=True)


print("Number sims complete: ", len(completeJupiter))
print("Number sims processed: ", len(processedJupiter))
print("")
print("Number sims complete: ", len(completeSaturn))
print("Number sims processed: ", len(processedSaturn))




