from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import os
import numpy as np

# # Set logfile
logfileDir = os.path.join(utils.pythonRunnerLogfilesDir, "VNVCurrent")
utils.checkFolderExist(logfileDir)
logfilePath = os.path.join(logfileDir, utils.createLogfileName("VNVCurrent"))
utils.logger = utils.createLogger(logfilePath)


VnVSubDir_base = "currentThrustVnV_%s"
VnVSubDir_1a = VnVSubDir_base %"1a"
VnVSubDir_1b = VnVSubDir_base %"1b"
VnVSubDir_2a = VnVSubDir_base %"2a"


############################## Create relevant jsons to run, from base data for each case in utils ####################

allBaseVariables = [utils.baseVariables_1a, utils.baseVariables_1b, utils.baseVariables_2a]
utils.createCurrentVNVJsons(allBaseVariables)


############################ Do reference calculations for each set of base data #################################

currentVNV_reference_1a = utils.currentVNVCalcs(utils.baseVariables_1a)
currentVNV_reference_1b = utils.currentVNVCalcs(utils.baseVariables_1b)
currentVNV_reference_2a = utils.currentVNVCalcs(utils.baseVariables_2a)

################# Load relevant data from the simulations for each set of base data, and print ################################

currentVNV_sim_1a_outputVector = utils.getCurrentVNVDataLine(VnVSubDir_1a, "1a")
currentVNV_reference_1a.printAllVNVValues(VNVString="1a")

currentVNV_sim_1b_outputVector = utils.getCurrentVNVDataLine(VnVSubDir_1b, "1b")
currentVNV_reference_1b.printAllVNVValues(VNVString="1b")

currentVNV_sim_2a_outputVector = utils.getCurrentVNVDataLine(VnVSubDir_2a, "2a")
currentVNV_reference_2a.printAllVNVValues(VNVString="2a")

latexMatrix = np.zeros((6, len(currentVNV_sim_1a_outputVector)))
latexMatrix[0] = currentVNV_sim_1a_outputVector
latexMatrix[1] = currentVNV_reference_1a.outputVector
latexMatrix[2] = currentVNV_sim_1b_outputVector
latexMatrix[3] = currentVNV_reference_1b.outputVector
latexMatrix[4] = currentVNV_sim_2a_outputVector
latexMatrix[5] = currentVNV_reference_2a.outputVector

np.savetxt(os.path.join(utils.pyplots_dir, "LatexMatrixCurrentVNV.txt"), latexMatrix, delimiter=" & ", newline=" \\\ \n\hline\n 1S & ", fmt="%1.4f")


utils.logger.info("========== whaaaaaaa =================\n")

currentVNV_valid = utils.currentVNVCalcs(utils.baseVariables_valid)
currentVNV_valid.printAllVNVValues()




