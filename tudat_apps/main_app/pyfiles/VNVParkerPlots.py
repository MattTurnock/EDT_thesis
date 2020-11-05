import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils

matplotlib.rcParams.update({'font.size': 20})

VnVSubDir = "ParkerVnV"
VnVSaveFolder = os.path.join("pyplots", "VnV")
showing=True

allSimData = utils.getAllSimDataFromFolder(VnVSubDir)
magDataArray = allSimData[4]
bodyDataArray = allSimData[0]

print(bodyDataArray)

utils.plotMagData(magDataArray, logScaleX=False)
utils.plotTrajectoryData(bodyDataArray, figsize=[10,10], plotSun=True)
utils.plotMagData(magDataArray, bodyDataArray=bodyDataArray, plotType="radius-magnitude", saveFolder=VnVSaveFolder, savename="Parker-R-B")




if showing: plt.show()

