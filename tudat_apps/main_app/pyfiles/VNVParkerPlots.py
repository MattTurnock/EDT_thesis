import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils

matplotlib.rcParams.update({'font.size': 20})

########################################################################################################################
##################### Regular Parker from simulation plots #############################################################
########################################################################################################################

VnVSubDir = "ParkerVnV"
VnVSaveFolder = os.path.join("pyplots", "VnV")
showing=True
loadingSimData=True

fig10Labels = ["Test Data", "V1 (1977-1989)", "V1 (1990-2004)", "V2 (1977-1989)", "V2 (1990-2004)"]

if loadingSimData:

    allSimData = utils.getAllSimDataFromFolder(VnVSubDir)
    magDataArray = allSimData[4]
    bodyDataArray = allSimData[0]

    # utils.plotMagData(magDataArray, logScaleX=False)
    # utils.plotTrajectoryData(bodyDataArray, figsize=[10,10], plotSun=True)
    utils.plotMagData(magDataArray, bodyDataArray=bodyDataArray, plotType="radius-magnitude", fignumber=10)

########################################################################################################################
##################### Voyager Data plots #############################################################
########################################################################################################################

# Set directory paths for magdata
vy1magdata_7789_dirpath = os.path.join(utils.voyager1MagDataDir, "1977-1989")
vy1magdata_9004_dirpath = os.path.join(utils.voyager1MagDataDir, "1990-2004")
vy2magdata_7789_dirpath = os.path.join(utils.voyager2MagDataDir, "1977-1989")
vy2magdata_9004_dirpath = os.path.join(utils.voyager2MagDataDir, "1990-2004")

# Load data from magdata files
vy1magdata_7789 = utils.loadVoyagerMagfieldData(vy1magdata_7789_dirpath, "7789")
vy1magdata_9004 = utils.loadVoyagerMagfieldData(vy1magdata_9004_dirpath, "9004")
vy2magdata_7789 = utils.loadVoyagerMagfieldData(vy2magdata_7789_dirpath, "7789")
vy2magdata_9004 = utils.loadVoyagerMagfieldData(vy2magdata_9004_dirpath, "9004")

# Load data from trajectory files
vy1trajdata = utils.loadVoyagerMagfieldData(utils.voyager1TrajDataDir, "traj")
vy2trajdata = utils.loadVoyagerMagfieldData(utils.voyager2TrajDataDir, "traj")

# Refine data, in the format to be used for plotting
vy1magdata_7789_refined = utils.voyagerArrayToPlotArray(vy1magdata_7789, inputType="7789")
vy1magdata_9004_refined = utils.voyagerArrayToPlotArray(vy1magdata_9004, inputType="9004", trajectoryDataArray=vy1trajdata)
vy2magdata_7789_refined = utils.voyagerArrayToPlotArray(vy2magdata_7789, inputType="7789")
vy2magdata_9004_refined = utils.voyagerArrayToPlotArray(vy2magdata_9004, inputType="9004", trajectoryDataArray=vy2trajdata)

# Calculate rolling average over interval for each
vy1magdata_7789_avg = utils.averageMagDataArray(vy1magdata_7789_refined, interval=1/12)
vy1magdata_9004_avg = utils.averageMagDataArray(vy1magdata_9004_refined, interval=1/12)
vy2magdata_7789_avg = utils.averageMagDataArray(vy2magdata_7789_refined, interval=1/12)
vy2magdata_9004_avg = utils.averageMagDataArray(vy2magdata_9004_refined, interval=1/12)

utils.plotMagData(vy1magdata_7789_refined, scatter=True, arrayType="voyager", logScaleX=False, logScaleY=False)
utils.plotMagData(vy1magdata_7789_refined, scatter=True, arrayType="voyager", logScaleX=False, logScaleY=False, plotType="radius-magnitude")


utils.plotMagData(vy1magdata_7789_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=True, fignumber=10, plotType="radius-magnitude")
utils.plotMagData(vy1magdata_9004_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=True, fignumber=10, plotType="radius-magnitude")
utils.plotMagData(vy2magdata_7789_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=True, fignumber=10, plotType="radius-magnitude")
utils.plotMagData(vy2magdata_9004_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=True, fignumber=10, plotType="radius-magnitude",
                  legend=fig10Labels, saveFolder=VnVSaveFolder, savename="Parker-R-B", xlims=[1,100], ylims=[10E-3, 10E1])


# utils.decimalYearArray2YearArray(vy1magdata_9004[:,1])

# print(utils.voyagerTrajArrayToRefined(vy1trajdata))
# utils.interpolateArrays(utils.voyagerTrajArrayToRefined(vy1trajdata), [1978.2, 1979.3, 1980.4])


if showing: plt.show()

