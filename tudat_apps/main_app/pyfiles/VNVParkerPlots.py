import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils

matplotlib.rcParams.update({'font.size': 20})
nT = utils.nT

########################################################################################################################
##################### Regular Parker from simulation plots #############################################################
########################################################################################################################

VnVSubDir = "ParkerVnV"
VnVReferenceSubDir = "ParkerReferenceVnV"
VnVSaveFolder = os.path.join("pyplots", "VnV")
showing=True
doVoyager = True
doVoyagerSims = True
doVoyagerScatters = False
doReference = False

fig10VoyagerLabels = ["Tudat Simulation", "V1 (1977-1989)", "V1 (1990-2004)", "V1 (2009 - 2019)", "V2 (1977-1989)", "V2 (1990-2004)", "V2 (2009 - 2019)"]
fig11VoyagerLabels = ["Tudat Simulation", "V1 (2009 - 2019)", "V2 (2009 - 2019)"]

########################################################################################################################
##################### Voyager Data plots #############################################################
########################################################################################################################

if doVoyager:
    print("----- Doing Voyager Plots -----")

    ############################################### Magnitude plots ###################################################
    # Load and plot voyager simdata
    if doVoyagerSims:
        print("loading sim data")
        allSimData = utils.getAllSimDataFromFolder(VnVSubDir, todoList=["bodyData", "magData"])
        magDataArray = allSimData[4] / nT
        bodyDataArray = allSimData[0]

        plt.figure(1)
        # plt.plot(bodyDataArray[:, 0] / utils.AU, bodyDataArray[:, 1]/utils.AU)
        utils.plotTrajectoryData(bodyDataArray, "bodyData")
        # plt.show()


        utils.plotMagData(magDataArray, bodyDataArray=bodyDataArray, plotType="radius-magnitude", fignumber=10)
        utils.plotMagData(magDataArray, bodyDataArray=bodyDataArray, plotType="radius-magnitude", fignumber=11)

    # Set directory paths for magdata
    print("setting directory paths")
    vy1magdata_7789_dirpath = os.path.join(utils.voyager1MagDataDir, "1977-1989")
    vy1magdata_9004_dirpath = os.path.join(utils.voyager1MagDataDir, "1990-2004")
    vy1magdata_0919_dirpath =  utils.voyager1MagDataDir_48s
    vy2magdata_7789_dirpath = os.path.join(utils.voyager2MagDataDir, "1977-1989")
    vy2magdata_9004_dirpath = os.path.join(utils.voyager2MagDataDir, "1990-2004")
    vy2magdata_0919_dirpath =  utils.voyager2MagDataDir_48s

    # Load data from magdata files
    print("loading magdata")
    vy1magdata_7789 = utils.loadVoyagerMagfieldData(vy1magdata_7789_dirpath, "7789")
    vy1magdata_9004 = utils.loadVoyagerMagfieldData(vy1magdata_9004_dirpath, "9004")
    vy1magdata_0919 = utils.loadVoyagerMagfieldData(vy1magdata_0919_dirpath, "0919")
    vy2magdata_7789 = utils.loadVoyagerMagfieldData(vy2magdata_7789_dirpath, "7789")
    vy2magdata_9004 = utils.loadVoyagerMagfieldData(vy2magdata_9004_dirpath, "9004")
    vy2magdata_0919 = utils.loadVoyagerMagfieldData(vy2magdata_0919_dirpath, "0919")

    # print(vy1magdata_0919)

    # Load data from trajectory files
    print("loading trajectory data")
    vy1trajdata = utils.loadVoyagerMagfieldData(utils.voyager1TrajDataDir, "traj")
    vy2trajdata = utils.loadVoyagerMagfieldData(utils.voyager2TrajDataDir, "traj")

    # Refine data, in the format to be used for plotting
    print("refining data formats")
    vy1magdata_7789_refined = utils.voyagerArrayToPlotArray(vy1magdata_7789, inputType="7789")
    vy1magdata_9004_refined = utils.voyagerArrayToPlotArray(vy1magdata_9004, inputType="9004", trajectoryDataArray=vy1trajdata)
    vy1magdata_0919_refined = utils.voyagerArrayToPlotArray(vy1magdata_0919, inputType="0919", trajectoryDataArray=vy1trajdata)
    vy2magdata_7789_refined = utils.voyagerArrayToPlotArray(vy2magdata_7789, inputType="7789")
    vy2magdata_9004_refined = utils.voyagerArrayToPlotArray(vy2magdata_9004, inputType="9004", trajectoryDataArray=vy2trajdata)
    vy2magdata_0919_refined = utils.voyagerArrayToPlotArray(vy2magdata_0919, inputType="0919", trajectoryDataArray=vy2trajdata)

    # Calculate rolling average over interval for each
    vy1magdata_7789_avg = utils.averageMagDataArray(vy1magdata_7789_refined, interval=1/12)
    vy1magdata_9004_avg = utils.averageMagDataArray(vy1magdata_9004_refined, interval=1/12)
    vy1magdata_0919_avg = utils.averageMagDataArray(vy1magdata_0919_refined, interval=1/12)
    vy2magdata_7789_avg = utils.averageMagDataArray(vy2magdata_7789_refined, interval=1/12)
    vy2magdata_9004_avg = utils.averageMagDataArray(vy2magdata_9004_refined, interval=1/12)
    vy2magdata_0919_avg = utils.averageMagDataArray(vy2magdata_0919_refined, interval=1/12)

    if doVoyagerScatters:
        utils.plotMagData(vy1magdata_7789_refined, scatter=True, arrayType="voyager", logScaleX=False, logScaleY=False, plotType="time-magnitude")
        utils.plotMagData(vy1magdata_7789_refined, scatter=True, arrayType="voyager", logScaleX=False, logScaleY=False, plotType="radius-magnitude")


    utils.plotMagData(vy1magdata_7789_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=True, fignumber=10, plotType="radius-magnitude")
    utils.plotMagData(vy1magdata_9004_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=True, fignumber=10, plotType="radius-magnitude")
    utils.plotMagData(vy1magdata_0919_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=True, fignumber=10, plotType="radius-magnitude")
    utils.plotMagData(vy2magdata_7789_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=True, fignumber=10, plotType="radius-magnitude")
    utils.plotMagData(vy2magdata_9004_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=True, fignumber=10, plotType="radius-magnitude")
    utils.plotMagData(vy2magdata_0919_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=True, fignumber=10, plotType="radius-magnitude",
                      legend=fig10VoyagerLabels, saveFolder=VnVSaveFolder, savename="ParkerVoyager-R-B", xlims=[1, 200], ylims=[1E-2, 1E1], gridOn=True)

    # Plotting close up of IMF region
    utils.plotMagData(vy1magdata_0919_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=True, fignumber=11, plotType="radius-magnitude", color="red")
    utils.plotMagData(vy2magdata_0919_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=True, fignumber=11, plotType="radius-magnitude", color="C6",
                      legend=fig11VoyagerLabels, saveFolder=VnVSaveFolder, savename="ParkerVoyager-IMF-R-B", xlims=[75, 175], ylims=[0.05, 1E1], gridOn=True)


    ########################################## Directionality Plots ###############################################
    utils.plotMagData(vy1magdata_7789_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=False, fignumber=5, plotType="radius-BR")
    utils.plotMagData(vy1magdata_7789_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=False, fignumber=5, plotType="radius-BT")
    utils.plotMagData(vy1magdata_7789_avg, scatter=False, arrayType="voyager", logScaleX=False, logScaleY=False, fignumber=5, plotType="radius-BN")

    # plt.figure(100)
    # plt.plot(vy1magdata_7789[:, 4], vy1magdata_7789[:, 5])
    # plt.plot(vy2magdata_7789[:, 4], vy2magdata_7789[:, 5])



########################################################################################################################
##################### Reference Data plots #############################################################
########################################################################################################################

fig20ReferenceLabels = ["Tudat Simulation", "Reference Data"]

if doReference:
    print("----- Doing Reference Data Plots -----")
    # Load and plot reference sim data
    allSimDataReference = utils.getAllSimDataFromFolder(VnVReferenceSubDir)
    magDataArrayReference = allSimDataReference[4] / nT
    bodyDataArrayReference = allSimDataReference[0]

    # utils.plotMagData(magDataArray, logScaleX=False)
    # utils.plotTrajectoryData(bodyDataArray, figsize=[10,10], plotSun=True,fignumber=21)
    utils.plotMagData(magDataArrayReference, bodyDataArray=bodyDataArrayReference, plotType="radius-magnitude", fignumber=20)

    utils.plotMagData(utils.referenceParkerDataArray, scatter=False, arrayType="simple-radius-magnitude",
                      logScaleX=True, logScaleY=True, fignumber=20, plotType="radius-magnitude", xlims=[10, 200], ylims=[1E-2, 1E0],
                      legend=fig20ReferenceLabels, saveFolder=VnVSaveFolder, savename="ParkerReference-R-B", gridOn=True)

    utils.plotMagData(magDataArrayReference, bodyDataArray=bodyDataArrayReference, plotType="time-magnitude", logScaleY=True, logScaleX=False)




    # utils.decimalYearArray2YearArray(vy1magdata_9004[:,1])

    # print(utils.voyagerTrajArrayToRefined(vy1trajdata))
    # utils.interpolateArrays(utils.voyagerTrajArrayToRefined(vy1trajdata), [1978.2, 1979.3, 1980.4])


if showing: plt.show()

