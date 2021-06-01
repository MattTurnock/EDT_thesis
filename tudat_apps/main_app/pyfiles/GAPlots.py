import numpy as np
import os
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import glob
import natsort

from matplotlib import pyplot as plt
import matplotlib

year = 365.25*24*60*60
AU = 1.496E11

##################################### Set some common parameters ########################################################################
normalising = False
matplotlib.rcParams.update({'font.size': 20})
plotFolderBase = "pyplots/GA_%s_%s-%s"
figNumberCount=1

doingJupiter=False
doingSaturn=False
doingMars=False
doingJupiterStage2 = True
doingSaturnStage2 = True
doingProcessing = True
# justPrintingStage2Info = True

doingTitles = False
showing=False

plotTitles = {"Jupiter_DVTOFS_Syn": None,
              "Jupiter_yearsTOFS_Syn": None,
              "Jupiter_DVTOFS_Glob": None,
              "Jupiter_yearsTOFS_Glob": None,
              "Jupiter_porkchop_Glob": None,
              "Jupiter_porkchop_Syn": None,
              "Jupiter_porkchopClose_Glob": None,

              "Saturn_DVTOFS_Syn": None,
              "Saturn_yearsTOFS_Syn": None,
              "Saturn_DVTOFS_Glob": None,
              "Saturn_yearsTOFS_Glob": None,
              "Saturn_porkchop_Glob": None,
              "Saturn_porkchop_Syn": None,
              "Saturn_porkchopClose_Glob": None,

              }

if doingTitles:
    plotTitles["Jupiter_DVTOFS_Syn"] = "Jupiter DV-TOF Synodic"
    plotTitles["Jupiter_yearsTOFS_Syn"] = "Jupiter Launch Year-TOF Synodic"
    plotTitles["Jupiter_DVTOFS_Glob"] = "Jupiter DV-TOF Global"
    plotTitles["Jupiter_yearsTOFS_Glob"] = "Jupiter Launch Year-TOF Global"
    plotTitles["Jupiter_porkchop_Glob"] = "Jupiter Porkchop Plot Global"
    plotTitles["Jupiter_porkchop_Syn"] = "Jupiter Porkchop Plot Synodic"
    plotTitles["Jupiter_porkchopClose_Glob"] = "Jupiter Porkchop Plot Global Close-Up"

    plotTitles["Saturn_DVTOFS_Syn"] = "Saturn DV-TOF Synodic"
    plotTitles["Saturn_yearsTOFS_Syn"] = "Saturn Launch Year-TOF Synodic"
    plotTitles["Saturn_DVTOFS_Glob"] = "Saturn DV-TOF Global"
    plotTitles["Saturn_yearsTOFS_Glob"] = "Saturn Launch Year-TOF Global"
    plotTitles["Saturn_porkchop_Glob"] = "Saturn Porkchop Plot Global"
    plotTitles["Saturn_porkchop_Syn"] = "Saturn Porkchop Plot Synodic"
    plotTitles["Saturn_porkchopClose_Glob"] = "Saturn Porkchop Plot Global Close-Up"




scatterMarkerPorkchops = "x"
scatterLinewidthsPorkchops = 1
scatterPointSizePorkchops = 1.5
scatterColorPorkchops="C9"

if doingJupiter:
    ##################################### Do Jupiter Plots ################################################################

    # Set Jupiter-specific parameters
    quickConfigsJupiter = utils.quickConfigsJupiter
    plotFolderJupiter = plotFolderBase %(quickConfigsJupiter[0], quickConfigsJupiter[1], quickConfigsJupiter[2])
    JupiterSubfolderSynodic = "GAJupiterSynodic"
    # JupiterSubSubfolderBase = "GAJupiter/GACalculatorNominal_%s_%s-%s/"

    # Calculate allYears data and set numpy arrays to a name
    allYearsJupiterDataSynodic = utils.getAllYearsGA(quickConfigsJupiter, JupiterSubfolderSynodic)
    fitnessFileDVsJupiterSynodic, fitnessFileTOFSJupiterSynodic, launchYearsJupiterSynodic, arrivalYearsJupiterSynodic, startYearsJupiterSynodic, endYearsJupiterSynodic = allYearsJupiterDataSynodic[0:6]

    utils.plotManyDataGA(allYearsJupiterDataSynodic, figNumberCount, quickConfigsJupiter, plotType="DV-TOFS", saveFolder=plotFolderJupiter, savenameSuffix="_Synodic", plotTitle=plotTitles["Jupiter_DVTOFS_Syn"])
    figNumberCount += 1
    utils.plotManyDataGA(allYearsJupiterDataSynodic, figNumberCount, quickConfigsJupiter, plotType="launchYears-TOFS", saveFolder=plotFolderJupiter, savenameSuffix="_Synodic", plotTitle=plotTitles["Jupiter_yearsTOFS_Syn"])
    figNumberCount += 1

    # Set Jupiter-specific parameters
    quickConfigsJupiter = utils.quickConfigsJupiter
    plotFolderJupiter = plotFolderBase %(quickConfigsJupiter[0], quickConfigsJupiter[1], quickConfigsJupiter[2])
    JupiterSubfolderGlobal = "GAJupiterGlobal"
    # JupiterSubSubfolderBase = "GAJupiter/GACalculatorNominal_%s_%s-%s/"

    # Calculate allYears data and set numpy arrays to a name
    allYearsJupiterDataGlobal = utils.getAllYearsGA(quickConfigsJupiter, JupiterSubfolderGlobal)
    fitnessFileDVsJupiterGlobal, fitnessFileTOFSJupiterGlobal, launchYearsJupiterGlobal, arrivalYearsJupiterGlobal, startYearsJupiterGlobal, endYearsJupiterGlobal = allYearsJupiterDataGlobal[0:6]

    utils.plotManyDataGA(allYearsJupiterDataGlobal, figNumberCount, quickConfigsJupiter, plotType="DV-TOFS", saveFolder=plotFolderJupiter, savenameSuffix="_Global", plotTitle=plotTitles["Jupiter_DVTOFS_Glob"])
    figNumberCount += 1
    utils.plotManyDataGA(allYearsJupiterDataGlobal, figNumberCount, quickConfigsJupiter, plotType="launchYears-TOFS", saveFolder=plotFolderJupiter, savenameSuffix="_Global", plotTitle=plotTitles["Jupiter_yearsTOFS_Glob"])
    figNumberCount += 1

    jupiterPorkchopPathList = os.listdir(os.path.join(utils.simulation_output_dir, JupiterSubfolderGlobal))
    if len(jupiterPorkchopPathList) != 1:
        print("ERROR: should only be one entry in porkchop (global) folder")
        exit()
    else:
        jupiterPorkchopPath = os.path.join(utils.simulation_output_dir, JupiterSubfolderGlobal,  jupiterPorkchopPathList[0])

    # jupiterPorkchopSaveDirectory = os.path.join(utils.pyplots_dir, jupiterPorkchopPathList[0])
    # jupiterPorkchopSavename = "porkchopJupiter_2020-2050"
    # utils.porkchopPlot(jupiterPorkchopPath, "porkchop_GA_EJ", contourCount=5, contourLinewidths=0.1,
    #                    fignumber=5, saveDirectory=jupiterPorkchopSaveDirectory, saveName=jupiterPorkchopSavename, xlims=None)
    # utils.plotManyDataGA(allYearsJupiterDataSynodic, 5, quickConfigsJupiter, plotType="launchYears-TOFS", scatterPointSize=1,
    #                      saveFolder=jupiterPorkchopSaveDirectory, savenameOverride=jupiterPorkchopSavename, plotLegend=False)

    JupiterPorkchopSaveDirectory = os.path.join(utils.pyplots_dir, jupiterPorkchopPathList[0])
    # Full porkchop global opt
    JupiterPorkchopSavenameGlobal = "porkchopJupiterGlobal_2020-2050"
    utils.porkchopPlot(jupiterPorkchopPath, "porkchop_GA_EJ", contourCount=5, contourLinewidths=0.1,
                       fignumber=figNumberCount, saveDirectory=JupiterPorkchopSaveDirectory, saveName=JupiterPorkchopSavenameGlobal,
                       xlims=[2020, 2050])
    utils.plotManyDataGA(allYearsJupiterDataGlobal, figNumberCount, utils.quickConfigsJupiter, plotType="launchYears-TOFS", scatterPointSize=scatterPointSizePorkchops,
                         saveFolder=JupiterPorkchopSaveDirectory, savenameOverride=JupiterPorkchopSavenameGlobal, plotLegend=False,
                         scatterMarker=scatterMarkerPorkchops, scatterLinewidths=scatterLinewidthsPorkchops, scatterColour=scatterColorPorkchops, plotTitle=plotTitles["Jupiter_porkchop_Glob"])
    figNumberCount += 1

    # Full porkchop synodic opt
    JupiterPorkchopSavenameSynodic = "porkchopJupiterSynodic_2020-2050"
    utils.porkchopPlot(jupiterPorkchopPath, "porkchop_GA_EJ", contourCount=5, contourLinewidths=0.1,
                       fignumber=figNumberCount, saveDirectory=JupiterPorkchopSaveDirectory, saveName=JupiterPorkchopSavenameSynodic,
                       xlims=[2020, 2050])
    utils.plotManyDataGA(allYearsJupiterDataSynodic, figNumberCount, utils.quickConfigsJupiter, plotType="launchYears-TOFS", scatterPointSize=scatterPointSizePorkchops,
                         saveFolder=JupiterPorkchopSaveDirectory, savenameOverride=JupiterPorkchopSavenameSynodic, plotLegend=False,
                         scatterMarker=scatterMarkerPorkchops, scatterLinewidths=scatterLinewidthsPorkchops, scatterColour=None, plotTitle=plotTitles["Jupiter_porkchop_Syn"])
    figNumberCount += 1

    # Closeup porkchop global opt
    JupiterPorkchopSavenameGlobalClose = "porkchopJupiterGlobalClose_2020-2050"
    utils.porkchopPlot(jupiterPorkchopPath, "porkchop_GA_EJ", contourCount=5, contourLinewidths=0.1,
                       fignumber=figNumberCount, saveDirectory=JupiterPorkchopSaveDirectory, saveName=JupiterPorkchopSavenameGlobalClose,
                       xlims=[2029, 2033])
    utils.plotManyDataGA(allYearsJupiterDataGlobal, figNumberCount, utils.quickConfigsJupiter, plotType="launchYears-TOFS", scatterPointSize=scatterPointSizePorkchops,
                         saveFolder=JupiterPorkchopSaveDirectory, savenameOverride=JupiterPorkchopSavenameGlobalClose, plotLegend=False,
                         scatterMarker=scatterMarkerPorkchops, scatterLinewidths=scatterLinewidthsPorkchops, scatterColour=scatterColorPorkchops, plotTitle=plotTitles["Jupiter_porkchopClose_Glob"])
    figNumberCount += 1

if doingSaturn:
    ##################################### Do Saturn Plots ################################################################

    # Set Saturn-specific parameters
    plotFolderSaturn = plotFolderBase %(utils.quickConfigsSaturn[0], utils.quickConfigsSaturn[1], utils.quickConfigsSaturn[2])
    SaturnSubfolderSynodic = "GASaturnSynodic"


    # Calculate allYears data and set numpy arrays to a name
    allYearsSaturnDataSynodic = utils.getAllYearsGA(utils.quickConfigsSaturn, SaturnSubfolderSynodic)
    fitnessFileDVsSaturnSynodic, fitnessFileTOFSSaturnSynodic, launchYearsSaturnSynodic, arrivalYearsSaturnSynodic, startYearsSaturnSynodic, endYearsSaturnSynodic = allYearsSaturnDataSynodic[0:6]

    utils.plotManyDataGA(allYearsSaturnDataSynodic, figNumberCount, utils.quickConfigsSaturn, plotType="DV-TOFS", saveFolder=plotFolderSaturn, savenameSuffix="_Synodic", plotTitle=plotTitles["Saturn_DVTOFS_Syn"])
    figNumberCount += 1
    utils.plotManyDataGA(allYearsSaturnDataSynodic, figNumberCount, utils.quickConfigsSaturn, plotType="launchYears-TOFS", saveFolder=plotFolderSaturn, savenameSuffix="_Synodic", plotTitle=plotTitles["Saturn_yearsTOFS_Syn"])
    figNumberCount += 1

    # Set Saturn-specific parameters
    SaturnSubfolderGlobal = "GASaturnGlobal"


    # Calculate allYears data and set numpy arrays to a name
    allYearsSaturnDataGlobal = utils.getAllYearsGA(utils.quickConfigsSaturn, SaturnSubfolderGlobal)
    fitnessFileDVsSaturnGlobal, fitnessFileTOFSSaturnGlobal, launchYearsSaturnGlobal, arrivalYearsSaturnGlobal, startYearsSaturnGlobal, endYearsSaturnGlobal = allYearsSaturnDataGlobal[0:6]

    utils.plotManyDataGA(allYearsSaturnDataGlobal, figNumberCount, utils.quickConfigsSaturn, plotType="DV-TOFS", saveFolder=plotFolderSaturn, savenameSuffix="_Global", plotTitle=plotTitles["Saturn_DVTOFS_Glob"])
    figNumberCount += 1
    utils.plotManyDataGA(allYearsSaturnDataGlobal, figNumberCount, utils.quickConfigsSaturn, plotType="launchYears-TOFS", saveFolder=plotFolderSaturn, savenameSuffix="_Global", plotTitle=plotTitles["Saturn_yearsTOFS_Glob"])
    figNumberCount += 1

    SaturnPorkchopPathList = os.listdir(os.path.join(utils.simulation_output_dir, SaturnSubfolderGlobal))
    if len(SaturnPorkchopPathList) != 1:
        print("ERROR: should only be one entry in porkchop (global) folder")
        exit()
    else:
        SaturnPorkchopPath = os.path.join(utils.simulation_output_dir, SaturnSubfolderGlobal,  SaturnPorkchopPathList[0])

    SaturnPorkchopSaveDirectory = os.path.join(utils.pyplots_dir, SaturnPorkchopPathList[0])
    # Full porkchop global opt
    SaturnPorkchopSavenameGlobal = "porkchopSaturnGlobal_2020-2050"
    utils.porkchopPlot(SaturnPorkchopPath, "porkchop_GA_ES", contourCount=5, contourLinewidths=0.1,
                       fignumber=figNumberCount, saveDirectory=SaturnPorkchopSaveDirectory, saveName=SaturnPorkchopSavenameGlobal,
                       xlims=[2020, 2050])
    utils.plotManyDataGA(allYearsSaturnDataGlobal, figNumberCount, utils.quickConfigsSaturn, plotType="launchYears-TOFS", scatterPointSize=scatterPointSizePorkchops,
                         saveFolder=SaturnPorkchopSaveDirectory, savenameOverride=SaturnPorkchopSavenameGlobal, plotLegend=False,
                         scatterMarker=scatterMarkerPorkchops, scatterLinewidths=scatterLinewidthsPorkchops, scatterColour=scatterColorPorkchops, plotTitle=plotTitles["Saturn_porkchop_Glob"])
    figNumberCount += 1

    # Full porkchop synodic opt
    SaturnPorkchopSavenameSynodic = "porkchopSaturnSynodic_2020-2050"
    utils.porkchopPlot(SaturnPorkchopPath, "porkchop_GA_ES", contourCount=5, contourLinewidths=0.1,
                       fignumber=figNumberCount, saveDirectory=SaturnPorkchopSaveDirectory, saveName=SaturnPorkchopSavenameSynodic,
                       xlims=[2020, 2050])
    utils.plotManyDataGA(allYearsSaturnDataSynodic, figNumberCount, utils.quickConfigsSaturn, plotType="launchYears-TOFS", scatterPointSize=scatterPointSizePorkchops,
                         saveFolder=SaturnPorkchopSaveDirectory, savenameOverride=SaturnPorkchopSavenameSynodic, plotLegend=False,
                         scatterMarker=scatterMarkerPorkchops, scatterLinewidths=scatterLinewidthsPorkchops, scatterColour=None, plotTitle=plotTitles["Saturn_porkchop_Syn"])
    figNumberCount += 1

    # Closeup porkchop global opt
    SaturnPorkchopSavenameGlobalClose = "porkchopSaturnGlobalClose_2020-2050"
    utils.porkchopPlot(SaturnPorkchopPath, "porkchop_GA_ES", contourCount=5, contourLinewidths=0.1,
                       fignumber=figNumberCount, saveDirectory=SaturnPorkchopSaveDirectory, saveName=SaturnPorkchopSavenameGlobalClose,
                       xlims=[2028, 2034])
    utils.plotManyDataGA(allYearsSaturnDataGlobal, figNumberCount, utils.quickConfigsSaturn, plotType="launchYears-TOFS", scatterPointSize=scatterPointSizePorkchops,
                         saveFolder=SaturnPorkchopSaveDirectory, savenameOverride=SaturnPorkchopSavenameGlobalClose, plotLegend=False,
                         scatterMarker=scatterMarkerPorkchops, scatterLinewidths=scatterLinewidthsPorkchops, scatterColour=scatterColorPorkchops, plotTitle=plotTitles["Saturn_porkchopClose_Glob"])
    figNumberCount += 1

if doingMars:
    ##################################### Do Mars Plots ################################################################

    # Set Mars-specific parameters
    plotFolderMars = plotFolderBase %(utils.quickConfigsMars[0], utils.quickConfigsMars[1], utils.quickConfigsMars[2])
    MarsSubfolderSynodic = "GAMarsSynodic"


    # # Calculate allYears data and set numpy arrays to a name
    # allYearsMarsDataSynodic = utils.getAllYearsGA(utils.quickConfigsMars, MarsSubfolderSynodic)
    # fitnessFileDVsMarsSynodic, fitnessFileTOFSMarsSynodic, launchYearsMarsSynodic, arrivalYearsSaturnSynodic, startYearsSaturnSynodic, endYearsSaturnSynodic = allYearsSaturnDataSynodic[0:6]
    #
    # utils.plotManyDataGA(allYearsSaturnDataSynodic, figNumberCount, utils.quickConfigsSaturn, plotType="DV-TOFS", saveFolder=plotFolderSaturn, savenameSuffix="_Synodic")
    # figNumberCount += 1
    # utils.plotManyDataGA(allYearsSaturnDataSynodic, figNumberCount, utils.quickConfigsSaturn, plotType="launchYears-TOFS", saveFolder=plotFolderSaturn, savenameSuffix="_Synodic")
    # figNumberCount += 1
    #
    # Set Mars-specific parameters
    MarsSubfolderGlobal = "GAMarsGlobal"


    # Calculate allYears data and set numpy arrays to a name
    allYearsMarsDataGlobal = utils.getAllYearsGA(utils.quickConfigsMars, MarsSubfolderGlobal)
    fitnessFileDVsMarsGlobal, fitnessFileTOFSMarsGlobal, launchYearsMarsGlobal, arrivalYearsMarsGlobal, startYearsMarsGlobal, endYearsMarsGlobal = allYearsMarsDataGlobal[0:6]

    utils.plotManyDataGA(allYearsMarsDataGlobal, figNumberCount, utils.quickConfigsMars,
                         plotType="DV-TOFS", saveFolder=plotFolderMars, savenameSuffix="_Global", removeDominated=False,
                          TOFUnits="Days", ylims=[200,1000], xlims=[5.5, 10], printMinDV=True)
    figNumberCount += 1
    utils.plotManyDataGA(allYearsMarsDataGlobal, figNumberCount, utils.quickConfigsMars,
                         plotType="launchYears-TOFS", saveFolder=plotFolderMars, savenameSuffix="_Global", removeDominated=False,
                          TOFUnits="Days", ylims=[200,1000], printMinDV=True)
    figNumberCount += 1

    MarsPorkchopPathList = os.listdir(os.path.join(utils.simulation_output_dir, MarsSubfolderGlobal))
    if len(MarsPorkchopPathList) != 1:
        print("ERROR: should only be one entry in porkchop (global) folder")
        exit()
    else:
        MarsPorkchopPath = os.path.join(utils.simulation_output_dir, MarsSubfolderGlobal,  MarsPorkchopPathList[0])

    MarsPorkchopSaveDirectory = os.path.join(utils.pyplots_dir, MarsPorkchopPathList[0])
    # Full porkchop global opt
    MarsPorkchopSavenameGlobal = "porkchopMarsGlobal_2020-2050"
    utils.porkchopPlot(MarsPorkchopPath, "porkchop_GA_EM", contourCount=25, contourLinewidths=1,
                       fignumber=figNumberCount, saveDirectory=MarsPorkchopSaveDirectory, saveName=MarsPorkchopSavenameGlobal,
                       xlims=[2020, 2025], ylims=[200,1000], TOFUnit="Days")
    utils.plotManyDataGA(allYearsMarsDataGlobal, figNumberCount, utils.quickConfigsMars, plotType="launchYears-TOFS", scatterPointSize=1,
                         saveFolder=MarsPorkchopSaveDirectory, savenameOverride=MarsPorkchopSavenameGlobal, plotLegend=False, TOFUnits="Days",
                         scatterMarker=None, scatterLinewidths=None, scatterColour=scatterColorPorkchops, removeDominated=False)
    figNumberCount += 1



    # # Full porkchop synodic opt
    # MarsPorkchopSavenameSynodic = "porkchopMarsSynodic_2020-2050"
    # utils.porkchopPlot(MarsPorkchopPath, "porkchop_GA_ES", contourCount=5, contourLinewidths=0.1,
    #                    fignumber=figNumberCount, saveDirectory=MarsPorkchopSaveDirectory, saveName=MarsPorkchopSavenameSynodic,
    #                    xlims=[2020, 2050])
    # utils.plotManyDataGA(allYearsMarsDataSynodic, figNumberCount, utils.quickConfigsMars, plotType="launchYears-TOFS", scatterPointSize=1,
    #                      saveFolder=MarsPorkchopSaveDirectory, savenameOverride=MarsPorkchopSavenameSynodic, plotLegend=False)
    # figNumberCount += 1
    #
    # # Closeup porkchop global opt
    # SaturnPorkchopSavenameGlobalClose = "porkchopSaturnGlobalClose_2020-2050"
    # utils.porkchopPlot(SaturnPorkchopPath, "porkchop_GA_ES", contourCount=5, contourLinewidths=0.1,
    #                    fignumber=figNumberCount, saveDirectory=SaturnPorkchopSaveDirectory, saveName=SaturnPorkchopSavenameGlobalClose,
    #                    xlims=[2028, 2034])
    # utils.plotManyDataGA(allYearsSaturnDataGlobal, figNumberCount, utils.quickConfigsSaturn, plotType="launchYears-TOFS", scatterPointSize=1,
    #                      saveFolder=SaturnPorkchopSaveDirectory, savenameOverride=SaturnPorkchopSavenameGlobalClose, plotLegend=False)
    # figNumberCount += 1


def saveGAPlotNumpyData(GAPlotNumpyFolder, planet, emptyNumpyDirectory=False):

    utils.checkFolderExist(GAPlotNumpyFolder, emptyDirectory=emptyNumpyDirectory)

    if planet == "Jupiter":
        initialStateFilename = utils.GAJupiterInitialStateFilename
    elif planet == "Saturn":
        initialStateFilename = utils.GASaturnInitialStateFilename

    synodicDirPath = os.path.join(utils.simulation_output_dir, "GA%sSynodic" %planet)
    synodicSubDirList = utils.natsort.natsorted(glob.glob(synodicDirPath + "/*"))
    synodicSubDirListToUse = []

    synodicJsonDirPath = os.path.join(utils.jsonInputs_dir, "GATestVariables_%s_Stage2" %planet)
    synodicJsonDirList = utils.natsort.natsorted(glob.glob(synodicJsonDirPath + "/*"))

    listOfAllSimDataLists = []
    listOfCloseApproachArrays = []
    listOfLaunchDatesArrays = []
    listOfInitialStateLists = []


    bar = utils.IncrementalBar("Loading data ", max=len(synodicSubDirList), suffix='[%(percent).1f%%. ETA: %(eta)ds]   ')
    for i in range(len(synodicSubDirList)):

        subDirPath = synodicSubDirList[i]
        initialStates = np.genfromtxt(os.path.join(subDirPath, initialStateFilename), delimiter=",")

        stage2JsonsList = utils.natsort.natsorted(glob.glob(synodicJsonDirList[i] + "/*"))

        thisAllSimDataList = []
        closeApproachList = []
        launchDatesList = []
        initialStatesList = [] # These are the ones we will use in the future ok

        for j in range(len(stage2JsonsList)):

            allSimData = utils.getAllSimDataFromJson(stage2JsonsList[j], todoList=["depVarData"], printInfo=False)
            depVarData = allSimData[2]

            initialState = initialStates[j]
            deltaV = initialState[8]



            closeApproachValue = np.amin(depVarData[:, 14])
            launchDateValue = depVarData[0, 0]


            # Check if launch date and delta V within specs
            # print(utils.launchDateRange[0])
            # print(2000+ launchDateValue/ year)
            # print(utils.launchDateRange[1])
            # print(deltaV)
            # print(utils.launchDVLimit, "\n")
            if (utils.launchDateRange[0] < 2000 + launchDateValue/year < utils.launchDateRange[1]) and (deltaV / 1000 < utils.launchDVLimit):

                thisAllSimDataList.append(allSimData)
                closeApproachList.append(closeApproachValue)
                launchDatesList.append(launchDateValue)
                initialStatesList.append(initialState)

            # LOAD DATA FOR THE JSONS

        # Do not append empty arrays
        if len(thisAllSimDataList) != 0:
            listOfAllSimDataLists.append(thisAllSimDataList)
            listOfCloseApproachArrays.append(closeApproachList)
            listOfLaunchDatesArrays.append(launchDatesList)
            listOfInitialStateLists.append(initialStatesList)
            synodicSubDirListToUse.append(synodicSubDirList[i])

        bar.next()
    bar.finish()

    print("Saving numpy arrays")

    np.save(os.path.join(GAPlotNumpyFolder, "GA%sListOfAllSimDataLists.npy" %planet), listOfAllSimDataLists)
    np.save(os.path.join(GAPlotNumpyFolder, "GA%sListOfCloseApproachArrays.npy" %planet), listOfCloseApproachArrays)
    np.save(os.path.join(GAPlotNumpyFolder, "GA%sListOfLaunchDatesArrays.npy" %planet), listOfLaunchDatesArrays)
    np.save(os.path.join(GAPlotNumpyFolder, "GA%sListOfInitialStateLists.npy" %planet), listOfInitialStateLists)
    np.save(os.path.join(GAPlotNumpyFolder, "GA%sSynodicSubDirListToUse.npy" %planet), synodicSubDirListToUse)


def trueFlatten(inputArray):
    newList = []
    for i in range(len(inputArray)):
        # print(inputArray[i])
        newList = newList + list(inputArray[i])

    return np.array(newList)




GAPlotNumpyFolder = os.path.join(utils.numpyBinary_dir, "GAPlotsData/")
if doingProcessing:
    saveGAPlotNumpyData(GAPlotNumpyFolder, "Saturn", emptyNumpyDirectory=True)
    saveGAPlotNumpyData(GAPlotNumpyFolder, "Jupiter", emptyNumpyDirectory=False)



print("Loading numpy arrays")
GAJupiterListOfAllSimDataLists = np.load(os.path.join(GAPlotNumpyFolder, "GAJupiterListOfAllSimDataLists.npy"), allow_pickle=True)
GAJupiterListOfCloseApproachArrays = np.load(os.path.join(GAPlotNumpyFolder, "GAJupiterListOfCloseApproachArrays.npy"), allow_pickle=True)
GAJupiterListOfLaunchDatesArrays = np.load(os.path.join(GAPlotNumpyFolder, "GAJupiterListOfLaunchDatesArrays.npy"), allow_pickle=True)
GAJupiterListOfInitialStateLists = np.load(os.path.join(GAPlotNumpyFolder, "GAJupiterListOfInitialStateLists.npy"), allow_pickle=True)
GAJupiterSynodicSubDirListToUse = np.load(os.path.join(GAPlotNumpyFolder, "GAJupiterSynodicSubDirListToUse.npy"), allow_pickle=True)

GASaturnListOfAllSimDataLists = np.load(os.path.join(GAPlotNumpyFolder, "GASaturnListOfAllSimDataLists.npy"), allow_pickle=True)
GASaturnListOfCloseApproachArrays = np.load(os.path.join(GAPlotNumpyFolder, "GASaturnListOfCloseApproachArrays.npy"), allow_pickle=True)
GASaturnListOfLaunchDatesArrays = np.load(os.path.join(GAPlotNumpyFolder, "GASaturnListOfLaunchDatesArrays.npy"), allow_pickle=True)
GASaturnListOfInitialStateLists = np.load(os.path.join(GAPlotNumpyFolder, "GASaturnListOfInitialStateLists.npy"), allow_pickle=True)
GASaturnSynodicSubDirListToUse = np.load(os.path.join(GAPlotNumpyFolder, "GASaturnSynodicSubDirListToUse.npy"), allow_pickle=True)


stage2PyplotFolder = "pyplots/GA_Stage2/"

jupiterStage2PyplotSavenameBase = "GA_Jupiter_Stage2_%s.png"
saturnStage2PyplotSavenameBase = "GA_Saturn_Stage2_%s.png"

print(len(GAJupiterListOfCloseApproachArrays))
print("Largest Jupiter separation [AU]: ", np.amax(trueFlatten(GAJupiterListOfCloseApproachArrays)) / AU)
print("", len(GAJupiterListOfCloseApproachArrays), "----", len(GASaturnListOfCloseApproachArrays))

# for i in range(len(GAJupiterListOfCloseApproachArrays)):
#     print(len(GAJupiterListOfAllSimDataLists[i]), "    ", len(GAJupiterListOfCloseApproachArrays[i]), "    ", len(GAJupiterListOfLaunchDatesArrays[i]), "    ", len(GAJupiterListOfInitialStateLists[i]), "    ", len(GAJupiterSynodicSubDirListToUse[i]))
# print("")
# for i in range(len(GASaturnListOfCloseApproachArrays)):
#     print(len(GASaturnListOfAllSimDataLists[i]), "    ", len(GASaturnListOfCloseApproachArrays[i]), "    ", len(GASaturnListOfLaunchDatesArrays[i]), "    ", len(GASaturnListOfInitialStateLists[i]), "    ", len(GASaturnSynodicSubDirListToUse[i]))


print(len(GASaturnListOfCloseApproachArrays))
print("Largest Saturn separation [m]: ", np.amax(trueFlatten(GASaturnListOfCloseApproachArrays)) )

utils.plotManyStage2GAData(GAJupiterListOfAllSimDataLists, GAJupiterListOfInitialStateLists, plotType="DV-FinalAlt", removeDominated=True, plotParetoFront=False,
                           saveFolder=stage2PyplotFolder, savename=jupiterStage2PyplotSavenameBase %"Syn-DV-Alt-Pareto", scatterPointSize=20, plotLegend=False, fignumber=1, plotGrid=False)
utils.plotManyStage2GAData(GAJupiterListOfAllSimDataLists, GAJupiterListOfInitialStateLists, plotType="DV-FinalAlt", removeDominated=True, plotParetoFront=True,
                           saveFolder=stage2PyplotFolder, savename=jupiterStage2PyplotSavenameBase %"Syn-DV-Alt-Pareto", scatterPointSize=10, plotLegend=False, fignumber=1, scatterLinewidths=0.5)
utils.plotManyStage2GAData(GAJupiterListOfAllSimDataLists, GAJupiterListOfInitialStateLists, plotType="DV-FinalAlt", removeDominated=False, plotParetoFront=False,
                               saveFolder=stage2PyplotFolder, savename=jupiterStage2PyplotSavenameBase %"Syn-DV-Alt-Scatter", scatterPointSize=2, plotLegend=False, fignumber=2)

print("")

# utils.plotManyStage2GAData(GASaturnListOfAllSimDataLists, GASaturnListOfInitialStateLists, plotType="DV-FinalAlt", removeDominated=True, plotParetoFront=True,
#                            saveFolder=stage2PyplotFolder, savename=saturnStage2PyplotSavenameBase %"Syn-DV-Alt-Pareto", scatterPointSize=2, initialStateFilename=utils.GASaturnInitialStateFilename, plotLegend=False)

utils.plotManyStage2GAData(GASaturnListOfAllSimDataLists, GASaturnListOfInitialStateLists, plotType="DV-FinalAlt", removeDominated=True, plotParetoFront=False,
                           saveFolder=stage2PyplotFolder, savename=saturnStage2PyplotSavenameBase %"Syn-DV-Alt-Pareto", scatterPointSize=20, plotLegend=False, fignumber=3, plotGrid=False)
utils.plotManyStage2GAData(GASaturnListOfAllSimDataLists, GASaturnListOfInitialStateLists, plotType="DV-FinalAlt", removeDominated=True, plotParetoFront=True,
                           saveFolder=stage2PyplotFolder, savename=saturnStage2PyplotSavenameBase %"Syn-DV-Alt-Pareto", scatterPointSize=10, plotLegend=False, fignumber=3, scatterLinewidths=0.5)


utils.plotManyStage2GAData(GASaturnListOfAllSimDataLists, GASaturnListOfInitialStateLists, plotType="DV-FinalAlt", removeDominated=False, plotParetoFront=False,
                           saveFolder=stage2PyplotFolder, savename=saturnStage2PyplotSavenameBase %"Syn-DV-Alt-Scatter", scatterPointSize=2, initialStateFilename=utils.GASaturnInitialStateFilename, plotLegend=False, fignumber=4)



if showing:
    plt.show()


