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
plotFolder = "pyplots/InO"
utils.checkFolderExist(plotFolder)

showing = False
useExample = False
plotOptimalTrajectory = True

generationToUse = -1

if useExample:
    results_DirPathToUse = os.path.join(utils.numpyBinary_dir, "InO_EXAMPLE")
else:
    results_DirPathToUse = utils.InO_Results_DirPath
utils.checkFolderExist(results_DirPathToUse)

print("-- Loading Numpy Array --")
numpyArraysFileList = natsort.natsorted(glob.glob(results_DirPathToUse + "/*"))
numpyArraysFileList_Processed = natsort.natsorted(glob.glob(results_DirPathToUse + "/*Processed*"))
# print(numpyArraysFileList_Processed)

# InO_Populations
fitArrayList = []
popArrayList = []
for i in range(len(numpyArraysFileList_Processed)):
    filePath = numpyArraysFileList_Processed[i]

    if "Fit" in filePath:
        fitArrayList.append(np.load(filePath, allow_pickle=True))

    elif "Pop" in filePath:
        popArrayList.append(np.load(filePath, allow_pickle=True))

# print(fitArrayList[-1])
plt.figure()

if generationToUse != "all":

    fitArrayToUse = fitArrayList[generationToUse]
    popArrayToUse = popArrayList[generationToUse]
    # plt.scatter(fitArrayToUse[:, 1]/ 1000, fitArrayToUse[:, 0], 1 )
    plt.scatter(np.zeros(len(fitArrayToUse)), fitArrayToUse[:, 0], 1 )

else:
    for i in range(len(fitArrayList)):
        # plt.scatter(fitArrayList[i][:, 1] / 1000, fitArrayList[i][:, 0], 1)
        plt.scatter(np.zeros(len(fitArrayList[i])), fitArrayList[i][:, 0], 1)
        print(fitArrayList[i])
        print(popArrayList[i])
        print("")

plt.xlabel("N/A [km/s]")
plt.ylabel("Maximum Ap [AU]")


gen0StartYears = popArrayList[0][:, 0]
gen0InitialAp = popArrayList[0][:, 1]
gen0InitialPe = popArrayList[0][:, 2]
gen0MaxAps = np.reshape(fitArrayList[0], np.shape(gen0InitialPe))
gen0TargetPe = popArrayList[0][:, 3]

genFinalStartYears = popArrayList[-1][:, 0]
genFinalInitialAp = popArrayList[-1][:, 1]
genFinalInitialPe = popArrayList[-1][:, 2]
genFinalMaxAps = np.reshape(fitArrayList[-1], np.shape(genFinalInitialPe))
genFinalTargetPe = popArrayList[-1][:, 3]

scatterSize = 3

legend = ["First generation", "Final generation"]

# Launch year vs Max Ap
plt.figure(figsize=utils.figSizeDefault)
plt.scatter(gen0StartYears, gen0MaxAps, scatterSize)
plt.scatter(genFinalStartYears, genFinalMaxAps, scatterSize)
plt.xlabel("Start Year")
plt.ylabel("Maximum Ap [AU]")
plt.grid()
plt.legend(legend)
plt.savefig(os.path.join(utils.pyplots_dir, "InO/Year_MaxAp.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(utils.pyplots_dir, "InO/Year_MaxAp.png"), bbox_inches="tight")

# initial Pe vs max Ap
plt.figure(figsize=utils.figSizeDefault)
plt.scatter(gen0InitialPe, gen0MaxAps - gen0InitialAp, scatterSize)
plt.scatter(genFinalInitialPe, genFinalMaxAps - genFinalInitialAp, scatterSize)
plt.xlabel("Initial Pe [AU]")
plt.ylabel("Ap Change [AU]")
plt.grid()
plt.legend(legend)
plt.savefig(os.path.join(utils.pyplots_dir, "InO/InitialPe_MaxAp.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(utils.pyplots_dir, "InO/InitialPe_MaxAp.png"), bbox_inches="tight")

# target Pe vs max Ap
plt.figure(figsize=utils.figSizeDefault)
plt.scatter(gen0TargetPe, gen0MaxAps - gen0InitialAp, scatterSize)
plt.scatter(genFinalTargetPe, genFinalMaxAps - genFinalInitialAp, scatterSize)
plt.xlabel("Target Pe [AU]")
plt.ylabel("Ap Change [AU]")
plt.grid()
plt.legend(legend)
plt.savefig(os.path.join(utils.pyplots_dir, "InO/targetPe_MaxAp.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(utils.pyplots_dir, "InO/targetPe_MaxAp.png"), bbox_inches="tight")

if plotOptimalTrajectory:

    # Plot the ideal trajectory, or maybe all trajectories in a generation or something? #
    plotterTemplateJsonPath = os.path.join(utils.InO_JsonDirPath, "InO_Plotter.json")
    dummyProblemClass = utils.InO_Problem([0,0,0,0], [0,0,0,0], simLoadTodoList=["depVarData", "propData"], templateJsonPath=plotterTemplateJsonPath)

    # Gets the index for the one with the largest Ap #
    optimalIndex = utils.findNearestInArray(fitArrayToUse[:, 0], np.amax(fitArrayToUse[:, 0]) )[-1]
    optimalFit = fitArrayToUse[optimalIndex]
    optimalIndividual = popArrayToUse[optimalIndex]

    print(optimalIndividual)
    print(optimalFit)

    # dummyProblemClassNoEDT.fitness(optimalIndividual)
    # noEDTAllSimData = dummyProblemClassNoEDT.allSimData

    dummyProblemClass.fitness(optimalIndividual)
    optimalAllSimData = dummyProblemClass.allSimDataComplete
    optimalAllSimDataStage1 = dummyProblemClass.allSimData_stage1
    optimalAllSimDataStage2 = dummyProblemClass.allSimData_stage2


    InOTrajectoryLegendLabels = ["Trajectory Stage 2", "Trajectory Stage 1", "Sun", "4", "5"]
    # utils.plotTrajectoryData(optimalAllSimData[5], sameScale=True, fignumber=10, plotOnlyTrajectory=False, saveFolder=os.path.join(utils.pyplots_dir, "InO"), savename="optimalTrajectory", savePngAndPdf=True, plotSun=True, scatterSize=250 )
    utils.plotTrajectoryData(optimalAllSimDataStage2[5], sameScale=True, fignumber=11, plotOnlyTrajectory=True )
    utils.plotTrajectoryData(optimalAllSimDataStage1[5], sameScale=True, fignumber=11, plotOnlyTrajectory=True )
    utils.plotTrajectoryData(optimalAllSimDataStage2[5], sameScale=True, fignumber=11, plotOnlyTrajectory=False,
                             saveFolder=os.path.join(utils.pyplots_dir, "InO"), savename="optimalTrajectory", savePngAndPdf=True,
                             plotSun=True, scatterSize=250, doNotPlot=True, legendLabelsCustom=InOTrajectoryLegendLabels)
    # utils.plotTrajectoryData(noEDTAllSimData[5], sameScale=False, fignumber=11)

    depVarData = optimalAllSimData[2]
    times = depVarData[:, 0]
    SMAs = depVarData[:, 1]
    es = depVarData[:, 2]
    APs = SMAs * (1 + es)

    plt.figure()
    plt.plot(2000 + times / year, APs / AU)

    plt.figure()
    plt.plot(2000 + times/year, SMAs/AU)

    plt.figure()
    plt.plot(2000 + times/year, depVarData[:, 7]/AU)

if showing:
    plt.show()


