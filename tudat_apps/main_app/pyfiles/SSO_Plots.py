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
plotFolder = "pyplots/SSOP"
utils.checkFolderExist(plotFolder)

showing = True
useExample = False
plotOptimalTrajectory = True

generationToUse = -1

if useExample:
    results_DirPathToUse = os.path.join(utils.numpyBinary_dir, "SSOP_EXAMPLE")
else:
    results_DirPathToUse = utils.SSOP_Results_DirPath

print("-- Loading Numpy Array --")
numpyArraysFileList = natsort.natsorted(glob.glob(results_DirPathToUse + "/*"))
numpyArraysFileList_Processed = natsort.natsorted(glob.glob(results_DirPathToUse + "/*Processed*"))
# print(numpyArraysFileList_Processed)

# SSOP_Populations
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
gen0MaxAps = np.reshape(fitArrayList[0], np.shape(gen0InitialAp))

genFinalStartYears = popArrayList[-1][:, 0]
genFinalInitialAp = popArrayList[-1][:, 1]
genFinalMaxAps = np.reshape(fitArrayList[-1], np.shape(genFinalInitialAp))

scatterSize = 3

legend = ["First generation", "Final generation"]

plt.figure(figsize=utils.figSizeDefault)
plt.scatter(gen0StartYears, gen0MaxAps, scatterSize)
plt.scatter(genFinalStartYears, genFinalMaxAps, scatterSize)
plt.xlabel("Start Year")
plt.ylabel("Maximum Ap [AU]")
plt.grid()
plt.legend(legend)
plt.savefig(os.path.join(utils.pyplots_dir, "SSOP/Year_MaxAp.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(utils.pyplots_dir, "SSOP/Year_MaxAp.png"), bbox_inches="tight")

plt.figure(figsize=utils.figSizeDefault)
# print(gen0InitialAp)
# print("")
# print(gen0MaxAps[:,0])

plt.scatter(gen0InitialAp, gen0MaxAps - gen0InitialAp, scatterSize)
plt.scatter(genFinalInitialAp, genFinalMaxAps - genFinalInitialAp, scatterSize)
plt.xlabel("Initial Ap [AU]")
plt.ylabel("Ap Change [AU]")
plt.grid()
plt.legend(legend)
plt.savefig(os.path.join(utils.pyplots_dir, "SSOP/InitialAp_MaxAp.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(utils.pyplots_dir, "SSOP/InitialAp_MaxAp.png"), bbox_inches="tight")

if plotOptimalTrajectory:

    # Plot the ideal trajectory, or maybe all trajectories in a generation or something? #

    dummyProblemClass = utils.SSOP_Problem([0,0,0], [0,0,0], templateJsonPath=utils.SSOP_TemplateJsonPathPlotter)
    dummyProblemClassNoEDT = utils.SSOP_Problem([0,0,0], [0,0,0], templateJsonPath=utils.SSOP_TemplateJsonPathPlotter_NoEDT)

    # Gets the index for the one with the largest Ap #
    optimalIndex = utils.findNearestInArray(fitArrayToUse[:, 0], np.amax(fitArrayToUse[:, 0]) )[-1]
    optimalFit = fitArrayToUse[optimalIndex]
    optimalIndividual = popArrayToUse[optimalIndex]

    print(optimalIndividual)
    print(optimalFit)

    dummyProblemClassNoEDT.fitness(optimalIndividual)
    noEDTAllSimData = dummyProblemClassNoEDT.allSimData

    dummyProblemClass.fitness(optimalIndividual)
    optimalAllSimData = dummyProblemClass.allSimData



    utils.plotTrajectoryData(optimalAllSimData[5], sameScale=True, fignumber=10, plotOnlyTrajectory=False, saveFolder=os.path.join(utils.pyplots_dir, "SSOP"), savename="optimalTrajectory", savePngAndPdf=True, plotSun=True )
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


