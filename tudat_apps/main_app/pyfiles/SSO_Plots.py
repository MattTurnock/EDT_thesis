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
useExample = True
plotOptimalTrajectory = False

generationToUse = "all"

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
    plt.scatter(fitArrayToUse[:, 1]/ 1000, fitArrayToUse[:, 0], 1 )

else:
    for i in range(len(fitArrayList)):
        plt.scatter(fitArrayList[i][:, 1] / 1000, fitArrayList[i][:, 0], 1)
plt.xlabel("Kick DV [km/s]")
plt.ylabel("Maximum Ap [AU]")


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

    dummyProblemClass.fitness(optimalIndividual)
    optimalAllSimData = dummyProblemClass.allSimData

    dummyProblemClassNoEDT.fitness(optimalIndividual)
    noEDTAllSimData = dummyProblemClassNoEDT.allSimData

    utils.plotTrajectoryData(optimalAllSimData[5], sameScale=True, fignumber=10, plotOnlyTrajectory=True )
    utils.plotTrajectoryData(noEDTAllSimData[5], fignumber=11)

if showing:
    plt.show()


