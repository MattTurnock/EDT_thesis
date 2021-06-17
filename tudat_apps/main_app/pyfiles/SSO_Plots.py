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

showing = False
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
fitArrayListSeed2 = []
popArrayListSeed2 = []
for i in range(len(numpyArraysFileList_Processed)):
    filePath = numpyArraysFileList_Processed[i]

    if "Seed2" not in filePath:
        if "Fit" in filePath:
            fitArrayList.append(np.load(filePath, allow_pickle=True))

        elif "Pop" in filePath:
            popArrayList.append(np.load(filePath, allow_pickle=True))
    else:
        if "Fit" in filePath:
            fitArrayListSeed2.append(np.load(filePath, allow_pickle=True))

        elif "Pop" in filePath:
            popArrayListSeed2.append(np.load(filePath, allow_pickle=True))

# print(fitArrayList[-1])
plt.figure()

if generationToUse != "all":

    fitArrayToUse = fitArrayList[generationToUse]
    popArrayToUse = popArrayList[generationToUse]

    fitArrayToUseSeed2 = fitArrayListSeed2[generationToUse]
    popArrayToUseSeed2 = popArrayListSeed2[generationToUse]
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

    # Plot the ideal trajectory #
    dummyProblemClass = utils.SSOP_Problem([0,0,0], [0,0,0], templateJsonPath=utils.SSOP_TemplateJsonPathPlotter)
    # dummyProblemClassNoEDT = utils.SSOP_Problem([0,0,0], [0,0,0], templateJsonPath=utils.SSOP_TemplateJsonPathPlotter_NoEDT)

    # Gets the index for the one with the largest Ap #
    optimalIndex = utils.findNearestInArray(fitArrayToUse[:, 0], np.amax(fitArrayToUse[:, 0]) )[-1]
    optimalFit = fitArrayToUse[optimalIndex]
    optimalIndividual = popArrayToUse[optimalIndex]
    print(fitArrayToUse)

    # print(optimalIndividual)
    # print(optimalFit)
    print("---SSOP Results, Main Seed----")
    print("Optimal Launch Year: %s" %optimalIndividual[0])
    print("Optimal Initial Ap: %s" %optimalIndividual[1])
    print("Optimal Initial Pe: %s" %optimalIndividual[2])
    print("Highest Max Ap: %s\n" %optimalFit[0])


    # newfit = 1/dummyProblemClass.fitness(optimalIndividual)[0]
    # noEDTAllSimData = dummyProblemClassNoEDT.allSimData

    newfit = 1/dummyProblemClass.fitness(optimalIndividual)[0]
    optimalAllSimData = dummyProblemClass.allSimData
    print("NEWFIT: ", newfit)



    utils.plotTrajectoryData(optimalAllSimData[5], sameScale=True, fignumber=10, plotOnlyTrajectory=False, saveFolder=os.path.join(utils.pyplots_dir, "SSOP"), savename="optimalTrajectory", savePngAndPdf=True, plotSun=True, scatterSize=500 )
    # utils.plotTrajectoryData(noEDTAllSimData[5], sameScale=False, fignumber=11)

    depVarData = optimalAllSimData[2]
    times = depVarData[:, 0]
    SMAs = depVarData[:, 1]
    es = depVarData[:, 2]
    APs = SMAs * (1 + es)

    # plt.figure()
    # plt.plot(2000 + times / year, APs / AU)
    #
    # plt.figure()
    # plt.plot(2000 + times/year, SMAs/AU)
    #
    # plt.figure()
    # plt.plot(2000 + times/year, depVarData[:, 7]/AU)




    # Get data from seed2 #
    dummyProblemClass_Seed2 = utils.SSOP_Problem([0,0,0], [0,0,0], templateJsonPath=utils.SSOP_TemplateJsonPathPlotter)
    # dummyProblemClassNoEDT_Seed2 = utils.SSOP_Problem([0,0,0], [0,0,0], templateJsonPath=utils.SSOP_TemplateJsonPathPlotter_NoEDT)
    #
    # Gets the index for the one with the largest Ap #
    optimalIndex_Seed2 = utils.findNearestInArray(fitArrayToUseSeed2[:, 0], np.amax(fitArrayToUseSeed2[:, 0]) )[-1]
    optimalFit_Seed2 = fitArrayToUseSeed2[optimalIndex_Seed2]
    optimalIndividual_Seed2 = popArrayToUseSeed2[optimalIndex_Seed2]

    # print(optimalIndividual)
    # print(optimalFit)
    print("---SSOP Results, Seed2 ----")
    print("Optimal Launch Year: %s" %optimalIndividual_Seed2[0])
    print("Optimal Initial Ap: %s" %optimalIndividual_Seed2[1])
    print("Optimal Initial Pe: %s" %optimalIndividual_Seed2[2])
    print("Highest Max Ap: %s\n" %optimalFit_Seed2[0])




    # Plot the ideal trajectory, with perturbations on #
    dummyProblemClass_Perturbed = utils.SSOP_Problem([0,0,0], [0,0,0], templateJsonPath=utils.SSOP_TemplateJsonPathPlotter_Perturbed)

    # Gets the index for the one with the largest Ap #
    optimalIndex_Perturbed = utils.findNearestInArray(fitArrayToUse[:, 0], np.amax(fitArrayToUse[:, 0]) )[-1]
    # optimalFit_Perturbed = fitArrayToUse[optimalIndex_Perturbed]
    optimalIndividual_Perturbed = popArrayToUse[optimalIndex_Perturbed]

    # print(optimalIndividual)
    # print(optimalFit)





    optimalFit_Perturbed = dummyProblemClass_Perturbed.fitness(optimalIndividual_Perturbed)
    optimalAllSimData_Perturbed = dummyProblemClass.allSimData

    print("---SSOP Results, Perturbed----")
    print("Optimal Launch Year: %s" %optimalIndividual_Perturbed[0])
    print("Optimal Initial Ap: %s" %optimalIndividual_Perturbed[1])
    print("Optimal Initial Pe: %s" %optimalIndividual_Perturbed[2])
    print("Highest Max Ap: %s\n" %(1/optimalFit_Perturbed[0])  )


    utils.plotTrajectoryData(optimalAllSimData_Perturbed[5], sameScale=True, fignumber=11, plotOnlyTrajectory=False, saveFolder=os.path.join(utils.pyplots_dir, "SSOP"), savename="optimalTrajectory_Perturbed", savePngAndPdf=True, plotSun=True, scatterSize=500 )



if showing:
    plt.show()


