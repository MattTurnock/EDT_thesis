import os
import numpy as np
import json
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import copy
import sys
import pygmo as pg




## Structure of the variableArray: [Start Year, Initial Ap, Initial Pe, TargetPe]  w/ units: [Year, AU, AU, AU]##
inputLowerBounds = [utils.launchDateRange[0],
                    utils.InOInitialApRangeAU[0],
                    utils.InOInitialPeRangeAU[0],
                    utils.InOTargetPeRangeAU[0]]
inputUpperBounds = [utils.launchDateRange[1],
                    utils.InOInitialApRangeAU[1],
                    utils.InOInitialPeRangeAU[1],
                    utils.InOTargetPeRangeAU[1]]



def doInOOptimisation(results_PopulationFileBase, results_FitnessFileBase, results_PopulationProcessedBase, results_FitnessProcessedBase, seedNumber):

    utils.checkFolderExist(utils.InO_Results_DirPath, emptyDirectory=False)
    # utils.checkFolderExist(utils.InODir, emptyDirectory=True)


    # Create problem, with defined bounds #
    prob = pg.problem(utils.InO_Problem(inputLowerBounds, inputUpperBounds, printSetting=utils.InO_printSetting))

    # Create algorithm. Using moead with default values, and a single generation (since generations are done by evolution). Also a specific seed #
    algo = pg.algorithm(pg.de(gen=utils.InO_NoMoeadGenerations, seed=seedNumber))


    print("\n=== Running evolutions ===")


    bar = utils.IncrementalBar("EVOLUTION PROGRESS: ", max=utils.InO_NoEvolutions, suffix='[%(percent).1f%%. ETA: %(eta)ds]   ')
    for i in range(utils.InO_NoEvolutions):

        # Create initial population, and then continue to evolve
        if i == 0:
            pop = pg.population(prob, size=utils.InO_PopulationSize, seed=seedNumber)

        else:
            pop = algo.evolve(pop)

        # Get raw population and fitness arrays
        populationArray =  pop.get_x()
        fitnessArray = pop.get_f()

        # Process arraays (basically just make the Ap fitness correct
        populationArrayProcessed = populationArray
        fitnessArrayProcessed = fitnessArray
        fitnessArrayProcessed[:, 0] = 1/fitnessArray[:, 0]

        # Save raw arrays
        np.save(os.path.join(utils.InO_Results_DirPath, results_PopulationFileBase %i), populationArray)
        np.save(os.path.join(utils.InO_Results_DirPath, results_FitnessFileBase %i), fitnessArray)

        # Save processed arrays
        np.save(os.path.join(utils.InO_Results_DirPath, results_PopulationProcessedBase %i), populationArrayProcessed)
        np.save(os.path.join(utils.InO_Results_DirPath, results_FitnessProcessedBase %i), fitnessArrayProcessed)

        bar.next()
    bar.finish()


doSeed1 = True
doSeed2 = False

if doSeed1:
    print("seed1")
    doInOOptimisation(utils.InO_Results_PopulationFileBase, utils.InO_Results_FitnessFileBase, utils.InO_Results_PopulationProcessedBase, utils.InO_Results_FitnessProcessedBase, utils.InO_seed)

if doSeed2:
    print("seed2")
    doInOOptimisation(utils.InO_Results_PopulationFileBase_seed2, utils.InO_Results_FitnessFileBase_seed2, utils.InO_Results_PopulationProcessedBase_seed2, utils.InO_Results_FitnessProcessedBase_seed2, utils.InO_seed2)












