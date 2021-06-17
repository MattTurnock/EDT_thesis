import os
import numpy as np
import json
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import copy
import sys
import pygmo as pg


## Create

## Structure of the variableArray: [Start Year, Initial Ap, Initial Pe]  w/ units: [Year, AU, AU]##
inputLowerBounds = [utils.launchDateRange[0],
                    utils.SSOPInitialApRangeAU[0],
                    utils.SSOPInitialPeRangeAU[0]]
inputUpperBounds = [utils.launchDateRange[1],
                    utils.SSOPInitialApRangeAU[1],
                    utils.SSOPInitialPeRangeAU[1]]


def doSSOPOptimisation(results_PopulationFileBase, results_FitnessFileBase, results_PopulationProcessedBase, results_FitnessProcessedBase, seedNumber):

    utils.checkFolderExist(utils.SSOP_Results_DirPath, emptyDirectory=False)


    # Create problem, with defined bounds #
    prob = pg.problem(utils.SSOP_Problem(inputLowerBounds, inputUpperBounds))

    # Create algorithm. Using moead with default values, and a single generation (since generations are done by evolution). Also a specific seed #
    algo = pg.algorithm(pg.de(gen=utils.SSOP_NoMoeadGenerations, seed=seedNumber))


    print("\n=== Running evolutions ===")


    bar = utils.IncrementalBar("EVOLUTION PROGRESS: ", max=utils.SSOP_NoEvolutions, suffix='[%(percent).1f%%. ETA: %(eta)ds]   ')
    for i in range(utils.SSOP_NoEvolutions):

        # Create initial population, and then continue to evolve
        if i == 0:
            pop = pg.population(prob, size=utils.SSOP_PopulationSize, seed=seedNumber)

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
        np.save(os.path.join(utils.SSOP_Results_DirPath, results_PopulationFileBase %i), populationArray)
        np.save(os.path.join(utils.SSOP_Results_DirPath, results_FitnessFileBase %i), fitnessArray)

        # Save processed arrays
        np.save(os.path.join(utils.SSOP_Results_DirPath, results_PopulationProcessedBase %i), populationArrayProcessed)
        np.save(os.path.join(utils.SSOP_Results_DirPath, results_FitnessProcessedBase %i), fitnessArrayProcessed)

        bar.next()
    bar.finish()

doSeed1 = True
doSeed2 = True

if doSeed1:
    print("seed1")
    doSSOPOptimisation(utils.SSOP_Results_PopulationFileBase, utils.SSOP_Results_FitnessFileBase, utils.SSOP_Results_PopulationProcessedBase, utils.SSOP_Results_FitnessProcessedBase, utils.SSOP_seed)

if doSeed2:
    print("seed2")
    doSSOPOptimisation(utils.SSOP_Results_PopulationFileBase_seed2, utils.SSOP_Results_FitnessFileBase_seed2, utils.SSOP_Results_PopulationProcessedBase_seed2, utils.SSOP_Results_FitnessProcessedBase_seed2, utils.SSOP_seed2)









