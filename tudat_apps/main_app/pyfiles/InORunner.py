import os
import numpy as np
import json
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import copy
import sys
import pygmo as pg




## Structure of the variableArray: [Start Year, Initial Ap, Initial Pe, TargetPe]  w/ units: [Year, AU, AU, AU]##
inputLowerBounds = [utils.launchDateRange[0],
                    utils.SSOPInitialApRangeAU[0],
                    utils.SSOPInitialPeRangeAU[0]]
inputUpperBounds = [utils.launchDateRange[1],
                    utils.SSOPInitialApRangeAU[1],
                    utils.SSOPInitialPeRangeAU[1]]

SSOPChangeKeys =  [ ["GuidanceConfigs", "initialEphemerisYear"],
                    ["GuidanceConfigs", "terminationSettings", "timeTerminationYears"],
                    ["GuidanceConfigs", "vehicleInitialKeplerian", "a_au"],
                    ["GuidanceConfigs", "vehicleInitialKeplerian", "e"],
                    ["GuidanceConfigs", "vehicleInitialKeplerian", "i_deg"],
                    ["GuidanceConfigs", "vehicleInitialKeplerian", "aop_deg"],
                    ["GuidanceConfigs", "vehicleInitialKeplerian", "raan_deg"],
                    ["GuidanceConfigs", "vehicleInitialKeplerian", "ta_deg"],
                    ["saveDataConfigs", "outputSubFolder"],
                    ["saveDataConfigs", "baseFilename"]]

# SSOP_JsonDirPathBase = os.path.join(utils.jsonInputs_dir, "SSOP", "SSOP_Gen_%s")
# SSOP_JsonNameBase = "SSOP_Gen_%s_Pop_%s.json"
#
# SSOP_SimOutputSubFolderBase = "SSOP/SSOP_Gen_%s/SSOP_Gen_%s_Pop_%s/"
# SSOP_SimOutputFilenameBase = "SSOP_Gen_%s_Pop_%s-"

utils.checkFolderExist(utils.SSOP_Results_DirPath, emptyDirectory=True)

class SSOP_Problem:

    def __init__(self, inputLowerBounds, inputUpperBounds):
        ### Do any problem initialisation steps ###
        self.lowerBounds = inputLowerBounds
        self.upperBounds = inputUpperBounds



    def fitness(self, x):

        # Get variables from allowed ones #
        launchEpoch = x[0] * utils.year
        launchYear = x[0]
        initialApAU = x[1]
        initialPeAU = x[2]

        # Calculate intermediate Keplerian parameters to pass to json #
        a, e, i, AOP, RAAN, TA = utils.ApPeToKeplerian(initialApAU, initialPeAU)

        # Set SSOP change values, and create the new json #
        SSOPChangeValues = [ launchYear,
                             utils.SSOP_simulationLengthYears,
                             a,
                             e,
                             i,
                             AOP,
                             RAAN,
                             TA,
                             utils.SSOP_SimOutputSubFolderTemp,
                             utils.SSOP_SimOutputFilenameTemp]

        utils.createModifiedJson(utils.SSOP_JsonDirPathBaseFile, utils.SSOP_JsonDirPathTemp, utils.SSOP_JsonNameTemp, SSOPChangeKeys, SSOPChangeValues)

        # Run the simulation using the temporary json #
        utils.runAllSimulations(utils.SSOP_JsonSubDirPathTemp, runPath=utils.simulationRunPathDefault, printSetting=utils.SSOP_printSetting, runOnlyThisFile=utils.SSOP_JsonNameTemp, printProgress=False)

        # Load simulation data from file #
        allSimData = utils.getAllSimDataFromJson(os.path.join(utils.SSOP_JsonDirPathTemp, utils.SSOP_JsonNameTemp), todoList=["propData", "depVarData"], printInfo=False)

        # Use simulation data to get the Apogee and Perigee at the target TOF #
        maxAp, maxPe = utils.getApPe_At_TOF(allSimData, utils.SSOP_targetTOFYears)
        maxApAU = maxAp / utils.AU

        # Use known orbit information to find DV to put into that orbit, for second fitness parameter #
        VCircEarth = utils.getOrbitalV(utils.mu_Sun, 1 * utils.AU, 1*utils.AU)
        VKick = utils.getOrbitalV(utils.mu_Sun, 1*utils.AU, a*utils.AU)
        DVKick = abs(VKick - VCircEarth)

        return [1/maxApAU, DVKick]


    # Return number of objectives
    def get_nobj(self):
        return 2

    def get_bounds(self):
        ### Gets the simulation bounds for each input parameter ###
        return (self.lowerBounds, self.upperBounds)


    def get_name(self):
        return "SSO+ Problem"


    def get_extra_info(self):
        return "\tHuman readable lower bounds: [%s, %s, %s] " \
               "\n\tHuman readable upper bounds: [%s, %s, %s]" %(self.lowerBounds[0], self.lowerBounds[1], self.lowerBounds[2],
                                                                 self.upperBounds[0], self.upperBounds[1], self.upperBounds[2])

# Create problem, with defined bounds #
prob = pg.problem(SSOP_Problem(inputLowerBounds, inputUpperBounds))

# Create algorithm. Using moead with default values, and a single generation (since generations are done by evolution). Also a specific seed #
algo = pg.algorithm(pg.moead(gen=utils.SSOP_NoMoeadGenerations, seed=utils.SSOP_seed, neighbours=utils.SSOP_MoeadNeighbours))
# print(prob)
# print(algo)
# print("\n=== Creating initial population ===")
# Create and run initial population, using given population size. Also a specific seed #


print("\n=== Running evolutions ===")


bar = utils.IncrementalBar("EVOLUTION PROGRESS: ", max=utils.SSOP_NoEvolutions, suffix='[%(percent).1f%%. ETA: %(eta)ds]   ')
for i in range(utils.SSOP_NoEvolutions):

    # Create initial population, and then continue to evolve
    if i == 0:
        pop = pg.population(prob, size=utils.SSOP_PopulationSize, seed=utils.SSOP_seed)

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
    np.save(os.path.join(utils.SSOP_Results_DirPath, utils.SSOP_Results_PopulationFileBase %i), populationArray)
    np.save(os.path.join(utils.SSOP_Results_DirPath, utils.SSOP_Results_FitnessFileBase %i), fitnessArray)

    # Save processed arrays
    np.save(os.path.join(utils.SSOP_Results_DirPath, utils.SSOP_Results_PopulationProcessedBase %i), populationArrayProcessed)
    np.save(os.path.join(utils.SSOP_Results_DirPath, utils.SSOP_Results_FitnessProcessedBase %i), fitnessArrayProcessed)

    bar.next()
bar.finish()











