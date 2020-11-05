import os
import numpy as np
import subprocess
import sys
from matplotlib import pyplot as plt
import scipy.interpolate as si
import json
import copy
import re

#######################################################################################################################
############################## Set various project directories ############################################
#######################################################################################################################

pyfiles_dir = os.path.dirname(os.path.realpath(__file__))
main_app_dir = os.path.abspath(os.path.join(pyfiles_dir, os.pardir))
simulation_output_dir = os.path.abspath(os.path.join(main_app_dir, "SimulationOutput"))
tudatBundle_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(main_app_dir))))
cppApplications_dir = os.path.abspath(os.path.join(tudatBundle_dir, "tudatApplications", "EDT_thesis", "tudat_apps", "bin", "applications"))
jsonInputs_dir = os.path.join(main_app_dir, "JsonInputs")
pyplots_dir = os.path.join(pyfiles_dir, "pyplots")

#######################################################################################################################
############################## Set constant values for use in GACalculatorRunner.py and others ######################
#######################################################################################################################

quickConfigsJupiter = ["Jupiter", 2020.580, 2050, 90, 398.8/365.25, False] #[planetName, start year optimal, end year, generation, synodic period, ??]
JupiterInfoList = [2020.58, 398/365.25] # Info list contains [initial same side time, synodic period]

quickConfigsSaturn = ["Saturn", 2020.580, 2050, 90, 378/365.25, False] #[planetName, start year optimal, end year, generation, synodic period, ??]
SaturnInfoList = [2020.58, 378/365.25] # Info list contains [initial same side time, synodic period]

quickConfigsMars = ["Mars", 2020.789, 2050, 10, 780/365.25, False] #[planetName, start year optimal, end year, generation, synodic period, ??]
MarsInfoList = [2020.789, 780/365.25] # Info list contains [initial same side time, synodic period]

#######################################################################################################################
############################## Set constant values for VnV stuff ######################
#######################################################################################################################
AU = 1.496e11


#######################################################################################################################
############################# Some useful functions ####################################################
#######################################################################################################################

# Function to run a generic bash command. Command contained in first entry of argumentsList, other arguments listed in
# remaining entried
def runBashCommand(argumentsList, printSetting=2): #Printsettings: 0 = no printing; 1 = print in realtime; 2 = print all at end

    # Run command in different ways, depending on print setting
    if printSetting == 1:
        p = subprocess.Popen(argumentsList, stdout=subprocess.PIPE, bufsize=1)
        for line in iter(p.stdout.readline, b''):
            decodedLine = line.decode('utf-8')
            print(decodedLine)
        p.stdout.close()
        p.wait()

    else:
        out = subprocess.Popen(argumentsList,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
        stdout,stderr = out.communicate()
        stdout = stdout.decode('utf-8')
        if printSetting == 2:
            print(stdout)
            if stderr is not None:
                print("Command finished with error: ", stderr)

def normaliseValue(initialValue, lowerBound, upperBound, denormalise = False):
    if denormalise:
        denormalisedValue = initialValue * (upperBound - lowerBound) + lowerBound
        newValue = denormalisedValue
    else:
        normalisedValue = (initialValue - lowerBound) / (upperBound - lowerBound)
        newValue = normalisedValue

    return newValue

def checkFolderExist(folderToCheck):
    """
    Function checks if the folder exists, and creates it if it doesnt
    :param folderToCheck:
    :return:
    """

    if not os.path.exists(folderToCheck):
        print("Creating new folder %s" %folderToCheck)
        os.mkdir(folderToCheck)

def findNearestInArray(array, value):
    "Element in nd array closest to the scalar value "
    idx = np.abs(array - value).argmin()
    return (array.flat[idx], idx)



#######################################################################################################################
####################################### GA Related functions #########################################################
#######################################################################################################################

# Function to get array of paretos, given an input array
def getParetoArray(costs, returnParetoArray=True, sortOutput=True):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
    """

    is_efficient = np.ones(costs.shape[0], dtype = bool)
    for i, c in enumerate(costs):
        is_efficient[i] = np.all(np.any(costs[:i]>c, axis=1)) and np.all(np.any(costs[i+1:]>c, axis=1))

    if returnParetoArray:
        paretoArray = costs[is_efficient]

        if sortOutput:
        	paretoArray = paretoArray[paretoArray[:,0].argsort()]
        return paretoArray
        
    else:
        return is_efficient

# Converts between xy coords and numpy array. If second parameter not given, is array->coords, else other way around
def arrayCoordsConvert(inputParameter1, inputParameter2=None):
	
	if inputParameter2 is None:
		outputX = inputParameter1[:,0]
		outputY = inputParameter1[:,1]
		return outputX, outputY

	else:
		outputArray = np.zeros((len(inputParameter1), 2))
		outputArray[:,0] = inputParameter1
		outputArray[:,1] = inputParameter2
		return outputArray

# Function to get lists of pareto parameters x and y, given input x and y.
def getParetoLists(XInput, YInput, sortOutput=True):
    if len(XInput) != len(YInput):
        print("ERROR: Input lengths of lists for pareto are not equal")
        sys.exit()

    inputArray = np.zeros((len(XInput), 2))
    inputArray[:,0] = XInput
    inputArray[:,1] = YInput

    paretoArray = getParetoArray(inputArray, returnParetoArray=True, sortOutput=sortOutput)

    paretoXs = paretoArray[:,0]
    paretoYs = paretoArray[:,1]

    return paretoXs, paretoYs

# Function to load and normalise the GA stuff, using bounds for normalisation
# Outputs a tuple of variables:
def loadAndNormaliseGA(GAFolder, simulation_output_dir, generation=90, fileSuffix="GA_EJ", normalising=False, DVBounds=None, TOFBounds=None):

    GA_calculatorTestTemp = os.path.join(GAFolder, "unperturbed_fullProblem_leg_0.dat")
    GA_calculatorTestPerturbedTemp = os.path.join(GAFolder, "perturbed_fullProblem_leg_0.dat")
    GA_calculatorTestDepVarsTemp = os.path.join(GAFolder, "unperturbed_depVars_leg_0.dat")


    fitnessFileDir = os.path.abspath(os.path.join(simulation_output_dir, GAFolder, "fitness_%s_%s.dat" %(fileSuffix, generation)))
    fitnessFileArray = np.genfromtxt(fitnessFileDir, delimiter=",")[:,1:3]
    fitnessFileArrayTOFs = fitnessFileArray[:,1]
    fitnessFileArrayDVs = fitnessFileArray[:,0]

    popFileDir = os.path.abspath(os.path.join(simulation_output_dir, GAFolder, "population_%s_%s.dat" %(fileSuffix, generation)))
    popFileArray = np.genfromtxt(popFileDir, delimiter=",")[:,1:3]
    launchYears = popFileArray[:,0]
    dummyList = np.ones(np.size(launchYears))


    # Get rid of bad DV cases, and combine into full array:
    fitnessFileDVsListFINAL = []
    fitnessFileTOFSListFINAL = []
    launchYearsListFINAL = []
    arrivalYearsListFINAL = []
    for j in range(len(fitnessFileArrayDVs)):
        if fitnessFileArrayDVs[j] < 10000:
            DV = fitnessFileArrayDVs[j]
            TOF = fitnessFileArrayTOFs[j]

            if normalising:
                if (DVBounds is None) or (TOFBounds is None):
                    print("ERROR: To normalise, DVBounds and TOFBounds must be specified")
                    sys.exit()
                DV = normaliseValue(DV, DVBounds[0], DVBounds[1], denormalise=True)
                TOF = normaliseValue(TOF, TOFBounds[0], TOFBounds[1], denormalise=True)

            fitnessFileDVsListFINAL.append(DV)
            fitnessFileTOFSListFINAL.append(TOF)
            launchYear = launchYears[j] + 2000
            launchYearsListFINAL.append(launchYear)
            arrivalYearsListFINAL.append(launchYear + TOF)

    return (fitnessFileDVsListFINAL, fitnessFileTOFSListFINAL, launchYearsListFINAL, arrivalYearsListFINAL)

#
# Quick configs is: [PlanetName, startYear, endYear, generation, synodicPeriod, normalising (bool)]

def getAllYearsGA(quickConfigs, GASubfolder):
    """
    Function to run through various year folders, and save to some new variables.
    :param quickConfigs: [PlanetName, startYear, endYear, generation, synodicPeriod, normalising (bool)]
    :param GASubfolderBase: Of format "GAJupiter/GACalculatorNominal_%s_%s-%s/" where %s represent planet name, start and end years
    :return outputTuple: output of GA informations (fitnessFileDVsAll, fitnessFileTOFSAll, launchYearsAll, arrivalYearsAll)
    """


    planetName = quickConfigs[0]
    startYear = quickConfigs[1]
    endYear = quickConfigs[2]
    generation = quickConfigs[3]
    synodicPeriod = quickConfigs[4]
    normalising = quickConfigs[5]

    if planetName == "Jupiter":
        fileSuffix = "GA_EJ"
    elif planetName == "Saturn":
        fileSuffix = "GA_ES"
    elif planetName == "Mars":
        fileSuffix = "GA_EM"
    else:
        print("ERROR: planetName in quickConfigs not recognised")
        sys.exit()

    # startYearsRaw = np.arange(startYear, endYear+1, synodicPeriod)
    # endYearsRaw = startYearsRaw + synodicPeriod
    # startYears = np.round(startYearsRaw, 2)
    # endYears = np.round(endYearsRaw, 2)
    #
    fitnessFileDVsAllList = []
    fitnessFileTOFSAllList = []
    launchYearsAllList = []
    arrivalYearsAllList = []
    # startEndYearSetList = [] # Differs from launchYearsAllList, as this is the start year of the band of years assessed
    startYears=[]
    endYears=[]

    subFoldersToRun = os.listdir(os.path.join(simulation_output_dir, GASubfolder))
    subFoldersToRun.sort()

    for i in range(len(subFoldersToRun)):
        # startYear = startYears[i]
        # endYear = endYears[i]

        GA_SubfolderTemp = os.path.join(GASubfolder, subFoldersToRun[i])
        years = re.findall(r'\d+.\d+', GA_SubfolderTemp)
        startYears.append(float(years[0]))
        endYears.append(float(years[1]))

        # GA_SubfolderTemp = GASubfolderBase %(planetName, "%.2f" %startYear, "%.2f" %endYear)

        normalisedGAData = loadAndNormaliseGA(GA_SubfolderTemp, simulation_output_dir, generation=generation, fileSuffix=fileSuffix, normalising=normalising)
        fitnessFileDVsListTemp, fitnessFileTOFSListTemp, launchYearsListTemp, arrivalYearsListTemp = normalisedGAData[0:4]

        fitnessFileDVsAllList.append(fitnessFileDVsListTemp)
        fitnessFileTOFSAllList.append(fitnessFileTOFSListTemp)
        launchYearsAllList.append(launchYearsListTemp)
        arrivalYearsAllList.append(arrivalYearsListTemp)
        # startEndYearSetList.append([startYear, endYear])

    fitnessFileDVsAll = np.array(fitnessFileDVsAllList)
    fitnessFileTOFSAll = np.array(fitnessFileTOFSAllList)
    launchYearsAll = np.array(launchYearsAllList)
    arrivalYearsAll = np.array(arrivalYearsAllList)
    # startEndYearSet = np.array(startEndYearSetList)

    outputTuple = (fitnessFileDVsAll, fitnessFileTOFSAll, launchYearsAll, arrivalYearsAll, startYears, endYears)
    return outputTuple

def plotManyDataGA(allYearsGAData, fignumber, quickConfigs, plotType="DV-TOFS", yearSpacing=1, markerScale=10, legendSize=11,
                   figsize=[16,9], saveFolder=None, savenameSuffix=None, scatterPointSize=1, scatterColour=None, scatterMarker=".", scatterLinewidths=None,
                   savenameOverride=None, plotLegend=True, TOFUnits="Years", removeDominated=True, plotParetoFront=False,
                   xlims=None, ylims=None, printMinDV=False):
    """
    Function to plot multiple data sets to the same plot, using all data and the spacing
    :param allYearsGAData:
    :param fignumber:
    :param plotType:
    :param yearSpacing:
    :return figLabels:
    """

    fitnessFileDVsAll, fitnessFileTOFSAll, launchYearsAll, arrivalYearsAll, startYears, endYears = allYearsGAData[0:6]
    startYearsToPlot = startYears[::yearSpacing]
    generation = quickConfigs[3]

    plt.figure(fignumber, figsize=figsize)
    figLabels=[]

    for i in range(len(fitnessFileDVsAll)):
        startYear = startYears[i]
        endYear = endYears[i]



        if TOFUnits == "Days":
            TOFsToPlot = fitnessFileTOFSAll[i] * 365.25
        else:
            TOFsToPlot = fitnessFileTOFSAll[i]

        if printMinDV:
            print(min(fitnessFileDVsAll[i]))
            print(fitnessFileDVsAll[i])
            print(TOFsToPlot)

        if startYear in startYearsToPlot:
            plt.figure(fignumber)
            if plotType=="DV-TOFS":
                xToPlot = fitnessFileDVsAll[i]
                yToPlot = TOFsToPlot
                # plt.scatter(fitnessFileDVsAll[i], TOFsToPlot, scatterPointSize, c=scatterColour, marker=scatterMarker, linewidth=scatterLinewidths)
            elif plotType=="launchYears-TOFS":
                xToPlot = launchYearsAll[i]
                yToPlot = TOFsToPlot
                # plt.scatter(launchYearsAll[i], TOFsToPlot, scatterPointSize, c=scatterColour, marker=scatterMarker, linewidth=scatterLinewidths)
            elif plotType=="launchYears-DV":
                xToPlot = launchYearsAll[i]
                yToPlot = fitnessFileDVsAll[i]
                # plt.scatter(launchYearsAll[i], fitnessFileDVsAll[i], scatterPointSize, c=scatterColour, marker=scatterMarker, linewidth=scatterLinewidths)
            else:
                print("ERROR: plotType not recognised")
                sys.exit()

            if removeDominated and plotType=="DV-TOFS":
                xyArrayRaw = arrayCoordsConvert(inputParameter1=xToPlot, inputParameter2=yToPlot)
                xyArrayPareto = getParetoArray(xyArrayRaw, returnParetoArray=True, sortOutput=True)
                xToPlot, yToPlot = xyArrayPareto[:, 0], xyArrayPareto[:, 1]
            elif removeDominated:
                xyArrayRaw = arrayCoordsConvert(inputParameter1=xToPlot, inputParameter2=yToPlot)
                DVTOFArrayRaw = arrayCoordsConvert(inputParameter1=fitnessFileDVsAll[i], inputParameter2=TOFsToPlot)

                paretoMask = getParetoArray(DVTOFArrayRaw, returnParetoArray=False)
                xyArrayPareto = xyArrayRaw[paretoMask]
                xToPlot, yToPlot = xyArrayPareto[:, 0], xyArrayPareto[:, 1]


            if plotParetoFront:
                if removeDominated == False:
                    print("ERROR: Cannot plot pareto front without removing dominated points")
                    sys.exit()

                plt.plot(xToPlot, yToPlot)

            plt.scatter(xToPlot, yToPlot, scatterPointSize, c=scatterColour, marker=scatterMarker, linewidth=scatterLinewidths)

            if xlims is not None:
                plt.xlim(xlims)
            if ylims is not None:
                plt.ylim(ylims)

            figLabels.append("%s-%s" %(startYear, endYear))

    savenameBase = "%s_%s_gen%s_space%s"
    if plotType=="DV-TOFS":
        plt.xlabel("DVs [km/s]")
        plt.ylabel("Time of flight [%s]" %TOFUnits)
        savename = savenameBase %("DV", "TOF", generation, yearSpacing)

    elif plotType=="launchYears-TOFS":
        plt.xlabel("Launch year")
        plt.ylabel("Time of flight [%s]" %TOFUnits)
        savename = savenameBase %("LaunchYears", "TOF", generation, yearSpacing)

    elif plotType=="launchYears-DV":
        plt.xlabel("Launch year")
        plt.ylabel("DVs [km/s]")
        savename = savenameBase %("LaunchYears", "DV", generation, yearSpacing)

    if savenameSuffix is not None:
        savename = savename + savenameSuffix

    if savenameOverride is not None:
        savename = savenameOverride

    savename = savename + ".png"


    if plotLegend: plt.legend(figLabels, markerscale=markerScale, prop={'size':legendSize})
    plt.grid()
    if saveFolder is not None:
        checkFolderExist(saveFolder)
        plt.savefig(os.path.join(saveFolder, savename ))



def porkchopPlot(directoryPath, baseFilename,
                 fignumber=None, figsize=[16,9],
                 contourCount=10, contourLinewidths=0.25, contourColormap=plt.cm.viridis,
                 xlabel="Departure Date", ylabel="Travel Time [%s]", cbarLabel="$\Delta$V [km/s]",
                 title="", EarthMarsCorrection=False, TOFUnit="Years",
                 saveDirectory=None, saveName=None, saveFormat="png",
                 xlims=None, ylims=None, colorbarOffset=50):




    porkchop_x_data = np.genfromtxt(os.path.join(directoryPath, baseFilename + "_x_data.dat"))
    porkchop_y_data = np.genfromtxt(os.path.join(directoryPath, baseFilename + "_y_data.dat"))
    porkchopRaw = np.genfromtxt(os.path.join(directoryPath, baseFilename + ".dat"))
    porkchop = np.transpose(porkchopRaw)
    if EarthMarsCorrection:
        porkchop_x_data = ((porkchop_x_data-2451545)/365) + 2000
    else:
        porkchop_x_data = porkchop_x_data + 2000

    ylabel = ylabel %TOFUnit
    if TOFUnit == "Days":
        porkchop_y_data = porkchop_y_data * 365.25

    plt.figure(num=fignumber, figsize=figsize)
    plt.contour(porkchop_x_data, porkchop_y_data, porkchop,
                contourCount, linewidths=contourLinewidths, colors="black")
    plt.contourf(porkchop_x_data, porkchop_y_data, porkchop, contourCount, cmap=contourColormap)
    plt.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    cbar = plt.colorbar()
    cbar.set_label(cbarLabel, rotation=0, labelpad=colorbarOffset)
    if xlims is not None:
        plt.xlim(xlims)
    if ylims is not None:
        plt.ylim(ylims)

    def fmt(x, y):
        z = np.take(si.interp2d(porkchop_x_data, porkchop_y_data, porkchop)(x, y), 0)
        return 'x={x:.5f}  y={y:.5f}  z={z:.5f}'.format(x=x, y=y, z=z)

    plt.gca().format_coord = fmt

    if (saveName is not None) and (saveDirectory is not None):
        checkFolderExist(saveDirectory)
        plt.savefig(os.path.join(saveDirectory, saveName + "." + saveFormat))

    return (porkchop_x_data, porkchop_y_data, porkchop)

def runCppGASims(jsonName, jsonDirectory=jsonInputs_dir, GAApplicationName="application_GA_calculator", printing=0):
    """
    Function to run a single simulation from json
    :param jsonName:
    :param jsonDirectory:
    :param GAApplicationName:
    :param printing:
    :return:
    """

    # Make full path to json file, and the base commands / arguments list
    jsonPath = os.path.join(jsonDirectory, jsonName)
    baseGACommand = os.path.join(cppApplications_dir, GAApplicationName)
    GAArgumentsList = [baseGACommand, jsonName]


    # Open json and use to run sim, with some info printed out
    with open(jsonPath, 'r+') as f:
        # Load json data into a variable
        GAVariables = json.load(f)

        # Set some variables to print what's being simulated
        planetName = GAVariables["PlanetConfigs"]["planetToFlyby"]

        if planetName == "Jupiter":
            planetConfigString = "JupiterConfigs"
        elif planetName == "Saturn":
            planetConfigString = "SaturnConfigs"
        elif planetName == "Mars":
            planetConfigString = "MarsConfigs"
        else:
            print("WARNING: Planet not recognised")

        startYearLower = GAVariables["PlanetConfigs"][planetConfigString]["Bounds"]["StartYearLower"]
        startYearUpper = GAVariables["PlanetConfigs"][planetConfigString]["Bounds"]["StartYearUpper"]


        # Run simulations using relevant json file
        print("Running %s simulations for %s - %s" %(planetName, startYearLower, startYearUpper))
        runBashCommand( GAArgumentsList, printSetting=printing)

def createGARunnerJsons(quickConfigs, outputSubFolderBase, jsonSaveSubDir, jsonFilenameBase, inputStartYearsRange, infoList,
                        templateJsonPath=os.path.join(jsonInputs_dir, "GAConfigsNominal.json"), createSynodicJsons=True,
                        algorithmConfigs=None):

    print("======= Creating relevant json files ==========")

    # Create and / or set some directories
    jsonSaveDir = os.path.join(jsonInputs_dir, jsonSaveSubDir)
    jsonPathBase = os.path.join(jsonSaveDir, jsonFilenameBase)
    checkFolderExist(jsonSaveDir)

    # optional inputs TODO: check which of these wanted
    integratorSettings = 1 # Add if wanted / needed
    integratorSettingsPerturbed = 1 # Add if wanted / needed
    terminationTypeEtc = 1 # add if needed / wanted



    # Set values for start years, known opposing year and synodic period
    inputStartYearInit = inputStartYearsRange[0]
    inputStartYearFin = inputStartYearsRange[-1]
    knownOpposingYear = infoList[0]
    synodicPeriod = infoList[1]

    if createSynodicJsons:
        # Create list of possible starting years, running from far in the past to far in the future
        possibleStartYears = np.arange(knownOpposingYear - 100*synodicPeriod, 3000, synodicPeriod)

        # Find value in possible years nearest to the one wanted, and create list of actual start years to run with
        value, index = findNearestInArray(possibleStartYears, inputStartYearsRange[0])
        startYearsToRun = np.round(np.arange(possibleStartYears[index], inputStartYearsRange[1], infoList[1] ), 2)
        endYearsToRun = np.round( startYearsToRun + synodicPeriod, 2)
    else:
        startYearsToRun = [inputStartYearsRange[0]]
        endYearsToRun = [inputStartYearsRange[1]]


    # Set planet config string
    if quickConfigs[0] == "Jupiter":
        planetConfigString = "JupiterConfigs"
        planetName = "Jupiter"
    elif quickConfigs[0] == "Saturn":
        planetConfigString = "SaturnConfigs"
        planetName = "Saturn"
    elif quickConfigs[0] == "Mars":
        planetConfigString = "MarsConfigs"
        planetName = "Mars"
    else:
        print("WARNING: Planet not recognised")

    # Open template json file, and use it as the structure for the rest
    with open(templateJsonPath, 'r+') as f:
        templateJsonStructure = json.load(f)
        f.close()



    for i in range(len(startYearsToRun)):
        startYear = startYearsToRun[i]
        endYear = endYearsToRun[i]
        outputSubFolder = outputSubFolderBase + "_%s_%s-%s" %(planetName, startYear, endYear)

        jsonInfoTemp = copy.deepcopy(templateJsonStructure) # Makes a true copy of the template structure

        jsonInfoTemp["PlanetConfigs"][planetConfigString]["Bounds"]["StartYearLower"] = startYear
        jsonInfoTemp["PlanetConfigs"][planetConfigString]["Bounds"]["StartYearUpper"] = endYear
        jsonInfoTemp["PlanetConfigs"]["planetToFlyby"] = planetName
        jsonInfoTemp["saveDataConfigs"]["outputSubFolder"] = outputSubFolder
        jsonInfoTemp["saveDataConfigs"]["outputSubFolderBase"] = outputSubFolderBase

        if algorithmConfigs is not None:
            jsonInfoTemp["AlgorithmConfigs"]["algorithmName"] = algorithmConfigs[0]
            jsonInfoTemp["AlgorithmConfigs"]["islandSize"] = algorithmConfigs[1]
            jsonInfoTemp["AlgorithmConfigs"]["noGenerations"] = algorithmConfigs[2]
            jsonInfoTemp["AlgorithmConfigs"]["normaliseValues"] = algorithmConfigs[3]
            jsonInfoTemp["AlgorithmConfigs"]["doGridSearch"] = algorithmConfigs[4]
            jsonInfoTemp["AlgorithmConfigs"]["gridSearchSize"] = algorithmConfigs[5]
            jsonInfoTemp["AlgorithmConfigs"]["includeDepartureDV"] = algorithmConfigs[6]
            jsonInfoTemp["AlgorithmConfigs"]["includeArrivalDV"] = algorithmConfigs[7]

        jsonPathTemp = jsonPathBase %(planetName, startYear, endYear)

        # Dump new data into a json file
        with open(jsonPathTemp, 'w+') as f:
            f.seek(0)
            json.dump(jsonInfoTemp, f, indent=4)
            f.truncate()
            f.close()

def runAllSimulations(jsonSubDirectory, jsonInputsDir=jsonInputs_dir,
                      runPath=os.path.join(cppApplications_dir, "application_GA_calculator"),
                      printSetting=0):
    """
    Runs all simulations within a json directory
    :param jsonSubDirectory:
    :param jsonInputsDir:
    :param runPath:
    :param printSetting:
    :return:
    """

    jsonsToRunFilenames = os.listdir(os.path.join(jsonInputsDir, jsonSubDirectory))
    jsonsToRunFilenames.sort()
    jsonsToRunPaths = []
    applicationName = os.path.basename(os.path.normpath(runPath))

    for i in range(len(jsonsToRunFilenames)):
        jsonsToRunPaths.append(os.path.join(jsonSubDirectory, jsonsToRunFilenames[i]))

    for i in range(len(jsonsToRunPaths)):
        # Run simulations using relevant json file
        print("Running application %s, with json %s" %(applicationName, jsonsToRunFilenames[i]))
        argumentsList = [runPath, jsonsToRunPaths[i]]
        runBashCommand( argumentsList, printSetting=printSetting)


#######################################################################################################################
####################################### VNV Related functions #########################################################
#######################################################################################################################

def getAllSimDataFromFolder(dataSubdirectory, simulationDataDirectory=simulation_output_dir, dataFilenamePortions=[]):

    fullDataFilesDirectoryPath = os.path.join(simulation_output_dir, dataSubdirectory)
    dataToLoadFilenames = os.listdir(fullDataFilesDirectoryPath)

    for i in range(len(dataToLoadFilenames)):
        filename = dataToLoadFilenames[i]

        if "bodyData" in filename:
            bodyDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            bodyDataArray = np.genfromtxt(bodyDataFilename, delimiter=",")

        elif "currentData" in filename:
            currentDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            currentDataArray = np.genfromtxt(currentDataFilename, delimiter=",")

        elif "depVarData" in filename:
            depVarDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            depVarDataArray = np.genfromtxt(depVarDataFilename, delimiter=",")

        elif "ionoData" in filename:
            ionoDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            ionoDataArray = np.genfromtxt(ionoDataFilename, delimiter=",")

        elif "magData" in filename:
            magDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            magDataArray = np.genfromtxt(magDataFilename, delimiter=",")

        elif "propData" in filename:
            propDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            propDataArray = np.genfromtxt(propDataFilename, delimiter=",")

        elif "thrustData" in filename:
            thrustDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            thrustDataArray = np.genfromtxt(thrustDataFilename, delimiter=",")

        else:
            print("Data file does not have recognised type, ignoring: %s" %filename)


    return (bodyDataArray, currentDataArray, depVarDataArray, ionoDataArray, magDataArray, propDataArray, thrustDataArray)


def plotMagData(magDataArray, bodyDataArray=None, fignumber=None, plotType="time-magnitude", logScaleY=True, logScaleX=True,
figsize=[16,9], saveFolder=None, savename=None, xlims=None, ylims=None):

    times = magDataArray[:,0]
    timesYears = (times / (365.25 * 24 * 60 * 60)) + 2000
    magnitudes = magDataArray[:, 1]
    magfieldLocal = magDataArray[:, 2:5]
    magfieldInertial = magDataArray[:, 5:8]
    theta = magDataArray[:, 8]
    B0 = magDataArray[:, 9]

    if bodyDataArray is not None:
        bodyDataTimes = bodyDataArray[:, 0]
        print(len(times))
        print(len(bodyDataTimes))
        radii = bodyDataArray[:, 1]
        radiiAU = radii / AU
        stateVector = bodyDataArray[:, 2:8]
        stateVectorAU = stateVector / AU

        if len(bodyDataTimes) != len(times):
            print("WARNING: Body data times not equal to magdata times, terminating plot")
            return 1



    plt.figure(fignumber, figsize=figsize)

    if plotType == "time-magnitude":
        xToPlot = timesYears
        yToPlot = magnitudes
        xLabel = "Year"
        yLabel = "Magnetic field magnitude [T]"

    elif plotType == "radius-magnitude":
        if bodyDataArray is None:
            print("ERROR: Body data array input required, terminating plot")
            return 1
        xToPlot = radiiAU
        yToPlot = magnitudes
        xLabel = "Radius from Sun [AU]"
        yLabel = "Magnetic field magnitude [T]"

    else:
        print("ERROR: Plot type not recgonised, ending program")
        sys.exit()

    plt.plot(xToPlot, yToPlot)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.grid()
    if logScaleY: plt.yscale("log")
    if logScaleX:
        plt.xscale("log")
        plt.xticks([10E-6, 10E-5, 10E-4, 10E-3, 10E-2, 10E-1, 10E0, 10E1, 10E2])
    if xlims is not None: plt.xlim(xlims)
    if ylims is not None: plt.ylim(ylims)

    if saveFolder is not None:
        checkFolderExist(saveFolder)
        plt.savefig(os.path.join(saveFolder, savename ))

def plotTrajectoryData(bodyDataArray, fignumber=None, plotType="x-y", legendSize=11, plotSun=False,
                       figsize=[16,9], saveFolder=None, savename=None, xlims=None, ylims=None, sameScale=False):

    times = bodyDataArray[:,0]
    timesYears = times / (365.25 * 24 * 60 * 60)
    radii = bodyDataArray[:, 1]
    stateVector = bodyDataArray[:, 2:8]
    stateVectorAU = stateVector / AU
    legendLabels=[]

    plt.figure(fignumber, figsize=figsize)

    if plotType == "x-y":
        xToPlot = stateVectorAU[:, 0]
        yToPlot = stateVectorAU[:, 1]
        xLabel = "X coordinate [AU]"
        yLabel = "Y coordinate [AU]"
        legendLabels.append("Spacecraft trajectory")
    else:
        print("ERROR: Plot type not recognised, ending program")
        sys.exit()

    if plotSun:
        plt.scatter([0], [0], c="orange")
        legendLabels.append("Sun")

    plt.plot(xToPlot, yToPlot)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.grid()
    if sameScale: plt.axis('scaled')
    plt.legend(legendLabels)
    if xlims is not None: plt.xlim(xlims)
    if ylims is not None: plt.ylim(ylims)

    if saveFolder is not None:
        checkFolderExist(saveFolder)
        plt.savefig(os.path.join(saveFolder, savename ))