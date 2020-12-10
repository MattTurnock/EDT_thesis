import os
import numpy as np
import subprocess
import sys
from matplotlib import pyplot as plt
import scipy.interpolate as si
import json
import copy
import re
from scipy.interpolate import interp1d
import pandas as pd
from datetime import datetime
import logging

#######################################################################################################################
############################## Set various project directories ############################################
#######################################################################################################################

pyfiles_dir = os.path.dirname(os.path.realpath(__file__))
main_app_dir = os.path.abspath(os.path.join(pyfiles_dir, os.pardir))
tudat_apps_dir = os.path.abspath(os.path.join(main_app_dir, os.pardir))
pyplots_dir = os.path.join(pyfiles_dir, "pyplots")

simulation_output_dir = os.path.abspath(os.path.join(main_app_dir, "SimulationOutput"))
tudatBundle_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(main_app_dir))))
cppApplications_dir = os.path.abspath(os.path.join(tudatBundle_dir, "tudatApplications", "EDT_thesis", "tudat_apps", "bin", "applications"))
jsonInputs_dir = os.path.join(main_app_dir, "JsonInputs")

deepSpaceMagfieldDataDir = os.path.join(tudat_apps_dir, "deepSpaceMagfieldData")
voyagerMagDataSubDirBase = "Voyagers/pub/data/voyager/voyager%s/magnetic_fields/ip_1hour_ascii"
voyagerMagDataSubDirBase_48s = "Voyagers/pub/data/voyager/voyager%s/magnetic_fields/VIM_48s_mag_ascii"
voyagerTrajDataSubDirBase = "Voyagers/pub/data/voyager/voyager%s/plasma/sedr"
voyager1MagDataDir = os.path.join(deepSpaceMagfieldDataDir, voyagerMagDataSubDirBase % "1")
voyager2MagDataDir = os.path.join(deepSpaceMagfieldDataDir, voyagerMagDataSubDirBase % "2")
voyager1MagDataDir_48s = os.path.join(deepSpaceMagfieldDataDir, voyagerMagDataSubDirBase_48s % "1")
voyager2MagDataDir_48s = os.path.join(deepSpaceMagfieldDataDir, voyagerMagDataSubDirBase_48s % "2")
voyager1TrajDataDir = os.path.join(deepSpaceMagfieldDataDir, voyagerTrajDataSubDirBase % "1")
voyager2TrajDataDir = os.path.join(deepSpaceMagfieldDataDir, voyagerTrajDataSubDirBase % "2")

pythonRunnerLogfilesDir = os.path.join(pyfiles_dir, "cppRunnerLogfiles")


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
c = 299792458
W0 = 1367.0

indicatorString_9004 = "SN      DecimalYear     F1Mag   ElevAng AzimAng F2Mag"
indicatorString_0919 = " SC	Yr		DOY	  F1		   Br		   Bt		   Bn"
referenceParkerRsAU = [10,
                       11.1727129303638,
                       12.4017913168933,
                       13.6765745494801,
                       15.0823930695153,
                       16.7415642533113,
                       18.462433831912,
                       20.3601920250856,
                       22.4530212577845,
                       24.7609729309714,
                       27.3061595341309,
                       30.112966505075,
                       33.2082858668696,
                       36.6217738803607,
                       40.386135180867,
                       44.8288990982471,
                       49.4368728534518,
                       54.5185013839415,
                       60.1224717825836,
                       66.3024757007037,
                       73.5962225610273,
                       80.6335128892091,
                       88.9218518548644,
                       98.7038919812137,
                       108.141984321039,
                       119.257925948224,
                       131.516477996652,
                       145.984231852466,
                       160.989987572133,
                       177.5381886769]
referenceParkerBsnT = [0.523960135300264,
                       0.473629338502,
                       0.428133239871939,
                       0.387007425813149,
                       0.349832093577504,
                       0.316227766016838,
                       0.291683780824842,
                       0.263665089873036,
                       0.238337830856296,
                       0.215443469003188,
                       0.194748303990876,
                       0.176041084386555,
                       0.159130851241951,
                       0.143844988828766,
                       0.132680474971472,
                       0.119935394620923,
                       0.108414586893583,
                       0.09800045006277,
                       0.090394155316908,
                       0.081711033154572,
                       0.073861998220794,
                       0.066766929391876,
                       0.061584821106603,
                       0.054555947811685,
                       0.049315388199506,
                       0.045487779470038,
                       0.040296113202004,
                       0.036425331154496,
                       0.033598182862838,
                       0.031622776601684]
referenceParkerDataArray = np.zeros( (len(referenceParkerBsnT), 2) )
referenceParkerDataArray[:,0] = referenceParkerRsAU
referenceParkerDataArray[:,1] = referenceParkerBsnT


#######################################################################################################################
############################# Some useful functions ####################################################
#######################################################################################################################

# Function to get true anomaly from a, e, r (NOTE: hyperbolic only, i think). Returns in degrees
def getTA(a, e, r):
    inside = ( a*(1-e**2) - r )/(e*r)
    TA = 1/(np.cos(inside))

    return np.rad2deg(TA)

# a = -1000
# e = 1.3
# r = 1
#
# print(getTA(a, e, r))

# Function to run a generic bash command. Command contained in first entry of argumentsList, other arguments listed in
# remaining entried
def runBashCommand(argumentsList, printSetting=2): #Printsettings: 0 = no printing; 1 = print in realtime; 2 = print all at end

    # Run command in different ways, depending on print setting
    if printSetting == 1:
        p = subprocess.Popen(argumentsList, stdout=subprocess.PIPE, bufsize=1)
        for line in iter(p.stdout.readline, b''):
            decodedLine = line.decode('utf-8')
            decodedLine = decodedLine.strip("\n")
            logger.info(decodedLine)
        p.stdout.close()
        p.wait()

    else:
        out = subprocess.Popen(argumentsList,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
        stdout,stderr = out.communicate()
        stdout = stdout.decode('utf-8')
        if printSetting == 2:
            logger.info(stdout)
            if stderr is not None:
                logger.info("Command finished with error: ", stderr)

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
        logger.info("Creating new folder %s" %folderToCheck)
        os.mkdir(folderToCheck)

def findNearestInArray(array, value):
    "Element in nd array closest to the scalar value "
    idx = np.abs(array - value).argmin()
    return (array.flat[idx], idx)

# Lists all files in a directory, with a specific extension
def list_files(directory, extension=None, sort=True, exceptions=None):

    fileListFull = os.listdir(directory)

    if extension is None:
        fileList = fileListFull

    else:
        fileList = []
        for file in fileListFull:
            if file.endswith(".%s" %extension):
                fileList.append(file)

    if sort:
        fileList.sort()

    if exceptions is not None:
        newFileList = []
        for file in fileList:
            addFile = True
            for exception in exceptions:
                if exception in file:
                    addFile = False

            if addFile:
                newFileList.append(file)

        fileList = newFileList

    return fileList

def print_full(x):
    pd.set_option('display.max_rows', len(x))
    logger.info(x)
    pd.reset_option('display.max_rows')

    # Function to reshape arrays by interpolation. inputArray is the base data, reframeArray is the array of indices to reshape along. Must be on 0 column!
def interpolateArrays(inputArray, reframeArray):



    inputDataFrame = pd.DataFrame(inputArray[:,1:], index=inputArray[:,0])
    outputDataFrame = inputDataFrame[~inputDataFrame.index.duplicated(keep='first')]

    outputDataFrame = outputDataFrame.reindex(outputDataFrame.index.values.tolist()+list(reframeArray))

    outputDataFrame.sort_index(inplace=True)

    outputDataFrame.interpolate(inplace=True)

    outputDataFrame = inputDataFrame[~inputDataFrame.index.duplicated(keep='first')]
    outputDataFrame.sort_index(inplace=True)

    # print_full(outputDataFrame)
    outputDataFrame = outputDataFrame.reindex(reframeArray, method="nearest")
    



    outputIndices = outputDataFrame.index.to_numpy()
    outputData = outputDataFrame.to_numpy()

    outputArray = np.zeros(( len(outputIndices), np.shape(outputData)[1]+1 ))

    outputArray[:,0] = outputIndices
    outputArray[:,1:] = outputData

    return outputArray




# "saveDataConfigs": {
#     "outputSubFolder": "SSO-CHB-Test-custom-2/",
#     "baseFilename": "SSO-CH-Test-out-E6-"
# },
# "scConfigs": {
#     "SRP": {
#         "endmassArea1": 1.2,
#         "endmassArea2": 1.2,
#         "endmassRadiationCoefficient": 0.5,
#         "tetherRadiationCoefficient": 0.5,
#         "rotationFactor": 0.7
#     }
# },



def createModifiedJson(jsonTemplatePath, jsonSaveToDir, jsonSaveToFilename, listOfChangeKeys, listOfChangeValues):

    # Create save folder if not existing, and set path of save file
    checkFolderExist(jsonSaveToDir)
    jsonSaveToPath = os.path.join(jsonSaveToDir, jsonSaveToFilename)

    # Open template json file, and use it as the structure for the rest
    with open(jsonTemplatePath, 'r+') as f:
        templateJson = json.load(f)
        f.close()

    # Make copy of template json, to use as new json
    newJson = copy.deepcopy(templateJson)

    # Check if lists are of same length
    if len(listOfChangeKeys) != len(listOfChangeValues):
        logger.info("ERROR: Length of keys list must equal length of values list, skipping changes")
        return 1

    # Make the relevant changes to each key here
    for i in range(len(listOfChangeKeys)):
        
        # If-else statements to change the value. Works up to a key nesting of 5
        tempKeys = listOfChangeKeys[i]
        keyNesting = len(tempKeys)
        if keyNesting == 1:
            newJson [tempKeys[0]] = listOfChangeValues[i]
        elif keyNesting == 2:
            newJson [tempKeys[0]] [tempKeys[1]] = listOfChangeValues[i]
        elif keyNesting == 3:
            newJson [tempKeys[0]] [tempKeys[1]] [tempKeys[2]] = listOfChangeValues[i]
        elif keyNesting == 4:
            newJson [tempKeys[0]] [tempKeys[1]] [tempKeys[2]] [tempKeys[3]] = listOfChangeValues[i]
        elif keyNesting == 5:
            newJson [tempKeys[0]] [tempKeys[1]] [tempKeys[2]] [tempKeys[3]] [tempKeys[4]] = listOfChangeValues[i]
        else:
            logger.info("ERROR: Function only supports key nesting up to 5, skipping changes")
            return 1

    # Dump new data into a json file
    with open(jsonSaveToPath, 'w+') as f:
        f.seek(0)
        json.dump(newJson, f, indent=4)
        f.truncate()
        f.close()

    return 0

#Function to create name for a logfile using the time
def createLogfileName(baseName):

    timeNow = datetime.now()
    timeNowString = timeNow.strftime("%d-%m-%Y-%H-%M-%S")
    logfileName = "%s-%s.log" %(baseName, timeNowString)

    return logfileName

# Function to create the logger for a file, must replace print statements with logger.info statements

def createLogger(logfileName, makeDummy=False):
    logger = logging.getLogger('scope.name')
    logger.handlers=[]

    if makeDummy is False:
        file_log_handler = logging.FileHandler(logfileName)
        logger.addHandler(file_log_handler)

    stderr_log_handler = logging.StreamHandler()
    logger.addHandler(stderr_log_handler)
    logger.setLevel('DEBUG')

    # nice output format
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    if makeDummy is False: file_log_handler.setFormatter(formatter)
    stderr_log_handler.setFormatter(formatter)

    return logger

# Dummy log variable for when none made
logger = createLogger("", makeDummy=True)


# testVariablesJsonPath = "/home/matt/LinkToEDT_thesis/tudat_apps/main_app/JsonInputs/testVariables.json"
# saveDir = "/home/matt/LinkToEDT_thesis/tudat_apps/main_app/JsonInputs/testDir/"
# saveName = "test.json"
# exampleChangeKeys = [["saveDataConfigs","outputSubFolder"], ["scConfigs", "SRP", "endmassArea1"]]
# exampleChangeValues = ["badabingBadaboom/", 999.1]
#
# createModifiedJson(testVariablesJsonPath, saveDir, saveName, exampleChangeKeys, exampleChangeValues)

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
        logger.info("ERROR: Input lengths of lists for pareto are not equal")
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
                    logger.info("ERROR: To normalise, DVBounds and TOFBounds must be specified")
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
        logger.info("ERROR: planetName in quickConfigs not recognised")
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
            logger.info(min(fitnessFileDVsAll[i]))
            logger.info(fitnessFileDVsAll[i])
            logger.info(TOFsToPlot)

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
                logger.info("ERROR: plotType not recognised")
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
                    logger.info("ERROR: Cannot plot pareto front without removing dominated points")
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
            logger.info("WARNING: Planet not recognised")

        startYearLower = GAVariables["PlanetConfigs"][planetConfigString]["Bounds"]["StartYearLower"]
        startYearUpper = GAVariables["PlanetConfigs"][planetConfigString]["Bounds"]["StartYearUpper"]


        # Run simulations using relevant json file
        logger.info("Running %s simulations for %s - %s" %(planetName, startYearLower, startYearUpper))
        runBashCommand( GAArgumentsList, printSetting=printing)

def createGARunnerJsons(quickConfigs, outputSubFolderBase, jsonSaveSubDir, jsonFilenameBase, inputStartYearsRange, infoList,
                        templateJsonPath=os.path.join(jsonInputs_dir, "GAConfigsNominal.json"), createSynodicJsons=True,
                        algorithmConfigs=None):



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
        logger.info("WARNING: Planet not recognised")

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

        return 0

def runAllSimulations(jsonSubDirectory, jsonInputsDir=jsonInputs_dir,
                      runPath=os.path.join(cppApplications_dir, "application_GA_calculator"),
                      printSetting=0, fileIgnores=[]):
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
        jsonFilenameTemp = jsonsToRunFilenames[i]
        toAppend = True

        for ignoreString in fileIgnores:
            if ignoreString in jsonFilenameTemp:
                toAppend = False

        if toAppend:
            jsonsToRunPaths.append(os.path.join(jsonSubDirectory, jsonFilenameTemp))
        else:
            logger.info("Ignoring file: %s" %jsonFilenameTemp )


    for i in range(len(jsonsToRunPaths)):
        # Run simulations using relevant json file
        logger.info("\n\nRunning application %s, with json %s\n" %(applicationName, jsonsToRunPaths[i].split("/")[-1]))
        argumentsList = [runPath, jsonsToRunPaths[i]]
        runBashCommand( argumentsList, printSetting=printSetting)


#######################################################################################################################
####################################### VNV Related functions #########################################################
#######################################################################################################################

def getAllSimDataFromFolder(dataSubdirectory, simulationDataDirectory=simulation_output_dir, dataFilenamePortions=[],
                            todoList=["bodyData", "currentData", "depVarData", "ionoData", "magData", "propData", "thrustData"]):

    fullDataFilesDirectoryPath = os.path.join(simulation_output_dir, dataSubdirectory)
    dataToLoadFilenames = os.listdir(fullDataFilesDirectoryPath)

    # Make dummy variables in case not wanted to load
    bodyDataArray = np.array([0,0,0])
    currentDataArray = np.array([0,0,0])
    depVarDataArray = np.array([0,0,0])
    ionoDataArray = np.array([0,0,0])
    magDataArray = np.array([0,0,0])
    propDataArray = np.array([0,0,0])
    thrustDataArray = np.array([0,0,0])

    for i in range(len(dataToLoadFilenames)):
        filename = dataToLoadFilenames[i]

        if ("bodyData" in filename) and ("bodyData" in todoList):
            bodyDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            bodyDataArray = np.genfromtxt(bodyDataFilename, delimiter=",")

        elif ("currentData" in filename) and ("currentData" in todoList):
            currentDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            currentDataArray = np.genfromtxt(currentDataFilename, delimiter=",")

        elif ("depVarData" in filename) and ("depVarData" in todoList):
            depVarDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            depVarDataArray = np.genfromtxt(depVarDataFilename, delimiter=",")

        elif ("ionoData" in filename) and ("ionoData" in todoList):
            ionoDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            ionoDataArray = np.genfromtxt(ionoDataFilename, delimiter=",")

        elif ("magData" in filename) and ("magData" in todoList):
            magDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            magDataArray = np.genfromtxt(magDataFilename, delimiter=",")

        elif ("propData" in filename) and ("propData" in todoList):
            propDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            propDataArray = np.genfromtxt(propDataFilename, delimiter=",")

        elif ("thrustData" in filename) and ("thrustData" in todoList):
            thrustDataFilename = os.path.join(fullDataFilesDirectoryPath, filename)
            thrustDataArray = np.genfromtxt(thrustDataFilename, delimiter=",")

        else:
            logger.info("Data file does not have recognised type (or is purposely being ignored), ignoring: %s" %filename)


    return (bodyDataArray, currentDataArray, depVarDataArray, ionoDataArray, magDataArray, propDataArray, thrustDataArray)


def plotMagData(magDataArray, arrayType="simulation", bodyDataArray=None, fignumber=None, plotType="time-magnitude", logScaleY=True, logScaleX=True,
figsize=[16,9], saveFolder=None, savename=None, xlims=None, ylims=None, scatter=False, scatterPointSize=1, legend=None, gridOn=False, color=None):

    if arrayType == "simulation":
        times = magDataArray[:,0]
        timesYears = (times / (365.25 * 24 * 60 * 60)) + 2000
        magnitudes = magDataArray[:, 1]
        magfieldLocal = magDataArray[:, 2:5]
        magfieldInertial = magDataArray[:, 5:8]
        theta = magDataArray[:, 8]
        B0 = magDataArray[:, 9]


        if bodyDataArray is not None:
            bodyDataTimes = bodyDataArray[:, 0]
            radii = bodyDataArray[:, 1]
            radiiAU = radii / AU
            stateVector = bodyDataArray[:, 2:8]
            stateVectorAU = stateVector / AU

            if len(bodyDataTimes) != len(times):
                logger.info("WARNING: Body data times not equal to magdata times, terminating plot")
                return 1

    elif arrayType == "voyager":
        timesYears = magDataArray[:,0]
        radiiAU = magDataArray[:, 1]
        magnitudes = magDataArray[:,2]
        deltasDeg = magDataArray[:,3]
        lambdasDeg = magDataArray[:,4]
        BRs = magDataArray[:, 5]
        BTs = magDataArray[:, 6]
        BNs = magDataArray[:, 7]

    elif arrayType == "simple-radius-magnitude":
        radiiAU = magDataArray[:,0]
        magnitudes = magDataArray[:,1]

    plt.figure(fignumber, figsize=figsize)

    if plotType == "time-magnitude":
        xToPlot = timesYears
        yToPlot = magnitudes
        xLabel = "Year"
        yLabel = "Magnetic field magnitude [nT]"

    elif plotType == "radius-magnitude":
        xToPlot = radiiAU
        yToPlot = magnitudes
        xLabel = "Radius from Sun [AU]"
        yLabel = "Magnetic field magnitude [nT]"

    elif plotType == "radius-BR":
        xToPlot = radiiAU
        yToPlot = BRs
        xLabel = "Radius from Sun [AU]"
        yLabel = "Magnetic field component strength [nT]"
    elif plotType == "radius-BT":
        xToPlot = radiiAU
        yToPlot = BTs
        xLabel = "Radius from Sun [AU]"
        yLabel = "Magnetic field component strength [nT]"
    elif plotType == "radius-BN":
        xToPlot = radiiAU
        yToPlot = BNs
        xLabel = "Radius from Sun [AU]"
        yLabel = "Magnetic field component strength [nT]"

    else:
        logger.info("ERROR: Plot type not recgonised, ending program")
        sys.exit()

    if scatter:
        plt.scatter(xToPlot, yToPlot, scatterPointSize, c=color)
    else:
        plt.plot(xToPlot, yToPlot, color=color)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    if gridOn: plt.grid(which="both")
    if logScaleY: plt.yscale("log")
    if logScaleX:
        plt.xscale("log")
        plt.xticks([10E-6, 10E-5, 10E-4, 10E-3, 10E-2, 10E-1, 10E0, 10E1, 10E2])
    if xlims is not None: plt.xlim(xlims)
    if ylims is not None: plt.ylim(ylims)
    if legend is not None: plt.legend(legend)

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
        logger.info("ERROR: Plot type not recognised, ending program")
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

# Function to load magfield data from voyager datafiles. load type is based on year coverage: can be "7789" or "9004"
def loadVoyagerMagfieldData(dirPath, loadType="7789"):

    if (loadType == "7789") :

        fileList = list_files(dirPath, extension="asc", sort=True, exceptions=None)
        outputArray = np.genfromtxt(os.path.join(dirPath, fileList[0]))
        outputArray = np.delete(outputArray, np.where (outputArray==0) [0], 0)

    elif (loadType == "9004") or (loadType == "0919")or (loadType=="traj"):

        if loadType == "9004":
            indicatorString = indicatorString_9004
            extension = "txt"
            exceptions = ["00readme", "format"]
        elif loadType == "0919":
            indicatorString = indicatorString_0919
            extension = "dat"
            exceptions = ["V2_48s_2019_001-241_sans-masked", "V2_48s_2015_001-031", "V2_48s_2012_301-331", "V2_48s_2011_091-121_141024_corr",
                          "V2_48s_2009_301-331", "V2_48s_2009_091-121", "V2_48s_2009_061-091"]
        elif loadType == "traj":
            indicatorString = None
            extension = "asc"
            exceptions = ["v2_sedr_1977_2030"]

        fileList = list_files(dirPath, extension=extension, sort=True, exceptions=exceptions)
        fileDataList = []
        for file in fileList:
            filePath = os.path.join(dirPath, file)

            if indicatorString is not None:
                with open(filePath) as search:
                    for num, line in enumerate(search, 1):
                        # line = search[lineNo]
                        line = line.rstrip()  # remove '\n' at end of line
                        if line == indicatorString:
                            skipHeaderNum = num
            else:
                skipHeaderNum=0

            fileData = np.genfromtxt( filePath, skip_header=skipHeaderNum )
            fileData = np.delete(fileData, np.where (fileData==999) [0], 0)

            # If statement checks for odd files without SN col
            if (fileData[0,0] != 1.0) and (fileData[0,0] != 2.0) and (loadType=="9004"):
                ReplacementSNCol = np.ones((np.shape(fileData)[0], 1))
                fileData = np.concatenate(( ReplacementSNCol, fileData ), axis=1)


            fileDataList.append(fileData)

        for i in range(len(fileDataList)):
            if i==0:
                fullDataArray = fileDataList[i]
            else:
                fullDataArray = np.concatenate((fullDataArray, fileDataList[i]), axis=0)

        outputArray = fullDataArray

    return outputArray

# Function to average data by year
def averageMagDataArray(inputArray, interval=1):

    startTime = inputArray[0,0]
    endTime = inputArray[-2,0]
    newTimes = np.arange(startTime, endTime, interval)
    # outputArray = np.zeros( ( len(newTimes), np.shape(inputArray)[1] ) )
    # outputArray[:,0] = newTimes

    intervalMeans = []
    for i in range(len(newTimes)-1):
        intervalSubData = inputArray[ np.where(np.logical_and(inputArray[:,0] >= newTimes[i], inputArray[:,0] < newTimes[i+1])) ]
        intervalMean = np.mean(intervalSubData, axis=0)
        intervalMeans.append(intervalMean)

    outputArray= np.array(intervalMeans)
    return outputArray


# Function converts decimal year to regular year (eg 99.2 --> 1999.2). FOr one-dimensional arrays
def decimalYearArray2YearArray(inputYear, limit=50):
    # np.delete(fileData, np.where (fileData==999) [0], 0)
    pre20s = 1900 + inputYear[ np.where(inputYear >= limit)]
    post20s = 2000 + inputYear[ np.where(inputYear <= limit)]
    outputArray = np.concatenate(( pre20s, post20s ))
    return outputArray


# Function to convert array of format [year, day, hour] to simple years. Day type has Jan1 or Jan0, depending if Jan 1st is counted as 1 or 0
def yearDayHour2Years(inputArray, yearType="19", dayType="Jan1", hoursType="nominal"):

    if yearType == "19":
        years = 1900 + inputArray[:,0]
    elif yearType == "20":
        years = 2000 + inputArray[:,0]
    elif yearType == "YYYY":
        years = inputArray[:,0]


    if dayType == "Jan1":
        days = (inputArray[:,1] - 1) / 365.25
    elif dayType == "Jan0":
        days = inputArray[:,1] / 365.25
    elif dayType is None:
        days = np.zeros(len(inputArray))

    if hoursType == "nominal":
        hours = inputArray[:,2] / (365.25 * 24)
    elif hoursType is None:
        hours = np.zeros(len(inputArray))

    return years + days + hours



# FUnction to convert regfular trajectory array to anew one of format [Year, Radius]
def voyagerTrajArrayToRefined(inputArray):

    ## OLD FILE TYPE LOADING ##
    # outputArray = np.zeros( ( np.shape(inputArray)[0], 2 ) )
    #
    # times = yearDayHour2Years(inputArray[:,0:2], yearType="YYYY", dayType="Jan1", hoursType=None)
    #
    # outputArray[:, 0] = times
    # outputArray[:, 1] = inputArray[:,2]

    outputArray = np.zeros( ( np.shape(inputArray)[0], 2 ) )

    times = inputArray[:, 15] - 100 + 2000
    radii = inputArray[:, 12]

    outputArray[:,0] = times
    outputArray[:,1] = radii

    return outputArray

def voyagerArrayToPlotArray(inputArray, inputType="7789", trajectoryDataArray=None):

    outputArray = np.zeros( (np.shape(inputArray)[0], 8)   )


    if inputType == "7789":

        times = yearDayHour2Years(inputArray[:, 1:4], yearType="19", dayType="Jan1")

        outputArray[:,0] = times # Times column
        outputArray[:,1] = inputArray[:, 7] # Radii column
        outputArray[:, 2] = inputArray[:, 9] # Magfield column (F2)
        outputArray[:, 3] = inputArray[:, 10] # Delta column
        outputArray[:, 4] = inputArray[:, 11] # Lambda column

    elif (inputType == "9004") or (inputType == "0919"):
        if trajectoryDataArray is None:
            logger.info("ERROR: trajectory data array required for 1990-2004 data, ending plotting")
            return 1


        if inputType == "9004":
            times = decimalYearArray2YearArray(inputArray[:, 1])
        elif inputType == "0919":
            times = yearDayHour2Years(inputArray[:, 1:3], yearType="20", dayType="Jan1", hoursType=None)

        trajectoryDataArrayRefined = voyagerTrajArrayToRefined(trajectoryDataArray)
        trajectoryDataArrayRefinedInterpolated = interpolateArrays(trajectoryDataArrayRefined, times)

        if inputType == "9004":
            outputArray[:, 0] = times # times column
            outputArray[:, 1] = trajectoryDataArrayRefinedInterpolated[:,1] # Radii column
            outputArray[:, 2] = inputArray[:, 5] # Magfield column (F2)
            outputArray[:, 3] = inputArray[:, 3] # Delta column
            outputArray[:, 4] = inputArray[:, 4] # Lambda column
        elif inputType == "0919":
            outputArray[:, 0] = times # times column
            outputArray[:, 1] = trajectoryDataArrayRefinedInterpolated[:,1] # Radii column
            outputArray[:, 2] = inputArray[:, 3] # Magfield column (F1, F2 unavailable)

            # NOTE: no deltas or lambdas for these later data
            # outputArray[:, 3] = inputArray[:, 3] # Delta column
            # outputArray[:, 4] = inputArray[:, 4] # Lambda column

    elif inputType == "0919":
        if trajectoryDataArray is None:
            logger.info("ERROR: trajectory data array required for 1990-2004 data, ending plotting")
            return 1



    # Calculate components of B and add to the array. Formulae taken from voyager data itself
    BComponentsArray = np.zeros( (np.shape(outputArray)[0], 3) )
    F2sArray = outputArray[:, 2]
    deltasArray = outputArray[:, 3]
    lambdasArray = outputArray[:, 4]

    BComponentsArray[:, 0] = F2sArray * np.cos(np.deg2rad(lambdasArray)) * np.cos(np.deg2rad(deltasArray)) # Br
    BComponentsArray[:, 1] = F2sArray * np.sin(np.deg2rad(lambdasArray)) * np.cos(np.deg2rad(deltasArray)) # Bt
    BComponentsArray[:, 2] = F2sArray * np.sin(np.deg2rad(deltasArray)) # Bn

    outputArray[:, 5:8] = BComponentsArray


    return outputArray







