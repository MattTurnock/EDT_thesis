import numpy as np
import os
from matplotlib import pyplot as plt
import matplotlib.animation as animation

AU = 1.496e11 # AU in m

# Set various project directories
pyfiles_dir = os.path.dirname(os.path.realpath(__file__))
main_app_dir = os.path.abspath(os.path.join(pyfiles_dir, os.pardir))
simulation_output_dir = os.path.abspath(os.path.join(main_app_dir, "SimulationOutput"))

def normaliseValue(initialValue, lowerBound, upperBound, denormalise = False):
    if denormalise:
        denormalisedValue = initialValue * (upperBound - lowerBound) + lowerBound
        newValue = denormalisedValue
    else:
        normalisedValue = (initialValue - lowerBound) / (upperBound - lowerBound)
        newValue = normalisedValue

    return newValue





################### TESTING OUTPUT FROM SSO-CHB TEST 1 #####################

"SSO-CHB-Test/SSO-CH-Test-in-0-000325.dat"
year = 365*24*60*60
AU = 1.496E11
def testPlot(filename, parName, parVal, scale=10, plotType = "propData", filename2="", filename3="", savefig=False, savename=""):



    testDataFileDir = os.path.abspath(os.path.join(simulation_output_dir, filename))
    testDataArray = np.genfromtxt(testDataFileDir, delimiter=",")
    times = testDataArray[:,0]/year

    if filename2 != "":
        testDataFileDir2 = os.path.abspath(os.path.join(simulation_output_dir, filename2))
        testDataArray2 = np.genfromtxt(testDataFileDir2, delimiter=",")
        times2 = testDataArray2[:,0]/year

    if filename3 != "":
        testDataFileDir3 = os.path.abspath(os.path.join(simulation_output_dir, filename3))
        testDataArray3 = np.genfromtxt(testDataFileDir3, delimiter=",")
        times3 = testDataArray3[:,0]/year

    if plotType == "propData":

        # scale = 10
        legendList = []

        # Create figure, title, scale etc
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        plt.title("x-y plot, 10 year, %s = %s" %(parName, parVal))
        plt.xlabel("x coord [AU]")
        plt.ylabel("y coord [AU]")
        plt.xlim([-scale, scale])
        plt.ylim([-scale, scale])
        plt.grid()

        # Add Earth orbit circle
        an = np.linspace(0, 2 * np.pi, 100)
        plt.plot(1*np.cos(an), 1*np.sin(an))
        legendList.append("Earth Orbit")
        plt.plot(0.4*np.cos(an), 0.4*np.sin(an))
        legendList.append("Mercury Orbit")
        plt.plot(4*np.cos(an), 4*np.sin(an))
        legendList.append("Jupiter Orbit")
        plt.plot(70*np.cos(an), 70*np.sin(an))
        legendList.append("Termination Shock")

        # Add spacecraft plot
        plt.plot(testDataArray[:, 1]/AU, testDataArray[:, 2]/AU)
        legendList.append("Spacecraft Orbit")

        # Add legend
        plt.legend(legendList)

        plt.savefig("pyplots/%s-%s-%s.png" %(parName, parVal, scale))

    elif plotType == "posTheta":
        legendList = []

        # Create figure, title, scale etc
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        plt.title("posTheta plot")
        plt.xlabel("time (years)")
        plt.ylabel("theta (deg)")
        plt.grid()

        # Add theta, cos(theta) and sin(theta) plots
        plt.plot(testDataArray[:, 0]/year, np.rad2deg(testDataArray[:, 2]))
        legendList.append("Theta")
        # plt.plot(testDataArray[:, 0]/year, np.cos(testDataArray[:, 2]))
        # legendList.append("cos(theta)")
        # plt.plot(testDataArray[:, 0]/year, np.sin(testDataArray[:, 2]))
        # legendList.append("sin(theta)")

        # Add legend
        plt.legend(legendList)


    elif plotType == "magField":
        legendList = []

        # Create figure, title, scale etc
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        plt.title("Magnetic field strength plot")
        plt.xlabel("time (years)")
        plt.ylabel("Magfield strength (T)")
        plt.grid()

        # Add direct saved, local and inertial magfield plots

        # Do direct data
        plt.plot(testDataArray[:, 0]/year, testDataArray[:, 1])
        legendList.append("B (local, direct)")


        # Calculate local strength
        localStrength = np.zeros((len(testDataArray[:, 2:]), 2))
        localStrength[:,0] = testDataArray[:, 0]
        for i in range(len(localStrength)):
            localStrength[i,1] = np.linalg.norm(testDataArray[i, 3:6])

        plt.plot(testDataArray[:, 0]/year, localStrength[:,1])
        legendList.append("B (local)")

        # Calculate inertial strength
        inertialStrength = np.zeros((len(testDataArray[:, 6:]), 2))
        inertialStrength[:,0] = testDataArray[:, 0]
        for i in range(len(inertialStrength)):
            inertialStrength[i,1] = np.linalg.norm(testDataArray[i, 6:])

        plt.plot(testDataArray[:, 0]/year, inertialStrength[:,1])
        legendList.append("B (inertial)")

        # Add legend
        plt.legend(legendList)


    elif plotType == "magField2":
        legendList = []

        # Create figure, title, scale etc
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        plt.title("Magnetic field local and inertial plot")
        plt.xlabel("time (years)")
        plt.ylabel("Magfield strength (T)")
        plt.grid()

        # Add local and inertial magfield plots
        plt.plot(testDataArray[:, 0]/year, testDataArray[:, 3])
        legendList.append("Bx (local)")

        plt.plot(testDataArray[:, 0]/year, testDataArray[:, 4])
        legendList.append("By (local)")

        plt.plot(testDataArray[:, 0]/year, testDataArray[:, 5])
        legendList.append("Bz (local)")

        plt.plot(testDataArray[:, 0]/year, testDataArray[:, 6])
        legendList.append("Bx (inertial)")

        plt.plot(testDataArray[:, 0]/year, testDataArray[:, 7])
        legendList.append("By (inertial)")

        plt.plot(testDataArray[:, 0]/year, testDataArray[:, 8])
        legendList.append("Bz (inertial)")

        # Add legend
        plt.legend(legendList)

    elif plotType == "magField3":
        legendList = []

        # Create figure, title, scale etc
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        plt.title("Magnetic field B0 plot")
        plt.xlabel("time (years)")
        plt.ylabel("Magfield strength (T)")
        plt.grid()

        # Add local and inertial magfield plots
        plt.plot(testDataArray[:, 0]/year, testDataArray[:, 9])
        legendList.append("B0")

        # Add legend
        plt.legend(legendList)

    elif plotType == "keplerPlot1":
        legendList = []

        # Create figure, title, scale etc
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        plt.title("Kepler plot 1 (Pe, Ap, SMA)")
        plt.xlabel("time (years)")
        plt.ylabel("Distance (AU)")
        plt.grid()

        # Get kepler element lists

        SMAs = testDataArray[:,1]/AU
        ECCs = testDataArray[:,2]
        INCs = testDataArray[:,3]
        AOPs = testDataArray[:,4]
        RAANs = testDataArray[:,5]
        TAs = testDataArray[:,6]

        # Calculate other elements from base kepler elements
        Aps = SMAs * (1+ECCs)
        Pes = SMAs * (1-ECCs)

        # Plot SMA, Apo and Peri
        plt.plot(times, SMAs)
        legendList.append("SMA")

        plt.plot(times, Aps)
        legendList.append("Ap")

        plt.plot(times, Pes)
        legendList.append("Pe")

        # Add legend
        plt.legend(legendList)


    elif plotType == "altitudePlot":
        legendList = []

        # Create figure, title, scale etc
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        plt.title("Altitude plot")
        plt.xlabel("time (years)")
        plt.ylabel("Distance (AU)")
        plt.grid()

        # Create altitude list (from dependent variable data)
        positionArray = testDataArray[:, 1:4]
        altitudes = np.zeros(len(testDataArray[:,0]))
        # print(altitudes)
        for i in range(len(testDataArray[:,0])):
            altitudes[i] = np.linalg.norm(positionArray[i])
        altitudes = altitudes/AU

        # Plot altitudes
        plt.plot(times, altitudes)
        legendList.append("alt")

        # plt.plot(times2, testDataArray2[:, 1]/AU)
        # legendList.append("R")

        # Add legend
        plt.legend(legendList)

    elif plotType == "thrustPlot":
        legendList = []

        # Create figure, title, scale etc
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        plt.title("Thrust magniutde plot")
        plt.xlabel("time (years)")
        plt.ylabel("Thrust (N)")
        plt.grid()

        # Create thrust magnitude lists (for local and inertial)
        inertialThrustArray = testDataArray[:, 2:5]
        inertialThrustMagnitudes = np.zeros(len(testDataArray[:,0]))
        for i in range(len(testDataArray[:,0])):
            inertialThrustMagnitudes[i] = np.linalg.norm(inertialThrustArray[i])

        localThrustArray = testDataArray[:, 5:]
        localThrustMagnitudes = np.zeros(len(testDataArray[:,0]))
        for i in range(len(testDataArray[:,0])):
            localThrustMagnitudes[i] = np.linalg.norm(localThrustArray[i])


        # Plot thrusts
        plt.plot(times, inertialThrustMagnitudes)
        legendList.append("Inertial thrust magnitude")

        plt.plot(times, localThrustMagnitudes)
        legendList.append("Local thrust magnitude")

        # Add legend
        plt.legend(legendList)

    elif plotType == "statePlot":
        legendList = []

        # Create figure, title, scale etc
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        plt.title("State plot")
        plt.xlabel("time (years)")
        plt.ylabel("Distance (AU)")
        plt.grid()

        # Plot bodyData and integrationResult states
        plt.plot(times, testDataArray[:,2]/AU)
        legendList.append("bodyData x")
        plt.plot(times, testDataArray[:,3]/AU)
        legendList.append("bodyData y")
        plt.plot(times, testDataArray[:,4]/AU)
        legendList.append("bodyData z")

        # plt.plot(times2, testDataArray2[:,1]/AU)
        # legendList.append("propData x")
        # plt.plot(times2, testDataArray2[:,2]/AU)
        # legendList.append("propData y")
        # plt.plot(times2, testDataArray2[:,3]/AU)
        # legendList.append("propData z")


        # Add legend
        plt.legend(legendList)

    elif plotType == "propDataGA":

        # scale = 10
        legendList = []

        # Create figure, title, scale etc
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        # plt.title("x-y plot, Earth to Planet")
        plt.xlabel("x coord [AU]")
        plt.ylabel("y coord [AU]")
        plt.xlim([-scale, scale])
        plt.ylim([-scale, scale])
        plt.grid()

        # Add Earth orbit circle
        an = np.linspace(0, 2 * np.pi, 100)
        plt.plot(1*np.cos(an), 1*np.sin(an), color='blue')
        # plt.plot(testDataArray3[:,1]/AU, testDataArray3[:,2]/AU, color='blue')
        legendList.append("Earth Orbit")
        # plt.plot(0.4*np.cos(an), 0.4*np.sin(an))
        # legendList.append("Mercury Orbit")
        plt.plot(5.2*np.cos(an), 5.2*np.sin(an), color='green')
        # plt.plot(testDataArray3[:,4]/AU, testDataArray3[:,5]/AU, color='green')
        legendList.append("Jupiter Orbit")
        # plt.plot(9.58*np.cos(an), 9.58*np.sin(an))
        # legendList.append("Saturn Orbit")
        # plt.plot(70*np.cos(an), 70*np.sin(an))
        # legendList.append("Termination Shock")

        depVarsFirstRow = testDataArray3[0,:]
        depVarsFinalRow = testDataArray3[-1,:]
        EarthStart = [depVarsFirstRow[1], depVarsFirstRow[2]]
        PlanetEnd = [depVarsFinalRow[4], depVarsFinalRow[5]]

        plt.plot(EarthStart[0]/AU, EarthStart[1]/AU, 'ro', color='blue') # TODO: Make me modular
        legendList.append("Earth at start")
        # plt.plot(78934550712.8347/AU, -778090576795.287/AU, 'ro') # TODO: Make me modular
        # legendList.append("Earth at end")
        # plt.plot(-38693.254000109/AU, -4496.40733001221/AU, 'ro') # TODO: Make me modular
        # legendList.append("Jupiter at start")
        plt.plot(PlanetEnd[0]/AU, PlanetEnd[1]/AU, 'ro', color='green') # TODO: Make me modular
        legendList.append("Jupiter at end")
        plt.plot(0,0, 'ro', color='orange')
        legendList.append('Sun')

        # Add spacecraft plot
        plt.plot(testDataArray[:, 1]/AU, testDataArray[:, 2]/AU)
        legendList.append("Lambert Trajectory")

        if filename2 != "":
            plt.plot(testDataArray2[:, 1]/AU, testDataArray2[:, 2]/AU)
            legendList.append("Perturbed Trajectory")

        # Add legend
        plt.legend(legendList)

        initialEpoch = testDataArray[0,0]
        finalEpoch = testDataArray[-1,0]
        TOFSseconds = finalEpoch - initialEpoch
        print("TOF: %s" %(TOFSseconds/60/60/24/365.25))
        print("Launch date (since 2000): %s" %(initialEpoch/60/60/24/365.25))

        # plt.savefig("pyplots/%s-%s-%s.png" %(parName, parVal, scale))


    # elif plotType == "AnimatedPC":
    #     #NOTE: uses filename as the main trajectory file, and depvars as filename2
    #     x_earth = testDataArray2[:,1]
    #     y_earth = testDataArray2[:,2]
    #     x_target = testDataArray2[:,4]
    #     y_target = testDataArray2[:,5]
    #
    #     x = testDataArray[:,1]/AU
    #     y = testDataArray[:,2]/AU
    #
    #     fig, ax = plt.subplots()
    #     line, = ax.plot(x, y, color='k')
    #
    #     def update(num, x, y, line):
    #         line.set_data(x[:num], y[:num])
    #         # line.axes.axis([0, 10, 0, 1])
    #         return line,
    #
    #     ani = animation.FuncAnimation(fig, update, len(x), fargs=[x, y, line],
    #                                   interval=10, blit=False)


    if savefig == True:
        if savename == "": plotSaveName = plotType
        else: plotSaveName = savename
        plt.savefig("pyplots/%s.png" %(plotSaveName))


# testPlot("SSO-CHB-Test/SSO-CH-Test-out-0-0.dat", "0ms2",1.5)
# testPlot("SSO-CHB-Test/SSO-CH-Test-out-0-000325.dat", "0-000325ms2",75)
# testPlot("SSO-CHB-Test/SSO-CH-Test-out-0-000325.dat", "0-000325ms2",8)
# testPlot("SSO-CHB-Test/SSO-CH-Test-in-0-000325.dat", "-0-000325ms2",1.1)
# testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-0-000325.dat", "0-000325ms2",75)
# testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-expo-0-000325.dat", "0-000325ms2",75)

# testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-expo-0-000325.dat", "constAccel", "0-000325ms2",scale=75)
# testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-expo-enviro__E-7.dat", "constCurrent", "10^-7 T",scale=4)
# testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-expo-enviro__E-6.dat", "constCurrent", "10^-6 T",scale=10)

baseDirectE6 = "SSO-CHB-Test-custom-2/SSO-CH-Test-out-E6-%s.dat"

# testPlot(baseDirectE6 %"propData", "const current", "subheading", plotType="propData")
# testPlot(baseDirectE6 %"magData", "Theta over time", "subheading", plotType="posTheta")
# testPlot(baseDirectE6 %"magData", "Magfield over time", "subheading", plotType="magField")
# testPlot(baseDirectE6 %"magData", "Magfield over time", "subheading", plotType="magField2")
# testPlot(baseDirectE6 %"depVarData", "Kepler plot", "subheading", plotType="keplerPlot1")
# testPlot(baseDirectE6 %"depVarData", "Altitude plot", "subheading", plotType="altitudePlot", filename2=baseDirectE6 %"bodyData")
# testPlot(baseDirectE6 %"thrustData", "Thrust magnitude plot", "subheading", plotType="thrustPlot")
# testPlot(baseDirectE6 %"bodyData", "States plot", "subheading", plotType="statePlot", filename2=baseDirectE6 %"propData")
# testPlot(baseDirectE6 %"magData", "Magfield B0 over time", "subheading", plotType="magField3")
GA_Subfolder = "GA_calculator_2020-2050/"
plotFolder = "pyplots/" + GA_Subfolder
if not os.path.exists(plotFolder):
    print("Creating new folder")
    os.mkdir(plotFolder)

GA_calculatorTest = GA_Subfolder + "unperturbed_fullProblem_leg_0.dat"
GA_calculatorTestPerturbed = GA_Subfolder + "perturbed_fullProblem_leg_0.dat"
GA_calculatorTestDepVars = GA_Subfolder + "unperturbed_depVars_leg_0.dat"
# GA_calculatorDepVars = "GA_calculator/fullProblemInterplanetaryTrajectoryDepVars_0_leg_0.dat"

# SOme new testplots
# testPlot(GA_calculatorTest, "", "", plotType="AnimatedPC", filename2=GA_calculatorDepVars)
testPlot(GA_calculatorTest, "", "", plotType="propDataGA", savefig=True, savename="GA_test",
         filename2=GA_calculatorTestPerturbed, filename3=GA_calculatorTestDepVars, scale=6)
# testPlot(GA_calculatorTestPerturbed, "", "", plotType="propDataGA", savefig=True, savename="GA_test")

normalising = True
fitnessFile0Dir = os.path.abspath(os.path.join(simulation_output_dir, GA_Subfolder + "fitness_GA_EJ_0.dat"))
fitnessFile0Array = np.genfromtxt(fitnessFile0Dir, delimiter=",")[:,1:3]
fitnessFile0ArrayTOFs = fitnessFile0Array[:,1]
minTOF0 = min(fitnessFile0ArrayTOFs)
avgTOF0 = np.mean(fitnessFile0ArrayTOFs)

theIs = [10,20,30,40,50,60,70,80,90]
for i in theIs:
    fitnessFile9Dir = os.path.abspath(os.path.join(simulation_output_dir, GA_Subfolder + "fitness_GA_EJ_%s.dat" %i))
    fitnessFile9Array = np.genfromtxt(fitnessFile9Dir, delimiter=",")[:,1:3]
    fitnessFile9ArrayTOFs = fitnessFile9Array[:,1]
    fitnessFile9ArrayDVs = fitnessFile9Array[:,0]
    minTOF9 = min(fitnessFile9ArrayTOFs)
    avgTOF9 = np.mean(fitnessFile9ArrayTOFs)
    avgDV9 = np.mean(fitnessFile9ArrayDVs)
    # print("Fitness 0 min: %s" %(minTOF0/year))
    # print("Fitness 0 mean: %s" %(avgTOF0/year))
    # print("Fitness 9 min: %s" %(minTOF9/year))
    # print("Fitness 9 mean: %s" %(avgTOF9/year))
    # print("Fitness 9 mean DV: %s" %(avgDV9))

    popFile9Dir = os.path.abspath(os.path.join(simulation_output_dir, GA_Subfolder + "population_GA_EJ_%s.dat" %i))
    popFile9Array = np.genfromtxt(popFile9Dir, delimiter=",")[:,1:3]
    launchYears = popFile9Array[:,0]
    dummyList = np.ones(np.size(launchYears))


    # Get rid of bad DV cases, and combine into full array:
    fitnessFile9DVsListFINAL = []
    fitnessFile9TOFSListFINAL = []
    launchYearsListFINAL = []
    arrivalYearsListFINAL = []
    for j in range(len(fitnessFile9ArrayDVs)):
        if fitnessFile9ArrayDVs[j] < 10000:
            DV = fitnessFile9ArrayDVs[j]
            TOF = fitnessFile9ArrayTOFs[j]

            if normalising: #TODO: Modularise me!
                DV = normaliseValue(DV, 8.5, 9.68, denormalise=True)
                TOF = normaliseValue(TOF, 2.5, 10, denormalise=True)


            fitnessFile9DVsListFINAL.append(DV)
            fitnessFile9TOFSListFINAL.append(TOF)
            launchYear = launchYears[j] + 2000
            launchYearsListFINAL.append(launchYear)
            arrivalYearsListFINAL.append(launchYear + TOF)



    JupiterAps = np.arange(2005.267, 2060, 11.9)


    plt.figure()
    plt.scatter(launchYearsListFINAL, fitnessFile9TOFSListFINAL, 1)
    plt.xlabel("Launch year")
    plt.ylabel("Time of flight (years)")
    # plt.ylim([0,5])
    plt.xlim([2010, 2060])
    plt.grid()
    plt.savefig(plotFolder + "TOF_launch_10k_%s.png" %i)

    plt.figure()
    plt.scatter(fitnessFile9DVsListFINAL, fitnessFile9TOFSListFINAL, 1)
    plt.xlabel("DVs (km/s)")
    plt.ylabel("Time of flight (years)")
    # plt.ylim([0,5])
    plt.grid()
    plt.savefig(plotFolder + "TOF_DV_10k_%s.png" %i)

    # plt.figure()
    # plt.scatter(arrivalYearsListFINAL, fitnessFile9TOFSListFINAL, 1)
    # plt.xlabel("Arrival year")
    # plt.ylabel("Time of flight (years)")
    # for i in JupiterAps:
    #     plt.axvline(x=i)
    # # plt.ylim([0,5])
    # plt.xlim([2010, 2060])
    # plt.grid()
    # plt.savefig("pyplots/TOF_arrival_10k.png")


# theIs = [10,20,30,40,50,60,70,80,90]
# for j in theIs:
#     fitnessFile9Dir = os.path.abspath(os.path.join(simulation_output_dir, "GA_calculator/fitness_GA_EJ_%s.dat" %j))
#     fitnessFile9Array = np.genfromtxt(fitnessFile9Dir, delimiter=",")[:,1:3]
#     fitnessFile9ArrayTOFs = fitnessFile9Array[:,1]
#     fitnessFile9ArrayDVs = fitnessFile9Array[:,0]
#     minTOF9 = min(fitnessFile9ArrayTOFs)
#     avgTOF9 = np.mean(fitnessFile9ArrayTOFs)
#     avgDV9 = np.mean(fitnessFile9ArrayDVs)
#     # print("Fitness 0 min: %s" %(minTOF0/year))
#     # print("Fitness 0 mean: %s" %(avgTOF0/year))
#     # print("Fitness 9 min: %s" %(minTOF9/year))
#     # print("Fitness 9 mean: %s" %(avgTOF9/year))
#     # print("Fitness 9 mean DV: %s" %(avgDV9))
#
#     popFile9Dir = os.path.abspath(os.path.join(simulation_output_dir, "GA_calculator/population_GA_EJ_%s.dat" %j))
#     popFile9Array = np.genfromtxt(popFile9Dir, delimiter=",")[:,1:3]
#     launchYears = popFile9Array[:,0]/year
#     dummyList = np.ones(np.size(launchYears))
#
#
#     # Get rid of bad DV cases, and combine into full array:
#     fitnessFile9DVsListFINAL = []
#     fitnessFile9TOFSListFINAL = []
#     launchYearsListFINAL = []
#     arrivalYearsListFINAL = []
#     for i in range(len(fitnessFile9ArrayDVs)):
#         if fitnessFile9ArrayDVs[i] < 12000:
#             fitnessFile9DVsListFINAL.append(fitnessFile9ArrayDVs[i]/1000)
#             TOF = fitnessFile9ArrayTOFs[i]/year
#             fitnessFile9TOFSListFINAL.append(TOF)
#             launchYear = launchYears[i] + 2000
#             launchYearsListFINAL.append(launchYear)
#             arrivalYearsListFINAL.append(launchYear + TOF)
#
#
#
#     JupiterAps = np.arange(2005.267, 2060, 11.9)
#
#
#     plt.figure()
#     plt.scatter(launchYearsListFINAL, fitnessFile9TOFSListFINAL, 1)
#     plt.xlabel("Launch year")
#     plt.ylabel("Time of flight (years)")
#     # plt.ylim([0,5])
#     plt.xlim([2010, 2060])
#     plt.grid()
#     plt.savefig("pyplots/TOF_launch_12k_%s.png" %j)
#
#     plt.figure()
#     plt.scatter(fitnessFile9DVsListFINAL, fitnessFile9TOFSListFINAL, 1)
#     plt.xlabel("DVs (km/s)")
#     plt.ylabel("Time of flight (years)")
#     # plt.ylim([0,5])
#     plt.grid()
#     plt.savefig("pyplots/TOF_DV_12k_%s.png" %j)
#
#     plt.figure()
#     plt.scatter(arrivalYearsListFINAL, fitnessFile9TOFSListFINAL, 1)
#     plt.xlabel("Arrival year")
#     plt.ylabel("Time of flight (years)")
#     for i in JupiterAps:
#         plt.axvline(x=i, color='red')
#     plt.legend(["Jupiter Aphelion Times"])
#     # plt.ylim([0,5])
#     plt.xlim([2010, 2060])
#     plt.grid()
#     plt.savefig("pyplots/TOF_arrival_12k_%s.png" %j)

# # grid search look
# DVDir = os.path.abspath(os.path.join("/home/matt/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/SimulationOutput/THIS_IS_A_TEST.dat"))
# DVArray = np.genfromtxt(DVDir, dtype=float, delimiter='\t')[:,1:]
# launchDir =os.path.abspath(os.path.join("/home/matt/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/SimulationOutput/THIS_IS_A_TEST_x_data.dat"))
# launchArray = np.genfromtxt(launchDir, dtype=float, delimiter='\t')[:,1:]
# TOFDir =os.path.abspath(os.path.join("/home/matt/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/SimulationOutput/THIS_IS_A_TEST_y_data.dat"))
# TOFArray = np.genfromtxt(TOFDir, dtype=float, delimiter='\t')[:,1:]
# print(np.min(DVArray))
#
#
# plt.figure()
# for i in range(len(DVArray[:,0])):
#     plt.scatter(DVArray[:,i]/1000, TOFArray/year, 1)
# plt.xlabel("DVs (km/s)")
# plt.ylabel("Time of flight (years)")
# plt.xlim([2.5,10])
# plt.ylim([8.5,9.68])
# plt.grid()
# plt.title("Gride")
# # plt.savefig("pyplots/TOF_DV.png")


plt.show()


















# plt.figure(figsize=(10,10))
# plt.title("x against z")
# plt.xlabel("x coord [AU]")
# plt.ylabel("z coord [AU]")
# plt.xlim([-1, 1])
# plt.ylim([-1, 1])
# plt.grid()
# plt.plot(testDataArray[:, 1]/AU, testDataArray[:, 3]/AU)


# plt.figure(figsize=(10,10))
# plt.xlabel("y coord [AU]")
# plt.ylabel("z coord [AU]")
# plt.xlim([-1, 1])
# plt.ylim([-1, 1])
# plt.grid()
# plt.plot(testDataArray[:, 2]/AU, testDataArray[:, 3]/AU)
