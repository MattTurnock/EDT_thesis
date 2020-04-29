import numpy as np
import os
from matplotlib import pyplot as plt

AU = 1.496e11 # AU in m

# Set various project directories
pyfiles_dir = os.path.dirname(os.path.realpath(__file__))
main_app_dir = os.path.abspath(os.path.join(pyfiles_dir, os.pardir))
simulation_output_dir = os.path.abspath(os.path.join(main_app_dir, "SimulationOutput"))

################### TESTING OUTPUT FROM SSO-CHB TEST 1 #####################

"SSO-CHB-Test/SSO-CH-Test-in-0-000325.dat"

def testPlot(filename, parName, parVal, scale=10, plotType = "propData", filename2=""):

    year = 365*24*60*60
    AU = 1.496E11

    testDataFileDir = os.path.abspath(os.path.join(simulation_output_dir, filename))
    testDataArray = np.genfromtxt(testDataFileDir, delimiter=",")

    if filename2 != "":
        testDataFileDir2 = os.path.abspath(os.path.join(simulation_output_dir, filename2))
        testDataArray2 = np.genfromtxt(testDataFileDir, delimiter=",")

    times = testDataArray[:,0]/year

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

        # plt.savefig("pyplots/%s-%s-%s.png" %(parName, parVal, scale))

    elif plotType == "magField":
        legendList = []

        # Create figure, title, scale etc
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        plt.title("Magnetic field strength plot")
        plt.xlabel("time (years)")
        plt.ylabel("Magfield strength (T)")
        plt.grid()

        # Add local and inertial magfield plots
        plt.plot(testDataArray[:, 0]/year, testDataArray[:, 1])
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

        # plt.savefig("pyplots/%s-%s-%s.png" %(parName, parVal, scale))

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

        # plt.savefig("pyplots/%s-%s-%s.png" %(parName, parVal, scale))

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

        # plt.savefig("pyplots/%s-%s-%s.png" %(parName, parVal, scale))

    elif plotType == "altitudePlot":
        legendList = []

        # Create figure, title, scale etc
        plt.figure(figsize=(10,10))
        ax = plt.gca()
        plt.title("Altitude plot")
        plt.xlabel("time (years)")
        plt.ylabel("Distance (AU)")
        plt.grid()

        # Create altitude list
        positionArray = testDataArray[:, 1:4]
        altitudes = np.zeros(len(testDataArray[:,0]))
        print(altitudes)
        for i in range(len(testDataArray[:,0])):
            altitudes[i] = np.linalg.norm(positionArray[i])
        altitudes = altitudes/AU

        # Plot altitudes
        plt.plot(times, altitudes)
        legendList.append("alt")

        # Add legend
        plt.legend(legendList)

        # plt.savefig("pyplots/%s-%s-%s.png" %(parName, parVal, scale))


# testPlot("SSO-CHB-Test/SSO-CH-Test-out-0-0.dat", "0ms2",1.5)
# testPlot("SSO-CHB-Test/SSO-CH-Test-out-0-000325.dat", "0-000325ms2",75)
# testPlot("SSO-CHB-Test/SSO-CH-Test-out-0-000325.dat", "0-000325ms2",8)
# testPlot("SSO-CHB-Test/SSO-CH-Test-in-0-000325.dat", "-0-000325ms2",1.1)
# testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-0-000325.dat", "0-000325ms2",75)
# testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-expo-0-000325.dat", "0-000325ms2",75)

# testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-expo-0-000325.dat", "constAccel", "0-000325ms2",scale=75)
# testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-expo-enviro__E-7.dat", "constCurrent", "10^-7 T",scale=4)
# testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-expo-enviro__E-6.dat", "constCurrent", "10^-6 T",scale=10)

testPlot("SSO-CHB-Test-custom-2/SSO-CH-Test-out-expo-enviro__E-6_propData.dat", "const current", "subheading", plotType="propData")
testPlot("SSO-CHB-Test-custom-2/SSO-CH-Test-out-expo-enviro__E-6_magData.dat", "Theta over time", "subheading", plotType="posTheta")
testPlot("SSO-CHB-Test-custom-2/SSO-CH-Test-out-expo-enviro__E-6_magData.dat", "Magfield over time", "subheading", plotType="magField")
testPlot("SSO-CHB-Test-custom-2/SSO-CH-Test-out-expo-enviro__E-6_magData.dat", "Magfield over time", "subheading", plotType="magField2")
testPlot("SSO-CHB-Test-custom-2/SSO-CH-Test-out-expo-enviro__E-6_depVarData.dat", "Kepler plot", "subheading", plotType="keplerPlot1")
testPlot("SSO-CHB-Test-custom-2/SSO-CH-Test-out-expo-enviro__E-6_depVarData.dat", "Altitude plot", "subheading", plotType="altitudePlot")


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
