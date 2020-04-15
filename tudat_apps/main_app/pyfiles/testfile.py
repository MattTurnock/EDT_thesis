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

def testPlot(filename, parName, parVal, scale):

    testDataFileDir = os.path.abspath(os.path.join(simulation_output_dir, filename))

    testDataArray = np.genfromtxt(testDataFileDir, delimiter=",")
    # print(testDataArray[:,0])

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


# testPlot("SSO-CHB-Test/SSO-CH-Test-out-0-0.dat", "0ms2",1.5)
# testPlot("SSO-CHB-Test/SSO-CH-Test-out-0-000325.dat", "0-000325ms2",75)
# testPlot("SSO-CHB-Test/SSO-CH-Test-out-0-000325.dat", "0-000325ms2",8)
# testPlot("SSO-CHB-Test/SSO-CH-Test-in-0-000325.dat", "-0-000325ms2",1.1)
# testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-0-000325.dat", "0-000325ms2",75)
# testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-expo-0-000325.dat", "0-000325ms2",75)
testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-expo-0-000325.dat", "constAccel", "0-000325ms2",75)
testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-expo-enviro__E-7.dat", "constCurrent", "10^-7 T",4)
testPlot("SSO-CHB-Test-custom/SSO-CH-Test-out-expo-enviro__E-6.dat", "constCurrent", "10^-6 T",10)

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
