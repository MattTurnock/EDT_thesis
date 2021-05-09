import datetime

import numpy as np
import os
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils

matplotlib.rcParams.update({'font.size': 20})

year = 365*24*60*60
AU = 1.496E11


##################### GA Plots ########################
DEFAULTSIZE = None
SaveFolder = "pyplots/finalSimsTemp"


dataSubDir_SSO = "SSO"
allSimData_SSO_raw = utils.getAllSimDataFromFolder(dataSubDir_SSO)
allSimData_SSO = utils.interpolateAllDataArrays(allSimData_SSO_raw, dataRange=[0,0.1], forcedArrayLength=100)
propData_SSO = allSimData_SSO[5]
speed = np.linalg.norm(propData_SSO[:, 4:7], axis=1)


propDataTimes = propData_SSO[:,0]

thrustData_SSO = allSimData_SSO[6]
magData_SSO = allSimData_SSO[4]
currentVNVData_SSO = allSimData_SSO[7]

thrustVectorX = thrustData_SSO[:,2]
thrustVectorY = thrustData_SSO[:,3]
thrustMagnitude = thrustData_SSO[:,1]

magDataX = magData_SSO[:,6]
magDataY = magData_SSO[:,7]




















# Set parameters that are to be passed into the eventual function
fignumber = 5
figsize = [20,20]
plotSun = True
saveFolder = None
sameScale = True
xlims = [-10, 10]
ylims = [-10, 10]
planetsToPlot = ["Earth", "Jupiter"]
plotOnlyTrajectory = False
runtime = 10
dataTextBoxDimensions = [0.0, 0.9, 0.1, 0.05]  # Dimensions of the text box containing live data (eg date etc). Format is [left, bottom, width, height]



# Create figure and axes
plt.figure(fignumber, figsize=figsize)
fig = plt.gcf()
ax = fig.gca()


# add another axes at the top left corner of the figure for text data
axtext = fig.add_axes(dataTextBoxDimensions)
# turn the axis labels/spines/ticks off
axtext.axis("off")
# place the text to the other axes
dataText = axtext.text(0.5, 0.5, str(0), ha="left", va="top")



# Create x and y data for the animation point plot
x = propData_SSO[:,1] / AU
y = propData_SSO[:,2] / AU


# create the first plot - ie the first point and the whole trajectory path
point, = ax.plot([x[0]], [y[0]], 'o')
line, = ax.plot(x, y, label='SOKGA Nominal')

# Add some generic plotting stuff
legendLabels = []
if sameScale: ax.axis("scaled")
ax.set_xlabel("X coordinate [AU]")
ax.set_ylabel("Y coordinate [AU]")
ax.grid()
ax.set_xlim(xlims)
ax.set_ylim(ylims)


#### Add arrow plotting info ####
# Arrow coordinates
Qx = x[:, None]
Qy = y[:, None]

# Arrow pointing vector
UCoords = thrustVectorX * 3
VCoords = thrustVectorY * 3

U = UCoords[:, None]
V = VCoords[:, None]

# For frame 1, plot the 0th set of arrows
s = np.s_[0, :]
qax = ax.quiver(Qx[s], Qy[s], U[s], V[s],
                facecolor="red", scale=10)



# Add planets to plot where wanted
if len(planetsToPlot) != 0:
    for i in range(len(planetsToPlot)):
        currentPlanet = planetsToPlot[i]
        if currentPlanet == "Jupiter":
            circleRadius = 5.2
            circleColour = 'r'
            legendName = "Jupiter"
        elif currentPlanet == "Earth":
            circleRadius = 1.0
            circleColour = 'g'
            legendName = "Earth"
        elif currentPlanet == "Mercury":
            circleRadius = 0.4
            circleColour = 'brown'
            legendName = "Mercury"

        circleToPlot = plt.Circle((0, 0), circleRadius, color=circleColour, fill=False)
        ax.add_patch(circleToPlot)
        legendLabels.append(legendName)

# Add Sun spot
if plotSun:
    ax.scatter([0], [0], c="orange")
    legendLabels.append("Sun")

# Add legend to plot
ax.legend(legendLabels)

dataTextToAdd = ["Year", "Theta", "Thrust", "Speed", "EM", "I0", "ic", "Iavg"]
dataTextToAddArrays = [ (propData_SSO[:,0]/year) + 2000,
                        np.rad2deg(magData_SSO[:,2]),
                        currentVNVData_SSO[:, 1],
                        currentVNVData_SSO[:, 2],
                        currentVNVData_SSO[:, 3],
                        currentVNVData_SSO[:, 6]]

# Create text data to show over plot
timeList = (propData_SSO[:,0]/year) + 2000
thetaData = np.rad2deg(magData_SSO[:,2])
Ems = currentVNVData_SSO[:, 1]
I0s = currentVNVData_SSO[:,2]
ics = currentVNVData_SSO[:,3]
Iavg = currentVNVData_SSO[:, 6]

timeStringBase = "Date: %s\n" \
                 "Theta: %s\n" \
                 "Thrust [N]: %s\n" \
                 "Speed: [km/s]: %s\n" \
                 "EMs: %s\n" \
                 "I0s: %s\n" \
                 "ics: %s\n" \
                 "Iavg: %s"

decimalRounding = 3


# Define plotting function
def update_point(n, x, y):

    # Perform Spacecraft Point animation
    point.set_data(np.array([x[n], y[n]]))


    # Perform arrow animation
    # Update to frame i
    s = np.s_[n, :]
    # Change direction of arrows
    qax.set_UVC(U[s], V[s])
    # Change base position of arrows
    qax.set_offsets(np.c_[Qx[s].flatten(), Qy[s].flatten()])

    # Add data text to plot
    dataText.set_text(timeStringBase % (np.around(timeList[n], decimals=decimalRounding),
                                        np.around(thetaData[n], decimals=decimalRounding),
                                        np.around(thrustMagnitude[n], decimals=decimalRounding),
                                        np.around(speed[n]/1000, decimals=decimalRounding),
                                        np.around(Ems[n], decimals=decimalRounding),
                                        np.around(I0s[n], decimals=decimalRounding),
                                        np.around(ics[n], decimals=decimalRounding),
                                        np.around(Iavg[n], decimals=decimalRounding)))

    return point, qax, dataText

interval = (runtime / len(x))*1000
print(len(x))



ani=animation.FuncAnimation(fig, update_point, frames=len(x), fargs=(x, y), interval=interval, repeat=False, blit=False)
print("saving")
writer = animation.FFMpegWriter(
    fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("test.mp4")

plt.show()


