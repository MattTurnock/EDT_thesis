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

# Plot the base case
dataSubDir_SSO = "SSO"
dataSubDir_SSO = "SSO-Configs-Sensitivity/SSO-currents-20"
allSimData_SSO = utils.getAllSimDataFromFolder(dataSubDir_SSO)
propData_SSO = allSimData_SSO[5]
speed = np.linalg.norm(propData_SSO[:, 4:7], axis=1)
print("speed: ", speed)
# bodyData_SSO = allSimData_SSO[0]

propDataTimes = propData_SSO[:,0]

thrustData_SSO = allSimData_SSO[6]
magData_SSO = allSimData_SSO[4]
currentVNVData_SSO = allSimData_SSO[7]
# thrustLocal = thrustData_SSO[:,5:8]
# timesThrust = thrustData_SSO[:,0]/year
# altitudes = bodyData_SSO[:, 1]/AU

utils.plotTrajectoryData(propData_SSO, sameScale=True, planetsToPlot=["Earth"], plotSun=True, fignumber=1, plotOnlyTrajectory=False, trajectoryLabel="SSO Spacecraft", saveFolder=SaveFolder, savename="SSO-test.png")


# newtimes = np.arange(propDataTimes[0], propData_SSO[1], 100000)
# newtimes = np.arange(propDataTimes[0], propDataTimes[-1], 1000000)
# print(len(propDataTimes))
# print(len(newtimes))

newtimes = np.linspace(propDataTimes[0], propDataTimes[-1], 2*len(propDataTimes))
propData_SSO_Interpolated = utils.interpolateArrays(propData_SSO, newtimes)

thrustData_SSO_Interpolated = utils.interpolateArrays(thrustData_SSO, newtimes)
thrustVectorX = thrustData_SSO[:,2]
thrustVectorY = thrustData_SSO[:,3]
thrustMagnitude = thrustData_SSO[:,1]

magData_SSO_Interpolated = utils.interpolateArrays(magData_SSO, newtimes)
magDataX = magData_SSO[:,6]
print(magDataX)
magDataY = magData_SSO[:,7]
thetaData = np.rad2deg(magData_SSO[:,2])

utils.plotTrajectoryData(propData_SSO_Interpolated, sameScale=True, planetsToPlot=["Earth"], plotSun=True, fignumber=2, plotOnlyTrajectory=False, trajectoryLabel="SSO Spacecraft")#, saveFolder=SaveFolder, savename="SSO-test.png")












# Plot the SOKGA case
dataSubDir_SOKGA = "SOKGA"
allSimData_SOKGA = utils.getAllSimDataFromFolder(dataSubDir_SOKGA)
propData_SOKGA = allSimData_SOKGA[5]
bodyData_SOKGA = allSimData_SOKGA[0]

thrustData_SOKGA = allSimData_SOKGA[6]
# thrustLocal = thrustData_SSO[:,5:8]
# timesThrust = thrustData_SSO[:,0]/year
# altitudes = bodyData_SSO[:, 1]/AU

# utils.plotTrajectoryData(propData_SOKGA, sameScale=True, planetsToPlot=["Earth", "Jupiter"], plotSun=True, fignumber=2, plotOnlyTrajectory=False, trajectoryLabel="SOKGA Spacecraft", saveFolder=SaveFolder, savename="SOKGA-test.png")



dataSubDir_SOKGA_stage2 = "SOKGA-Stage2"
allSimData_SOKGA_stage2 = utils.getAllSimDataFromFolder(dataSubDir_SOKGA_stage2)
propData_SOKGA_stage2 = allSimData_SOKGA_stage2[5]
bodyData_SOKGA_stage2 = allSimData_SOKGA_stage2[0]

thrustData_SOKGA_stage2 = allSimData_SOKGA_stage2[6]
# thrustLocal = thrustData_SSO[:,5:8]
# timesThrust = thrustData_SSO[:,0]/year
# altitudes = bodyData_SSO[:, 1]/AU

propData_SOKGA_combined = np.concatenate((propData_SOKGA, propData_SOKGA_stage2))

propDataTimes_SOKGA = propData_SOKGA_combined[:,0]
newtimes_SOKGA = np.linspace(propDataTimes_SOKGA[0], propDataTimes_SOKGA[-1], 2*len(propDataTimes_SOKGA))
propData_SOKGA_Interpolated = utils.interpolateArrays(propData_SOKGA_combined, newtimes_SOKGA)


utils.plotTrajectoryData(propData_SOKGA_Interpolated, figsize=[14,14], sameScale=True, planetsToPlot=["Earth", "Jupiter"], plotSun=True, fignumber=3, plotOnlyTrajectory=True)#, trajectoryLabel="SOKGA Spacecraft", saveFolder=SaveFolder, savename="SOKGA-test2.png")



dataSubDir_SOKGA_reference = "SOKGA-Reference"
allSimData_SOKGA_reference = utils.getAllSimDataFromFolder(dataSubDir_SOKGA_reference)
propData_SOKGA_reference = allSimData_SOKGA_reference[5]
bodyData_SOKGA_reference = allSimData_SOKGA_reference[0]

thrustData_SOKGA_reference = allSimData_SOKGA_reference[6]
# thrustLocal = thrustData_SSO[:,5:8]
# timesThrust = thrustData_SSO[:,0]/year
# altitudes = bodyData_SSO[:, 1]/AU

customLegend = ["SOKGA Nominal", "GA Reference", "Earth", "Jupiter", "Sun"]

propDataTimes_SOKGA_reference = propData_SOKGA_reference[:,0]
newtimes_SOKGA_reference = np.linspace(propDataTimes_SOKGA_reference[0], propDataTimes_SOKGA_reference[-1], 2*len(propDataTimes_SOKGA_reference))
propData_SOKGA_reference_Interpolated = utils.interpolateArrays(propData_SOKGA_reference, newtimes_SOKGA_reference)

utils.plotTrajectoryData(propData_SOKGA_reference_Interpolated,  sameScale=True, planetsToPlot=["Earth", "Jupiter"], plotSun=True, fignumber=4, plotOnlyTrajectory=False, trajectoryLabel="SOKGA Spacecraft Reference", saveFolder=SaveFolder, savename="SOKGA-test.pdf", legendLabelsCustom=customLegend)










fignumber = 5
figsize = [20,20]
plotSun = True
saveFolder = None
savename = None
xlims = None
ylims = None
sameScale = True
planetsToPlot = ["Earth", "Jupiter"]
plotOnlyTrajectory = False
runtime = 10





plt.figure(fignumber, figsize=figsize)
fig = plt.gcf()
ax = fig.gca()
# ax2 = plt.subplots(1,1,1)


# add another axes at the top left corner of the figure
axtext = fig.add_axes([0.0,0.9,0.1,0.05])
# turn the axis labels/spines/ticks off
axtext.axis("off")
# place the text to the other axes
time = axtext.text(0.5,0.5, str(0), ha="left", va="top")



# create the parametric curve
# t=np.arange(0, 2*np.pi, 2*np.pi/100)
nummer = int(len(propData_SOKGA_Interpolated)*0)
nummer2 = int(len(propData_SOKGA_Interpolated)*0.5)
x=propData_SOKGA_Interpolated[nummer:nummer2,1] / AU
y=propData_SOKGA_Interpolated[nummer:nummer2,2] / AU
timeList = (propData_SOKGA_Interpolated[nummer:, 0] / year) + 2000
timeList = (propData_SSO_Interpolated[:, 0] / year) + 2000
# z=t/(2.*np.pi)

x = propData_SSO_Interpolated[:,1] / AU
y = propData_SSO_Interpolated[:,2] / AU

x = propData_SSO[:,1] / AU
y = propData_SSO[:,2] / AU

timeList = (propData_SSO[:,0]/year) + 2000
# create the first plot
point, = ax.plot([x[0]], [y[0]], 'o')
line, = ax.plot(x, y, label='SOKGA Nominal')

legendLabels = []
ax.axis("scaled")
ax.set_xlabel("X coordinate [AU]")
ax.set_ylabel("Y coordinate [AU]")
ax.grid()
ax.set_xlim([-10,10])
ax.set_ylim([-10,10])

############### ADD AN ARROW INTO THE MIX ##################
xCoords = x
yCoords = y
Qx = xCoords[:, None]
Qy = yCoords[:, None]


UCoords = thrustVectorX * 3
VCoords = thrustVectorY * 3

# scaleFactor = (1/magDataX[0])/50
# UCoords = magDataX * scaleFactor
# VCoords = magDataY * scaleFactor
U = UCoords[:, None]
V = VCoords[:, None]

# For frame 1, plot the 0th set of arrows
s = np.s_[0, :]
qax = ax.quiver(Qx[s], Qy[s], U[s], V[s],
                facecolor="red", scale=10)

print(UCoords)
print(VCoords)


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

if plotSun:
    ax.scatter([0], [0], c="orange")
    legendLabels.append("Sun")


ax.legend()
# ax.set_xlim([-1.5, 1.5])
# ax.set_ylim([-1.5, 1.5])
# ax.set_zlim([-1.5, 1.5])
# time.set_text("reeeeee")
# second option - move the point position at every frame

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
print(speed)

decimalRounding = 3
def update_point(n, x, y):
    point.set_data(np.array([x[n], y[n]]))
    # point.set_3d_properties(z[n], 'z')
    currentYear = timeList[n]
    currentDate = (utils.decimalYearToDatetimeObject(currentYear)).date()
    theta = thetaData[n]


    time.set_text(timeStringBase %(np.around(currentYear, decimals=decimalRounding),
                                       np.around(theta, decimals=decimalRounding),
                                       np.around(thrustMagnitude[n], decimals=decimalRounding),
                                       np.around(speed[n]/1000, decimals=decimalRounding),
                                       np.around(Ems[n], decimals=decimalRounding),
                                       np.around(I0s[n], decimals=decimalRounding),
                                       np.around(ics[n], decimals=decimalRounding),
                                       np.around(Iavg[n], decimals=decimalRounding)))
    # return point, time


    ### ARROW ANIMATION ###
    # Update to frame i
    s = np.s_[n, :]
    # Change direction of arrows
    qax.set_UVC(U[s], V[s])
    # Change base position of arrows
    qax.set_offsets(np.c_[Qx[s].flatten(), Qy[s].flatten()])

interval = (runtime / len(x))*1000



ani=animation.FuncAnimation(fig, update_point, frames=len(x), fargs=(x, y), interval=interval, repeat=False, blit=False)
# ani.save("test.mp4")

plt.show()


