import numpy as np
import os
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils

from matplotlib import pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 20})

year = 365*24*60*60
AU = 1.496E11

##################### GA Plots ########################
DEFAULTSIZE = None
VnVSaveFolder = "pyplots/VnV"

# Reference tolerance
dataSubDir_GA_reference = "integratorGAReferenceVnV"
allSimData_GA_reference = utils.getAllSimDataFromFolder(dataSubDir_GA_reference, todoList=["propData"])
propData_GA_reference = allSimData_GA_reference[5]


utils.plotTrajectoryData(propData_GA_reference, sameScale=True, planetsToPlot=["Jupiter", "Earth"], plotSun=True, fignumber=1, plotOnlyTrajectory=True, trajectoryLabel="Reference Spacecraft")

# Nominal tolerance
dataSubDir_GA = "integratorGAVnV"
allSimData_GA = utils.getAllSimDataFromFolder(dataSubDir_GA, todoList=["propData"])
propData_GA = allSimData_GA[5]

propData_GA_reference_resized = utils.reverseInterpolateArrays(propData_GA, propData_GA_reference)
times_GA_reference_resized = propData_GA_reference_resized[:,0]
pos_GA_reference_resized = propData_GA_reference_resized[:, 1:4]
vel_GA_reference_resized = propData_GA_reference_resized[:, 4:]

times_GA = propData_GA[:,0]
times_GA_years = times_GA/365.25/24/60/60 + 2000
pos_GA = propData_GA[:, 1:4]
vel_GA = propData_GA[:, 4:]

imposedLegendGA= ["Reference Spacecraft", "Nominal Spacecraft", "Jupiter", "Earth", "Sun"]
utils.plotTrajectoryData(propData_GA, sameScale=True, planetsToPlot=["Jupiter", "Earth"], plotSun=True, fignumber=1,
                         plotOnlyTrajectory=False, trajectoryLabel="Nominal Spacecraft", saveFolder=VnVSaveFolder,
                         savename="Integrator_GA.pdf", legendLabelsCustom=imposedLegendGA, figsize=DEFAULTSIZE, newFontSize=15)
# utils.rescaleAndSavePlot(1, [-5.3,5.3], [-5.3,5.3], VnVSaveFolder, "Integrator_GA_Zoom.pdf")


# Difference plotting

# pos_GA_interpolated = utils.interpolateArrays(pos_GA, times_GA_reference)
pos_GA_difference = pos_GA - pos_GA_reference_resized
pos_GA_norm_diff = np.linalg.norm(pos_GA_difference, axis=1)


maxValue = np.amax(pos_GA_norm_diff)
maxIndex = np.where(pos_GA_norm_diff == maxValue)[0][0]
print("==================== SOME VALUES FOR GA DATA ============")
print("Maximum error [AU]: ", maxValue/AU)
print("At index: ", maxIndex)
print("Reference Data Coordinates [AU]: ", pos_GA_reference_resized[maxIndex]/AU)
print("Nominal Data Coordinates [AU]: ", pos_GA[maxIndex]/AU)

plt.figure(10, figsize=utils.figSizeDefault)
plt.plot(times_GA_years[0:-1], pos_GA_norm_diff[0:-1]/AU)
plt.ylabel("Position Difference [AU]")
plt.xlabel("Time [year]")
plt.grid()
plt.savefig(os.path.join(VnVSaveFolder, "Integrator_GA_Differences.pdf"))


##################### Inner SS Plot ########################

# Reference Tolerance
dataSubDir_inner_reference = "integratorInnerReferenceVnV"
allSimData_inner_reference = utils.getAllSimDataFromFolder(dataSubDir_inner_reference, todoList=["propData"])
propData_inner_reference = allSimData_inner_reference[5]


utils.plotTrajectoryData(propData_inner_reference, sameScale=True, planetsToPlot=["Earth", "Jupiter", "Mercury"], plotSun=True, fignumber=3, plotOnlyTrajectory=True, trajectoryLabel="Reference Spacecraft")


# Nominal Tolerance
dataSubDir_inner = "integratorInnerVnV"
allSimData_inner = utils.getAllSimDataFromFolder(dataSubDir_inner, todoList=["propData"])
propData_inner = allSimData_inner[5]
# propData_inner = utils.interpolateArrays(propData_inner_raw, times_inner_reference)

propData_inner_reference_resized = utils.reverseInterpolateArrays(propData_inner, propData_inner_reference)
times_inner_reference_resized = propData_inner_reference_resized[:,0]
pos_inner_reference_resized = propData_inner_reference_resized[:, 1:4]
vel_inner_reference_resized = propData_inner_reference_resized[:, 4:]

times_inner = propData_inner[:,0]
times_inner_years = times_inner/365.25/24/60/60 + 2000
pos_inner = propData_inner[:, 1:4]
vel_inner = propData_inner[:, 4:]


imposedLegendInner = ["Reference Spacecraft", "Nominal Spacecraft", "Earth", "Jupiter", "Mercury", "Sun"]
utils.plotTrajectoryData(propData_inner, sameScale=True, planetsToPlot=["Earth", "Jupiter", "Mercury"], plotSun=True,
                         fignumber=3, plotOnlyTrajectory=False, trajectoryLabel="Nominal SPacecraft",
                         saveFolder=VnVSaveFolder, savename="Integrator_inner.pdf", legendLabelsCustom=imposedLegendInner,
                         figsize=DEFAULTSIZE, ylims=[-10,10], legendLocation='upper left')
# plt.show()
utils.rescaleAndSavePlot(3, [-0.5,0.5], [-0.5,0.5], VnVSaveFolder, "Integrator_inner_Zoom.pdf")


# Difference plotting
pos_inner_difference = pos_inner - pos_inner_reference_resized
pos_inner_norm_diff = np.linalg.norm(pos_inner_difference, axis=1)

maxValue_Inner = np.amax(pos_inner_norm_diff)
maxIndex_Inner = np.where(pos_inner_norm_diff == maxValue_Inner)[0][0]
print("==================== SOME VALUES FOR INNER DATA ============")
print("Maximum error [AU]: ", maxValue_Inner/AU)
print("At index: ", maxIndex_Inner)
print("Reference Data Coordinates [AU]: ", pos_inner_reference_resized[maxIndex_Inner]/AU)
print("Nominal Data Coordinates [AU]: ", pos_inner[maxIndex_Inner]/AU)

plt.figure(30, figsize=utils.figSizeDefault)
plt.plot(times_inner_years, pos_inner_norm_diff/AU)
plt.ylabel("Position Difference [AU]")
plt.xlabel("Time [year]")
plt.grid()
plt.savefig(os.path.join(VnVSaveFolder, "Integrator_inner_Differences.pdf"))


plt.show()