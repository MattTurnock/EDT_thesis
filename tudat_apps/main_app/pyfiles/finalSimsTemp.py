import numpy as np
import os
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils

from matplotlib import pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 30})

year = 365*24*60*60
AU = 1.496E11

doingSSO = True
doingSOKGA = False
doingInO = False
doingCompare = False



DEFAULTSIZE = None
SaveFolder = "pyplots/finalSimsTempv2"

##################### SSO ########################

if doingSSO:

    # Plot the base case
    dataSubDir_SSO = "SSO"
    allSimData_SSO = utils.getAllSimDataFromFolder(dataSubDir_SSO)
    propData_SSO = allSimData_SSO[5]
    bodyData_SSO = allSimData_SSO[0]

    thrustData_SSO = allSimData_SSO[6]
    thrustLocal = thrustData_SSO[:,5:8]
    timesThrust = thrustData_SSO[:,0]/year
    altitudes = bodyData_SSO[:, 1]/AU


    # plt.figure()
    # plt.plot(timesThrust, thrustLocal[:,0])
    # plt.plot(timesThrust, thrustLocal[:,1])
    #
    # plt.figure()
    # plt.scatter(altitudes, thrustLocal[:,0])
    # plt.scatter(altitudes, thrustLocal[:,1])

    utils.plotTrajectoryData(propData_SSO, sameScale=True, planetsToPlot=["Earth"], plotSun=True, fignumber=1, plotOnlyTrajectory=False, trajectoryLabel="SSO Spacecraft", saveFolder=SaveFolder, savename="SSO-trajectory.pdf")

    utils.plotTrajectoryData(propData_SSO, fignumber=2, trajectoryLabel="BOOB", dataArrayType="propData", plotType="time-altitude", logScaleY=True, saveFolder=SaveFolder, savename="SSO-time-altitude.pdf")

    utils.plotTrajectoryData(propData_SSO, fignumber=3, trajectoryLabel="BOOB2", dataArrayType="propData", plotType="time-speed", logScaleY=True, saveFolder=SaveFolder, savename="SSO-time-speed.pdf")



######################### SOKGA ###########################################

if doingSOKGA:
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





    dataSubDir_SOKGA_reference = "SOKGA-Reference"
    allSimData_SOKGA_reference = utils.getAllSimDataFromFolder(dataSubDir_SOKGA_reference)
    propData_SOKGA_reference = allSimData_SOKGA_reference[5]
    bodyData_SOKGA_reference = allSimData_SOKGA_reference[0]

    thrustData_SOKGA_reference = allSimData_SOKGA_reference[6]
    # thrustLocal = thrustData_SSO[:,5:8]
    # timesThrust = thrustData_SSO[:,0]/year
    # altitudes = bodyData_SSO[:, 1]/AU





    V_reference = np.linalg.norm(propData_SOKGA_reference[-1, 4:])
    V_SOKGA = np.linalg.norm(propData_SOKGA_combined[-1, 4:])

    print("Reference V: ", V_reference)
    print("SOKGA V: ", V_SOKGA)


    customLegendTrajectory = ["SOKGA Nominal", "GA Reference", "Earth", "Jupiter", "Sun"]
    utils.plotTrajectoryData(propData_SOKGA_combined, figsize=[14,14], sameScale=True, planetsToPlot=["Earth", "Jupiter"], plotSun=True, fignumber=10, plotOnlyTrajectory=True)#, trajectoryLabel="SOKGA Spacecraft", saveFolder=SaveFolder, savename="SOKGA-test2.png")
    utils.plotTrajectoryData(propData_SOKGA_reference,  sameScale=True, planetsToPlot=["Earth", "Jupiter"], plotSun=True, fignumber=10, plotOnlyTrajectory=False, trajectoryLabel="SOKGA Spacecraft Reference", saveFolder=SaveFolder, savename="SOKGA-trajectory.pdf", legendLabelsCustom=customLegendTrajectory)


    customLegendTimeAltitude = ["SOKGA Nominal", "GA Reference"]
    utils.plotTrajectoryData(propData_SOKGA_combined, fignumber=11, trajectoryLabel="BOOB", dataArrayType="propData", plotType="time-altitude", logScaleY=True, plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_reference, fignumber=11, trajectoryLabel="PEEPEE", dataArrayType="propData", plotType="time-altitude", logScaleY=True, plotOnlyTrajectory=False, legendLabelsCustom=customLegendTimeAltitude, saveFolder=SaveFolder, savename="SOKGA-TimeAltitude.pdf")

    customLegendTimeSpeed = ["SOKGA Nominal", "GA Reference"]
    utils.plotTrajectoryData(propData_SOKGA_stage2, fignumber=12, trajectoryLabel="BOOB2", dataArrayType="propData", plotType="time-speed", logScaleY=True, plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_reference, fignumber=12, trajectoryLabel="PEEPEE2", dataArrayType="propData", plotType="time-speed", logScaleY=True, plotOnlyTrajectory=False, legendLabelsCustom=customLegendTimeSpeed, saveFolder=SaveFolder, savename="SOKGA-TimeSpeed.pdf")



##################################### InO #################################################

if doingInO:

    # Plot the InO
    dataSubDir_InO = "InO"
    allSimData_InO = utils.getAllSimDataFromFolder(dataSubDir_InO)
    propData_InO = allSimData_InO[5]
    bodyData_InO = allSimData_InO[0]

    thrustData_InO = allSimData_InO[6]
    thrustLocal = thrustData_InO[:,5:8]
    timesThrust = thrustData_InO[:,0]/year
    altitudes = bodyData_InO[:, 1]/AU



    print("Time of termination: ", timesThrust[-1] + 2000)



    dataSubDir_InO_stage2 = "InO-Stage2"
    allSimData_InO_stage2 = utils.getAllSimDataFromFolder(dataSubDir_InO_stage2)
    propData_InO_stage2 = allSimData_InO_stage2[5]
    bodyData_InO_stage2 = allSimData_InO_stage2[0]

    thrustData_InO_stage2 = allSimData_InO_stage2[6]
    thrustLocal = thrustData_InO_stage2[:,5:8]
    timesThrust = thrustData_InO_stage2[:,0]/year
    altitudes = bodyData_InO_stage2[:, 1]/AU


    propData_InO_Combined = np.concatenate((propData_InO, propData_InO_stage2))

    InOLegendLabelsCustom = ["InO Spacecraft Stage 1", "InO Spacecraft Stage 2", "Earth", "Mercury", "Sun"]
    utils.plotTrajectoryData(propData_InO, sameScale=True, planetsToPlot=["Earth", "Mercury"], plotSun=True, fignumber=30, plotOnlyTrajectory=True)#, trajectoryLabel="InO Spacecraft", saveFolder=SaveFolder, savename="InO-test.pdf")
    utils.plotTrajectoryData(propData_InO_stage2, sameScale=True, planetsToPlot=["Earth", "Mercury"], plotSun=True, fignumber=30, plotOnlyTrajectory=False, trajectoryLabel="InO Spacecraft", legendLabelsCustom=InOLegendLabelsCustom, saveFolder=SaveFolder, savename="InO-trajectory.pdf")

    utils.rescaleAndSavePlot(fignumber=30, xlims=[-10,1.2], ylims=[-4,1.5], saveFolder=SaveFolder, savename="InO-trajectory-zoomed.pdf")


    utils.plotTrajectoryData(propData_InO_Combined, fignumber=31, dataArrayType="propData", plotType="time-altitude", logScaleY=True, saveFolder=SaveFolder, savename="InO-time-altitude.pdf")

    utils.plotTrajectoryData(propData_InO_Combined, fignumber=32, dataArrayType="propData", plotType="time-speed", logScaleY=True, saveFolder=SaveFolder, savename="InO-time-speed.pdf")

if doingCompare:
    customLegendTrajectoriesCombined = ["SSO", "Nominal SOKGA", "GA Reference", "InO", "Earth", "Jupiter", "Sun"]
    utils.plotTrajectoryData(propData_SSO, fignumber=40, plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_combined, fignumber=40, plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_reference, fignumber=40, plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_InO_Combined, fignumber=40, plotOnlyTrajectory=False, sameScale=True, planetsToPlot=["Earth", "Jupiter"], plotSun=True, legendLabelsCustom=customLegendTrajectoriesCombined, saveFolder=SaveFolder, savename="Compared-trajectories.pdf")

    utils.rescaleAndSavePlot(fignumber=40, xlims=[-10,10], ylims=[-5.5,5.5], saveFolder=SaveFolder, savename="Compared-trajectories-zoomed.pdf")

    customLegendTimeAltitudeCombined = ["SSO", "Nominal SOKGA", "GA Reference", "InO"]
    utils.plotTrajectoryData(propData_SSO, fignumber=41, plotType="time-altitude", plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_combined, fignumber=41, plotType="time-altitude", plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_reference, fignumber=41, plotType="time-altitude", plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_InO_Combined, fignumber=41, plotType="time-altitude", plotOnlyTrajectory=False, logScaleY=True, legendLabelsCustom=customLegendTimeAltitudeCombined, saveFolder=SaveFolder, savename="Compared-time-altitudes.pdf")


    customLegendTimeSpeedCombined = ["SSO", "Nominal SOKGA", "GA Reference", "InO"]
    utils.plotTrajectoryData(propData_SSO, fignumber=42, plotType="time-speed", plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_combined, fignumber=42, plotType="time-speed", plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_SOKGA_reference, fignumber=42, plotType="time-speed", plotOnlyTrajectory=True)
    utils.plotTrajectoryData(propData_InO_Combined, fignumber=42, plotType="time-speed", plotOnlyTrajectory=False, logScaleY=True, legendLabelsCustom=customLegendTimeSpeedCombined, saveFolder=SaveFolder, savename="Compared-time-speeds.pdf")



plt.show()