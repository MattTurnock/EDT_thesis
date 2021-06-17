import os
import numpy as np
import json
from tudatApplications.EDT_thesis.tudat_apps.main_app.pyfiles import utils
import copy
import sys
import glob
import shutil
from matplotlib import pyplot as plt
import matplotlib

# # Set logfile
logfileDir = utils.pythonRunnerLogfilesDir
utils.checkFolderExist(logfileDir)
logfilePath = os.path.join(logfileDir, utils.createLogfileName("ProfileSensitivityRunner"))
utils.logger = utils.createLogger(logfilePath)

simulationRunPath = os.path.join(utils.cppApplications_dir, "application_simulation_SSO-CHB")
utils.checkFolderExist(utils.outputJsonDirPath_ProfileSensitivity, emptyDirectory=False)

loadOnly = True
showing = False
individualLoadingBars = False
doPlotting = True

matplotlib.rcParams.update({'font.size': 30})
markersize = 10
profileSensitivityPlotFolder = os.path.join(utils.pyplots_dir, "profileSensitivity")
utils.checkFolderExist(profileSensitivityPlotFolder)



##################### Actual Running ########################################


barOver = utils.IncrementalBar("Profile Sensitivity Overall Progress: " , max=len(utils.changeKeyLines_ProfileSensitivity), suffix='[%(percent).1f%%. ETA: %(eta)ds]   ')
# For loop to go over each profile
for i in range(len(utils.changeKeyLines_ProfileSensitivity)):
    changeKeyLine = utils.changeKeyLines_ProfileSensitivity[i]
    sensitivityRange = utils.sensitivityRanges_ProfileSensitivity[i]
    parameterName = utils.parameterNames_ProfileSensitivity[i]

    sensitivitySpace = np.linspace(sensitivityRange[0], sensitivityRange[1], utils.profileSensitivitySpacing)

    SSOP_Fitness, SSOP_AllSimData = utils.runProfileSensitivity(sensitivitySpace, profileName="SSOP", parameterName=parameterName,
                                changeKeyLine=changeKeyLine, errorTolerance=utils.profileSensitivityTolerance,
                                numpySavePath=utils.numpySavePath_ProfileSensitivity,
                                numpySaveNameBase=utils.numpySaveName_ProfileSensitivity,
                                loadOnly=loadOnly, loadingBars=individualLoadingBars)

    InO_Fitness, InO_AllSimData = utils.runProfileSensitivity(sensitivitySpace, profileName="InO", parameterName=parameterName,
                                changeKeyLine=changeKeyLine, errorTolerance=utils.profileSensitivityTolerance,
                                numpySavePath=utils.numpySavePath_ProfileSensitivity,
                                numpySaveNameBase=utils.numpySaveName_ProfileSensitivity,
                                loadOnly=loadOnly, loadingBars=individualLoadingBars)

    SOKGA_TOFs, SOKGA_Vels, SOKGA_AllSimData = utils.runProfileSensitivity_SOKGA(sensitivitySpace, parameterName=parameterName,
                                               changeKeyLine=changeKeyLine, errorTolerance=utils.profileSensitivityTolerance,
                                                                                 loadOnly=loadOnly)


##################### Plot Generation ########################################
    if doPlotting:

        ylabel = "Maximum Ap [AU]"

        # DO the SSOP case plots
        plt.figure(figsize=utils.figSizeDefault)
        ax2 = plt.gca()
        ax2.get_yaxis().get_major_formatter().set_useOffset(False)
        ax2.ticklabel_format(useOffset=False, style="plain")
        plt.xlabel(utils.parameterPlotNames_ProfileSensitivity[i])
        plt.ylabel(ylabel)
        plt.grid(which="both")

        closestNominalIndex = utils.findNearestInArray(sensitivitySpace, utils.nominalValues_ProfileSensitivity[i])[1]
        plt.plot(sensitivitySpace, SSOP_Fitness , c="C0")
        # plt.plot(utils.nominalValues_ProfileSensitivity[i], utils.SSOP_NominalFitness, 'o', markersize=markersize, c="C2")
        plt.plot(sensitivitySpace[closestNominalIndex], SSOP_Fitness[closestNominalIndex], 'o', markersize=markersize, c="C2")


        plt.legend(["Parameter Trend", "Nominal Case"])
        plt.savefig(os.path.join(profileSensitivityPlotFolder, "par-ap-%s-%s.png" %("SSOP", parameterName)), bbox_inches='tight')
        plt.savefig(os.path.join(profileSensitivityPlotFolder, "par-ap-%s-%s.pdf" %("SSOP", parameterName)), bbox_inches='tight')


        # DO the InO case plots
        plt.figure(figsize=utils.figSizeDefault)
        ax2 = plt.gca()
        ax2.get_yaxis().get_major_formatter().set_useOffset(False)
        ax2.ticklabel_format(useOffset=False, style="plain")
        plt.xlabel(utils.parameterPlotNames_ProfileSensitivity[i])
        plt.ylabel(ylabel)
        plt.grid(which="both")

        closestNominalIndex = utils.findNearestInArray(sensitivitySpace, utils.nominalValues_ProfileSensitivity[i])[1]
        plt.plot(sensitivitySpace, InO_Fitness , c="C0")
        # plt.plot(utils.nominalValues_ProfileSensitivity[i], utils.InO_NominalFitness, 'o', markersize=markersize, c="C2")
        plt.plot(sensitivitySpace[closestNominalIndex], InO_Fitness[closestNominalIndex], 'o', markersize=markersize, c="C2")

        plt.legend(["Parameter Trend", "Nominal Case"])
        plt.savefig(os.path.join(profileSensitivityPlotFolder, "par-ap-%s-%s.png" %("InO", parameterName)), bbox_inches='tight')
        plt.savefig(os.path.join(profileSensitivityPlotFolder, "par-ap-%s-%s.pdf" %("InO", parameterName)), bbox_inches='tight')


        # DO the SOKGA case plots
        plt.figure(figsize=utils.figSizeDefault)
        ax2 = plt.gca()
        ax2.get_yaxis().get_major_formatter().set_useOffset(False)
        ax2.ticklabel_format(useOffset=False, style="plain")
        plt.xlabel(utils.parameterPlotNames_ProfileSensitivity[i])
        plt.ylabel("TOF [years]")
        plt.grid(which="both")

        closestNominalIndex = utils.findNearestInArray(sensitivitySpace, utils.nominalValues_ProfileSensitivity[i])[1]
        plt.plot(sensitivitySpace, SOKGA_TOFs / utils.year , c="C0")
        # plt.plot(utils.nominalValues_ProfileSensitivity[i], utils.InO_NominalFitness, 'o', markersize=markersize, c="C2")
        plt.plot(sensitivitySpace[closestNominalIndex], SOKGA_TOFs[closestNominalIndex]/utils.year, 'o', markersize=markersize, c="C2")

        plt.legend(["Parameter Trend", "Nominal Case"])
        plt.savefig(os.path.join(profileSensitivityPlotFolder, "par-TOF-%s-%s.png" %("SOKGA", parameterName)), bbox_inches='tight')
        plt.savefig(os.path.join(profileSensitivityPlotFolder, "par-TOF-%s-%s.pdf" %("SOKGA", parameterName)), bbox_inches='tight')




        plt.figure(figsize=utils.figSizeDefault)
        ax2 = plt.gca()
        ax2.get_yaxis().get_major_formatter().set_useOffset(False)
        ax2.ticklabel_format(useOffset=False, style="plain")
        plt.xlabel(utils.parameterPlotNames_ProfileSensitivity[i])
        plt.ylabel("Velocity [km/s]")
        plt.grid(which="both")

        closestNominalIndex = utils.findNearestInArray(sensitivitySpace, utils.nominalValues_ProfileSensitivity[i])[1]
        plt.plot(sensitivitySpace, SOKGA_Vels / 1000 , c="C0")
        # plt.plot(utils.nominalValues_ProfileSensitivity[i], utils.InO_NominalFitness, 'o', markersize=markersize, c="C2")
        plt.plot(sensitivitySpace[closestNominalIndex], SOKGA_Vels[closestNominalIndex]/1000, 'o', markersize=markersize, c="C2")

        plt.legend(["Parameter Trend", "Nominal Case"])
        plt.savefig(os.path.join(profileSensitivityPlotFolder, "par-Vel-%s-%s.png" %("SOKGA", parameterName)), bbox_inches='tight')
        plt.savefig(os.path.join(profileSensitivityPlotFolder, "par-Vel-%s-%s.pdf" %("SOKGA", parameterName)), bbox_inches='tight')



    barOver.next()
barOver.finish()

if showing:
    plt.show()
