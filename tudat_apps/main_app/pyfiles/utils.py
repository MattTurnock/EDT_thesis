import os
import numpy as np
import subprocess

############################## Set various project directories ############################################
pyfiles_dir = os.path.dirname(os.path.realpath(__file__))
main_app_dir = os.path.abspath(os.path.join(pyfiles_dir, os.pardir))
simulation_output_dir = os.path.abspath(os.path.join(main_app_dir, "SimulationOutput"))
tudatBundle_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(main_app_dir))))
cppApplications_dir = os.path.abspath(os.path.join(tudatBundle_dir, "tudatApplications", "EDT_thesis", "tudat_apps", "bin", "applications"))
jsonInputs_dir = os.path.join(main_app_dir, "JsonInputs")
# print(os.path.join(cppApplications_dir, "application_GA_calculator"))
############################## Set constant values for use in GACalculatorRunner.py and others ######################
quickConfigs = ["Jupiter", 2020.580, 2050] #Quick setup configs to change in json on the fly. corresponds to: planetToFlyby, startYear, endYear TODO: Expand as more configs needed to be run
synodicPeriod = 398.8/365.25 #TODO: make synodic period automatically created
startYearsRaw = np.arange(quickConfigs[1], quickConfigs[2]+1, synodicPeriod)
endYearsRaw = startYearsRaw + synodicPeriod
startYears = np.round(startYearsRaw, 2)
endYears = np.round(endYearsRaw, 2)
############################# Some useful functions ####################################################

def runBashCommand(argumentsList, printSetting=2): #Printsettings: 0 = no printing; 1 = print in realtime; 2 = print all at end

    # # Define the full command, using arguments list and base command
    # fullCommand = baseCommand
    # if len(argumentsList) != 0:
    #     for arg in argumentsList:
    #         fullCommand = fullCommand + " " + arg
    #
    # print(fullCommand)
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
