import os
import numpy as np
import json
from utils import cppApplications_dir, jsonInputs_dir, startYears, endYears, quickConfigs, runBashCommand # This can come up redlined sometimes, but still works

jsonName = "GAConfigsNominal.json"
jsonPath = os.path.join(jsonInputs_dir, jsonName)
GAApplicationName = "application_GA_calculator"

baseGACommand = os.path.join(cppApplications_dir, GAApplicationName)
GAArgumentsList = [baseGACommand, jsonName]
printing=0

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
    print("WARNING: Planet not recognised")

with open(jsonPath, 'r+') as f:
    # Load json data into a variable
    GAVariables = json.load(f)

    for i in range(len(startYears)):

        # Set and change start / end year and other variables in the GA variables
        startYear = startYears[i]
        endYear = endYears[i]
        # outputSubFolderBase = "GA_calculator_%s" %(planetName)

        GAVariables["PlanetConfigs"][planetConfigString]["Bounds"]["StartYearLower"] = startYear
        GAVariables["PlanetConfigs"][planetConfigString]["Bounds"]["StartYearUpper"] = endYear
        GAVariables["PlanetConfigs"]["planetToFlyby"] = planetName
        # GAVariables["saveDataConfigs"]["outputSubFolder"] = outputSubFolder


        # Dump new data into a json file TODO: make a bunch of new json files in a separate directory, to track changes?
        f.seek(0)
        json.dump(GAVariables, f, indent=4)
        f.truncate()

        # Run simulations using relevant json file
        print("Running %s simulations for %s - %s" %(planetName, startYear, endYear))
        runBashCommand( GAArgumentsList, printSetting=printing)


