//
// Created by matt on 19/04/2020.
//

#ifndef TUDATBUNDLE_GENERAL_FUNCTIONS_H
#define TUDATBUNDLE_GENERAL_FUNCTIONS_H

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h>

using namespace tudat;
using namespace tudat::mathematical_constants;

namespace gen {

    //////////////////// Create a bunch of useful input-output strings /////////////////////////

    // Get .../tudatBundle/tudat/Tudat rootpath and create .../tudatBundle/ root path, for use in other paths
    std::string tudatTudatRootPath = input_output::getTudatRootPath();
    boost::filesystem::path tudatTudatRootPathBoost(tudatTudatRootPath);
    boost::filesystem::path tudatBundleRootPathBoost = tudatTudatRootPathBoost.parent_path().parent_path().parent_path();
    std::string tudatBundleRootPathTemp = tudatBundleRootPathBoost.string();
    std::string tudatBundleRootPath = tudatBundleRootPathTemp + "/";

    // tudat rootpath
    std::string tudatRootPath = tudatBundleRootPath + "tudat/";
    // EDT thesis project rootpath
    std::string EDT_thesisRootPath = tudatBundleRootPath + "tudatApplications/EDT_thesis/";
    // tudat_apps inside EDT project path
    std::string tudat_appsRootPath = EDT_thesisRootPath + "tudat_apps/";
    // main_app inside tudat_apps path
    std::string main_appRootPath = tudat_appsRootPath + "main_app/";
    // JsonInputs inside main_app path
    std::string jsonInputsRootPath = main_appRootPath + "JsonInputs/";
    //SimulationOutput inside main_app path
    std::string simulationOutputRootPath = main_appRootPath + "SimulationOutput/";
    // Custom kernels MUST be placed into tudat kernels directory to work properly. This is that path
    std::string tudatKernelsRootPath = tudatRootPath + "Tudat/External/SpiceInterface/Kernels/";

    // Function to read json file (from JsonInputs directory)
    nlohmann::json readJson(std::string filename, std::string fileDirectory=gen::jsonInputsRootPath){

        std::string pathToJson = fileDirectory + filename;
        std::ifstream jsonFile(pathToJson);
        nlohmann::json jsonObject = nlohmann::json::parse(jsonFile);

        return jsonObject;
    }

    double getCircleArea(double diameter){
        double radius = 0.5*diameter;
        double area =  PI * std::pow(radius, 2);

        return area;
    }

    double getDonutArea(double diameterInner, double diameterOuter){
        double areaInner = gen::getCircleArea(diameterInner);
        double areaOuter = gen::getCircleArea(diameterOuter);
        double donutArea = areaOuter - areaInner;

        return donutArea;
    }

    // Functions to convert between LVLH and Inertial reference frames given the rotation angle between them (ie theta_)
    Eigen::Vector3d LvlhToInertial(Eigen::Vector3d LvlhVector, double rotationAngle){
        double XLvlh = LvlhVector[0];
        double YLvlh = LvlhVector[1];
        double ZLvlh = LvlhVector[2];

        double Xinertial = XLvlh * sin(rotationAngle) + YLvlh * cos(rotationAngle);
        double Yinertial = XLvlh * cos(rotationAngle) - YLvlh * sin(rotationAngle);
        double Zinertial = 0;

        Eigen::Vector3d inertialVector;
        inertialVector << Xinertial, Yinertial, Zinertial;
        return inertialVector;
    }

    Eigen::Vector3d InertialToLvlh(Eigen::Vector3d inertialVector, double rotationAngle){
        // Note: conversion is the same both ways, but 2 functions written to be explicit
        Eigen::Vector3d lvlhVector = LvlhToInertial(inertialVector, rotationAngle);

        return lvlhVector;
    }
}

#endif //TUDATBUNDLE_GENERAL_FUNCTIONS_H
