//
// Created by matt on 19/04/2020.
//

#ifndef TUDATBUNDLE_GENERAL_FUNCTIONS_H
#define TUDATBUNDLE_GENERAL_FUNCTIONS_H

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h>
#include "applicationOutput.h"

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

    /////////////////////////// Some useful constants and conversion factors //////////////////////

    double pi = 3.14159;
    double deg2rad = pi / 180;
    double AU = 1.496e11;


    /////////////////////////// Some useful functions //////////////////////

    // Function to read json file (from JsonInputs directory)
    nlohmann::json readJson(std::string filename, std::string fileDirectory=gen::jsonInputsRootPath){

        std::string pathToJson = fileDirectory + filename;
        std::ifstream jsonFile(pathToJson);
        nlohmann::json jsonObject = nlohmann::json::parse(jsonFile);

        return jsonObject;
    }

    // Calculate circle area
    double getCircleArea(double diameter){
        double radius = 0.5*diameter;
        double area =  PI * std::pow(radius, 2);

        return area;
    }

    // Calculate donut area (of 2 circles)
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

    // Functions to convert vetween maglocal (ie Br Bphi) and lvlh coordinates
    Eigen::Vector3d MaglocalToLvlh(Eigen::Vector3d maglocalVector){
        Eigen::Vector3d lvlhVector;
        lvlhVector << maglocalVector[1], maglocalVector[0], maglocalVector[2];
        return lvlhVector;
    }

    Eigen::Vector3d LvlhToMaglocal(Eigen::Vector3d lvlhVector){
        Eigen::Vector3d maglocalVector;
        maglocalVector << lvlhVector[1], lvlhVector[0], lvlhVector[2];
        return maglocalVector;
    }

    // Convert tudat time (seconds since j2000) to the (decimal) year
    double tudatTime2DecimalYear(double tudatTime){
        // convert by starting at year 2000, and adding seconds as years
        double decimalYear = 2000 + tudatTime/60/60/24/365.25;
        return decimalYear;
    }

    // Function to return 2 sines superimposed
    double twoSines(double x, std::vector<double> sinParameters ){
        double a1 = sinParameters[0];
        double b1 = sinParameters[1];
        double c1 = sinParameters[2];
        double a2 = sinParameters[3];
        double b2 = sinParameters[4];
        double c2 = sinParameters[5];
        double d = sinParameters[6];

        double y = a1*sin(b1*x + c1) + a2*sin(b2*x + c2) + d;
        return y;
    }

    // Function to convert year to seconds
    double years2Seconds(double years){
        return years * 365.25*24*60*60;
    }

    // Function to convert year (eg 2020) to MJD
    double year2MJDSeconds(double year){
        return years2Seconds(year - 2000);
    }

    // Function to convert km/s to m/s
    double km2m(double kms){
        return kms*1000;
    }

    // Function to print pagmo population to file (based on one in saveOptimizationResults.h)
    void printPopulationToFile( const std::vector< std::vector< double > >& population,
                                const std::string filePrefix,
                                const std::string fileSuffix,
                                const std::string outputSubFolder,
                                const bool isFitness )
    {

        Eigen::MatrixXd matrixToPrint( population.size( ), population.at( 0 ).size( ) );
        for( unsigned int i = 0; i < population.size( ); i++ )
        {
            for( unsigned int j = 0; j < population.at( 0 ).size( ); j++ )
            {
                matrixToPrint( i, j ) = population.at( i ).at( j );
            }
        }

        if( !isFitness )
        {
            tudat::input_output::writeMatrixToFile( matrixToPrint, filePrefix + "population_" + fileSuffix + ".dat", 16,
                                                    tudat_applications::getOutputPath( ) + outputSubFolder );
        }
        else
        {
            tudat::input_output::writeMatrixToFile( matrixToPrint, filePrefix + "fitness_" + fileSuffix + ".dat", 16,
                                                    tudat_applications::getOutputPath( ) + outputSubFolder  );
        }
    }


    ////////////////////////// Some Equations and formulae for bare tethers ///////////////////////

    // Function to convert from true current to dimensionless current through EDT, or vice versa
    double trueCurrentDimensionlessConvert(double currentToConvert, double unitCurrent){
        return currentToConvert / unitCurrent;
    }

    // Function to get unit current I_0
    double getUnitCurrent(double conductivity, double motionalEMF, double tetherCrossSectionalArea){
        return conductivity * motionalEMF * tetherCrossSectionalArea;
    }

    // Function to get motional EMF Em from simulation parameters
    double getMotionalEMF(Eigen::Vector3d velocityWrtMagField, Eigen::Vector3d magFieldVector, Eigen::Vector3d EDTLengthUnitVector){
        return (velocityWrtMagField.cross(magFieldVector)).dot(EDTLengthUnitVector);
    }

    // Function to get average unit current i_avg
    double getAvgDimensionlessCurrent(double tetherLength, double dimensionlessVoltageA){
        return -(1/5*tetherLength)*pow(dimensionlessVoltageA, 5/2) + 0.5*pow(dimensionlessVoltageA, 3/2);
    }

    // Function to get dimensionless voltage at A lambda_A, using a known value for ic
    double getDimensionlessVoltageA(double dimensionlessCurrentC){
        double ic = dimensionlessCurrentC;
        return pow( (2*ic - pow(ic, 2)), 2/3 );
    }








}

#endif //TUDATBUNDLE_GENERAL_FUNCTIONS_H
