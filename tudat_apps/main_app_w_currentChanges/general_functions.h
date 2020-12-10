//
// Created by matt on 19/04/2020.
//

#ifndef TUDATBUNDLE_GENERAL_FUNCTIONS_H
#define TUDATBUNDLE_GENERAL_FUNCTIONS_H

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h>
#include "applicationOutput.h"

using namespace tudat;
using namespace tudat::mathematical_constants;
using namespace tudat::simulation_setup;
using namespace tudat::numerical_integrators;

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
//    double absoluteTerminationTimeYears = 2099;


    ////////////////////////// Functions for time-cnversions ///////////////////////////////////////

    // Convert tudat time (seconds since j2000) to the (decimal) year
    double tudatTime2DecimalYear(double tudatTime){
        // convert by starting at year 2000, and adding seconds as years
        double decimalYear = 2000 + tudatTime/60/60/24/365.25;
        return decimalYear;
    }

    // Function to convert year to seconds
    double years2Seconds(double years){
        return years * 365.25*24*60*60;
    }

    // FUnction to convert year to days
    double years2Days(double years){
        return years * 365.25;
    }

    // Function to convert year (eg 2020) to MJD, in seconds
    double year2MJDSeconds(double year){
        return years2Seconds(year - 2000);
    }

    // Function to convert year (eg 2020) to MJD, in days
    double year2MJDDays(double year){
        return years2Days(year - 2000);
    }

    // Function to convert year (eg 2020) to year, in MJD format (eg 20)
    double year2MJDYears(double year){
        return year-2000;
    }


    ////////////////////////// Functions for frame conversions ///////////////////////////////////////

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


    ////////////////////////// Functions for GA related things ///////////////////////////////////////

    // Function to normalise or denormlaise deltaV or TOF, based on given bounds
    double normaliseValue(double initialValue, double lowerBound, double upperBound, bool denormalise = false) {
        double newValue;
        // if statement true, then denormalised
        if (denormalise) {
            double denormalisedValue = initialValue * (upperBound - lowerBound) + lowerBound;
            newValue = denormalisedValue;
        }
            // otherwise normalised
        else {
            double normalisedValue = (initialValue - lowerBound) / (upperBound - lowerBound);
            newValue = normalisedValue;
        }

        return newValue;
    }

    // Function to convert from numerical <int> vector for flyby trajectory to named bodies
    std::vector< std::string > FlybySequence2TransferBodyTrajectory(std::vector< int > flybySequence){

        // Create list of transfer bodies, and list to put the transfer trajectory into. Note, sun is a dummy
        std::vector< std::string > transferBodies = {"Sun",
                                                     "Mercury",
                                                     "Venus",
                                                     "Earth",
                                                     "Mars",
                                                     "Jupiter",
                                                     "Saturn",
                                                     "Uranus",
                                                     "Neptune"};
        std::vector< std::string > transferBodyTrajectory;

        // For each entry in the sequence, add to the new transferBodyTrajectoryVector
        int noFlybys = flybySequence.size();
        for(int i = 0; i < noFlybys; i++){
            transferBodyTrajectory.push_back( transferBodies[ flybySequence[i] ] );
        }

        return transferBodyTrajectory;
    }


    ////////////////////////// Functions for Optimisation things ///////////////////////////////////////

    // Function to create (and run and save) a grid search for a problem. Based on applicationOutput.h in tudat examples
    void createGridSearch(
            pagmo::problem& problem,
            const std::vector< std::vector< double > >& bounds,
            const std::vector< int > numberOfPoints,
            const std::string& fileName,
            std::string outputFolder)
    {
        if( bounds.at( 0 ).size( ) != 2 )
        {
            std::cerr<<"Warning when plotting grid search, size of problem does not equal 2"<<std::endl;
        }
        Eigen::MatrixXd gridSearch = Eigen::MatrixXd( numberOfPoints.at( 0 ), numberOfPoints.at( 1 ) );

        double xSpacing = ( bounds[ 1 ][ 0 ] - bounds[ 0 ][ 0 ] ) / static_cast< double >( numberOfPoints.at( 0 ) - 1 );
        double ySpacing = ( bounds[ 1 ][ 1 ] - bounds[ 0 ][ 1 ] ) / static_cast< double >( numberOfPoints.at( 1 ) - 1 );

        std::vector< double > decisionVector;
        decisionVector.resize( 2 );

        std::vector< double > xDataPoints;
        for( int i = 0; i < numberOfPoints.at( 0 ); i++ )
        {
            xDataPoints.push_back( bounds[ 0 ][ 0 ] + static_cast< double >( i ) * xSpacing );
        }

        std::vector< double > yDataPoints;
        for( int j = 0; j < numberOfPoints.at( 1 ); j++ )
        {
            yDataPoints.push_back( bounds[ 0 ][ 1 ] + static_cast< double >( j ) * ySpacing );
        }

        for( int i = 0; i < numberOfPoints.at( 0 ); i++ )
        {
            std::cout<<"Grid search "<<i<<std::endl;
            for( int j = 0; j < numberOfPoints.at( 1 ); j++ )
            {
                decisionVector[ 0 ] = xDataPoints[ i ];
                decisionVector[ 1 ] = yDataPoints[ j ];

                gridSearch( i, j ) = problem.fitness( decisionVector ).at( 0 );
            }
        }
        std::cout << "outputDirectory" << tudat_applications::getOutputPath( ) << std::endl;
        tudat::input_output::writeMatrixToFile( gridSearch, fileName + ".dat" , 16, tudat_applications::getOutputPath() + outputFolder );
        tudat::input_output::writeMatrixToFile( tudat::utilities::convertStlVectorToEigenVector(
                xDataPoints ), fileName + "_x_data.dat", 16, tudat_applications::getOutputPath() + outputFolder );
        tudat::input_output::writeMatrixToFile( tudat::utilities::convertStlVectorToEigenVector(
                yDataPoints ), fileName + "_y_data.dat", 16, tudat_applications::getOutputPath() + outputFolder );

    }


    ////////////////////////// Functions for Simulation related things (eg integrators etc) ////////////////////////////

    // Return the integrator coefficient set from string input
    RungeKuttaCoefficients::CoefficientSets getIntegratorCoefficientSet(std::string integratorString){

        RungeKuttaCoefficients::CoefficientSets coefficientSetOutput;
        if (integratorString == "RK87DP"){ // TODO: add support for more integrators!
            coefficientSetOutput = numerical_integrators::RungeKuttaCoefficients::rungeKutta87DormandPrince;
        }
        else{
            std::cout << "INTEGRATOR NOT RECOGNISED: terminating simulation" << std::endl;
            exit(1);
        }

        return coefficientSetOutput;
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
                                                    tudat_applications::getOutputPath( ) + outputSubFolder,
                                                    ",");
        }
        else
        {
            tudat::input_output::writeMatrixToFile( matrixToPrint, filePrefix + "fitness_" + fileSuffix + ".dat", 16,
                                                    tudat_applications::getOutputPath( ) + outputSubFolder,
                                                    ",");
        }
    }

//    // Function used to facilitate end of simulation based on a year
//    bool absoluteTimeTermination(double time){
//        double absoluteTerminationTime = year2MJDSeconds(absoluteTerminationTimeYears);
//        bool terminateSim ;
//
//        if (time > absoluteTerminationTime){
//            terminateSim = true;
//        }
//        else{
//            terminateSim = false;
//        }
//
//        return terminateSim;
//    }

    // Function used to facilitate end of simulation based on a year
    bool absoluteTimeTermination(double time, double absoluteTerminationTimeYear){
        double absoluteTerminationTime = year2MJDSeconds(absoluteTerminationTimeYear);
        bool terminateSim ;

        if (time > absoluteTerminationTime){
            terminateSim = true;
        }
        else{
            terminateSim = false;
        }

        return terminateSim;
    }

//    std::function< bool(const double ) > absoluteTimeTerminationFunction = &absoluteTimeTermination;

    ////////////////////////// Cpp utility functions (eg json loading, vector type conversion...) /////////////////////

    // Function to read json file (from JsonInputs directory)
    nlohmann::json readJson(std::string filename, std::string fileDirectory=gen::jsonInputsRootPath){

        std::string pathToJson = fileDirectory + filename;
        std::ifstream jsonFile(pathToJson);
        nlohmann::json jsonObject = nlohmann::json::parse(jsonFile);

        return jsonObject;
    }

    // Function to convert double to string, with decimal precision
    std::string floatToString(double number, int decimalPlaces=2){
        std::stringstream stream;
        stream << std::fixed << std::setprecision(decimalPlaces) << number;
        std::string s = stream.str();

        return s;
    }

    // Function to multiply all values of a (double) vector by a constant (double)
    std::vector< double > vectorScaling(std::vector<double> inputVector, double scaleFactor){
        std::vector<double> outputVector;
        for( int i=0 ; i<inputVector.size() ; i++){
            outputVector.push_back( scaleFactor * inputVector[i] );
        }
        return outputVector;
    }

    // COnvert an eigen::vector to std::vector (for doubles)
    std::vector< double > EigenVectorXd2StdVectorDoubles(Eigen::VectorXd inputVector){
        std::vector< double > outputVector;
        outputVector.resize(inputVector.size());
        Eigen::VectorXd::Map(&outputVector[0], inputVector.size()) = inputVector;

        return outputVector;
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
        return -(1.0/5.0*tetherLength)*pow(dimensionlessVoltageA, 5.0/2.0) + 0.5*pow(dimensionlessVoltageA, 3.0/2.0);
    }

    // Function to get dimensionless voltage at A lambda_A, using a known value for ic
    double getDimensionlessVoltageA(double dimensionlessCurrentC){
        double ic = dimensionlessCurrentC;
        return pow( (2*ic - pow(ic, 2)), 2.0/3.0 );
    }

    /// FOllowing for general tethers ///

    // Calculate circle area TODO: Remove me
    double getCircleArea(double diameter){
        double radius = 0.5*diameter;
        double area =  PI * std::pow(radius, 2);

        return area;
    }

    // Calculate donut area (of 2 circles) TODO: Remove me
    double getDonutArea(double diameterInner, double diameterOuter){
        double areaInner = gen::getCircleArea(diameterInner);
        double areaOuter = gen::getCircleArea(diameterOuter);
        double donutArea = areaOuter - areaInner;

        return donutArea;
    }

    // Calculate Diameter of a circle given the area
    double calculateCircleDiameter(double area){

        double radius = std::pow( area/PI, 0.5 );
        return radius*2;
    }

    // Calculate Diameter of the tether core, given the 2 area



    ////////////////////////// Misc Functions  /////////////////////

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

    // Function to convert km/s to m/s
    double km2m(double kms){
        return kms*1000;
    }




}

#endif //TUDATBUNDLE_GENERAL_FUNCTIONS_H
