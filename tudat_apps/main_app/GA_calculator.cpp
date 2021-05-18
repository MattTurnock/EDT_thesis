//
// Created by matt on 07/08/2020.
//

/////////////////// Include relevant header files ////////////////////

#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>

#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/simulated_annealing.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/nsga2.hpp"


#include "earthPlanetTransfer.h"


#include "simulation_SSO-CHB.h"
#include "GA_calculator.h"
//#include "general_functions.h"
//#include "universal_settings.h"

#include "EDTGuidance.h"
#include "universal_settings.h"
#include "environment_settings.h"

#include "../../../../tudatExampleApplications/libraryExamples/PaGMOEx/Problems/getAlgorithm.h"
#include "../../../../tudatExampleApplications/libraryExamples/PaGMOEx/Problems/saveOptimizationResults.h"

int main(int argc, char *argv[])
{

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat::basic_astrodynamics;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;



    // Read json file for general use (name of json is passed as first argument, or default set here)
    std::string jsonName;
    if (argc != 1){
        jsonName = argv[1];
    }
    else{
        jsonName = "GAConfigsNominal.json";
    }
    std::cout<< "Loading Configs from: " << jsonName << std::endl;
    // Set the subjson folders to file
    nlohmann::json dummyVariablesJson = gen::readJson("testVariables.json");
    nlohmann::json simulationVariables = gen::readJson(jsonName);
    nlohmann::json PlanetConfigs = simulationVariables["PlanetConfigs"];
    nlohmann::json AlgorithmConfigs = simulationVariables["AlgorithmConfigs"];
    std::string planetToFlyby = PlanetConfigs["planetToFlyby"];
    nlohmann::json GAPlanetSpecificConfigs = PlanetConfigs[planetToFlyby + "Configs"];
    std::cout << "Planet specific config: " << GAPlanetSpecificConfigs << std::endl;

    // Load some variables directly from json that are needed throughout
    double startYearLower = GAPlanetSpecificConfigs["Bounds"]["StartYearLower"];
    double startYearUpper = GAPlanetSpecificConfigs["Bounds"]["StartYearUpper"];
    std::string outputSubFolderBase = simulationVariables["saveDataConfigs"]["outputSubFolderBase"];
    std::string outputSubFolder = outputSubFolderBase + "_"
                                  + planetToFlyby + "_"
                                  + gen::floatToString(startYearLower) + "-" + gen::floatToString(startYearUpper);
    // Load custom spice kernels
    std::string customKernelName = simulationVariables["Spice"]["customKernelName"];
    std::string customSpiceKernelPath = gen::tudatKernelsRootPath + customKernelName;
    const std::vector<std::string> spicePathVector = {customSpiceKernelPath};
    spice_interface::loadStandardSpiceKernels( spicePathVector );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            Do Lambert solver stage              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////// Set up the problem, with bounds //////////////////////////////////////

    // Set seed for reproducible results
    unsigned int seed=123;
    pagmo::random_device::set_seed(seed);

    ///////////////  Set up all variables that are to be set in the upcoming sections  ///////////////////////

    // Uses 3 decision variables (start dates, flight duration, delta V) each with lower and upper bounds (which are
    // implemented with a vector of vectors)
    std::vector< std::vector< double > > Bounds(2, std::vector< double > (2, 0.0));
    std::vector< double > DVBounds(2);// TODO: finish adding variables

    // Create the flyby sequence (ie Earth Saturn / Earth Jupiter)
    std::vector< int > FlybySequence;

    // Create file suffix for saving
    std::string fileSuffix;

    if (planetToFlyby == "Jupiter") {
        // Define flyby sequence
        FlybySequence.push_back(3);
        FlybySequence.push_back(5);

        // Define file suffix
        fileSuffix = "GA_EJ_";
    }
    else if (planetToFlyby == "Saturn"){
        // Define flyby sequence
        FlybySequence.push_back(3);
        FlybySequence.push_back(6);

        // Define file suffix
        fileSuffix = "GA_ES_";
    }
    else if (planetToFlyby == "Mars"){
        // Define flyby sequence
        FlybySequence.push_back(3);
        FlybySequence.push_back(4);

        // Define file suffix
        fileSuffix = "GA_EM_";
    }
    else{
        std::cout << "PLANET NOT SUPPORTED: terminating simulation" << std::endl;
        return 1;
    }

    // Define bounds, by using values defined in the json
    Bounds[0][0] = gen::year2MJDYears(startYearLower);
    Bounds[1][0] = gen::year2MJDYears(startYearUpper);
    Bounds[0][1] = GAPlanetSpecificConfigs["Bounds"]["FlightDurationYearsLower"];
    Bounds[1][1] = GAPlanetSpecificConfigs["Bounds"]["FlightDurationYearsUpper"];
    DVBounds[0] = GAPlanetSpecificConfigs["Bounds"]["DeltaVKmsLower"];
    DVBounds[1] = GAPlanetSpecificConfigs["Bounds"]["DeltaVKmsUpper"];


    // Create problem object for problem fitness
    bool normaliseValues = AlgorithmConfigs["normaliseValues"];
    bool includeDepartureDV = AlgorithmConfigs["includeDepartureDV"];
    bool includeArrivalDV = AlgorithmConfigs["includeArrivalDV"];
    EarthPlanetTransfer earthPlanetProbStruct = EarthPlanetTransfer(Bounds, DVBounds, FlybySequence,
                                                                    normaliseValues, includeDepartureDV, includeArrivalDV, outputSubFolder);
    pagmo::problem Prob{ earthPlanetProbStruct };

    // Do a cheeky grid search if specified in json
    bool doGridSearch = AlgorithmConfigs["doGridSearch"];
    int gridSearchSize = AlgorithmConfigs["gridSearchSize"];
    if (doGridSearch) {
        std::cout << "Doing grid search for porkchop" << std::endl;
        std::string gridSearchFilename = "porkchop_" + fileSuffix;
        gridSearchFilename = gridSearchFilename.substr(0, gridSearchFilename.length() -1 ); // Remove _ from end of suffix

        gen::createGridSearch(Prob, Bounds,
                              {gridSearchSize, gridSearchSize}, gridSearchFilename, outputSubFolder);
    }

    // Select algorithm for problem, based on the one in the json TODO: expand possibilites of algorithsm that can be uesed with a custom function
    // NOTE: indices are: [0,1,2] = [nsga2, moead, ihs]. Can add more with custom get algorithm function if desired
    std::string algorithmName = AlgorithmConfigs["algorithmName"];
    int algorithmIndex;
    if (algorithmName == "nsga2"){
        algorithmIndex = 0;
    }
    else if (algorithmName == "moead"){
        algorithmIndex = 1;
    }
    else if (algorithmName == "ihs"){
        algorithmIndex = 2;
    }
    else{
        std::cout << "ALGORITHM NOT RECOGNISED: terminating simulation" << std::endl;
        return 1;
    }
    algorithm algo = getMultiObjectiveAlgorithm(algorithmIndex);
    // TODO: Set algorithm configs in json, or simply use default values (^^) again
//    pagmo::algorithm algo{pagmo::moead(1u,
//                                   "grid",
//                                   "tchebycheff",
//                                   20u,
//                                   1.0,
//                                   0.5,
//                                   20,
//                                   0.9,
//                                   2u,
//                                   true,
//                                   seed)};
    std::string info = algo.get_extra_info();
    std::cout << "extra info: " << info << std::endl;
//    moead::moead(unsigned gen, std::string weight_generation, std::string decomposition, population::size_type neighbours,
//            double CR, double F, double eta_m, double realb, unsigned limit, bool preserve_diversity, unsigned seed)

    // Create islands with individuals defined by json TODO: make json numbers reasonable. Also create a vector to dump all the islands into for reference in second stage of optimisation
    island Isl{algo, Prob, AlgorithmConfigs["islandSize"]};
    std::vector< island > islandsList;

//    // Initialise populations and fitness lists to add to for later
//    std::vector< std::vector< std::vector< double > > > populationsList;
//    std::vector< std::vector< std::vector< double > > > fitnessList;

    // Evolve for number of generations defined by json TODO: make json numbers reasonable
    int noGenerations = AlgorithmConfigs["noGenerations"];
    for( int i=0 ; i<noGenerations ; i++){

        // Evolve
        Isl.evolve();
        while( Isl.status( ) != pagmo::evolve_status::idle &&
                Isl.status( ) != pagmo::evolve_status::idle_error ){
            Isl.wait( );
        }
        Isl.wait_check( ); // Raises errors


        std::string filePrefix = ""; //TODO: decide if an actual prefix needed
        // Write current iteration results to file for Jupiter NOTE: outputs as "filePrefix populatin_/fitness_ fileSuffix .dat
        gen::printPopulationToFile(Isl.get_population( ).get_x( ),filePrefix, fileSuffix + std::to_string( i ), outputSubFolder, false );
        gen::printPopulationToFile(Isl.get_population( ).get_f( ),filePrefix, fileSuffix + std::to_string( i ), outputSubFolder,true );

//        populationsList.push_back(Isl.get_population().get_x());
//        fitnessList.push_back(Isl.get_population().get_f());

        // Add current iteration to vector of all iterations TODO: consider fiddling with this to reduce memory usage
        islandsList.push_back(Isl);

        // Print current generation
        std::cout<<i<<std::endl;

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            Second stage of optimisation              //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////// Basically just grab the island populations and use ////////////////////////////////////
    ///////////////////////// them to recreate trajectories, in order to save the initial conditions ////////////////

    for (int i=0; i<islandsList.size(); i++) {

        // Grab the islands for this iteration and also set the populations / fitnesses
        island thisIsland = islandsList[i];
        std::vector< std::vector< double > > thisPopulation = thisIsland.get_population().get_x();
        std::vector< std::vector< double > > thisFitness = thisIsland.get_population().get_f();

        // Quick calculation of number of legs
        int numberOfLegs = matt_trajectories::getNumberOfLegs(FlybySequence);

        // Map to save state into
        std::map < int, Eigen::VectorXd > thisStateMap;

        // For loop for each entry of the population
        for (int j=0; j<thisPopulation.size(); j++) {

            // Set the current entry and initialise variable vector
            std::vector< double > thisEntry = thisPopulation[j];
            Eigen::VectorXd variableVector;
            variableVector.resize(numberOfLegs + 1);

            // Add each part of the entry (ie departure date and leg TOFs) to the variable vector
            for (int k=0; k<thisEntry.size(); k++){
                variableVector[k] = thisEntry[k];
            }

            variableVector[numberOfLegs] = 1; //dummy
            variableVector *= physical_constants::JULIAN_YEAR; // TOF given in years

            Trajectory thisTrajectory(numberOfLegs,
                                      earthPlanetProbStruct.getLegTypeVector(),
                                      earthPlanetProbStruct.getEphemerisVector(),
                                      earthPlanetProbStruct.getGravitationalParameterVector(),
                                      variableVector,
                                      matt_trajectories::sunGravitationalParameter,
                                      earthPlanetProbStruct.getMinimumPericenterRadii(),
                                      earthPlanetProbStruct.getSemiMajorAxes(),
                                      earthPlanetProbStruct.getEccentricities(),
                                      includeDepartureDV,
                                      includeArrivalDV);

            double thisDV;
            thisTrajectory.calculateTrajectory(thisDV);
            double launchDate = variableVector[0] / physical_constants::JULIAN_YEAR;

//            std::cout << "Launch date (nominal)  : " << launchDate << std::endl;
////            std::cout << "Launch date (from traj): " << variableVector[0] << std::endl;
//            std::cout << "DV (nominal)           : " << thisFitness[j][0] << std::endl;
//            std::cout << "DV (from traj)         : " << thisDV << std::endl;


            Eigen::Vector3d	departureBodyPosition;
            Eigen::Vector3d	departureBodyVelocity;
            Eigen::Vector3d velocityAfterDeparture;
            thisTrajectory.getLaunchConditions(departureBodyPosition, departureBodyVelocity, velocityAfterDeparture);

            Eigen::VectorXd departureState;
            departureState.resize(8);
            departureState[0] = launchDate;
            departureState[1] = departureBodyPosition[0];
            departureState[2] = departureBodyPosition[1];
            departureState[3] = departureBodyPosition[2];
            departureState[4] = velocityAfterDeparture[0];// + departureBodyVelocity[0];
            departureState[5] = velocityAfterDeparture[1];// + departureBodyVelocity[1];
            departureState[6] = velocityAfterDeparture[2];// + departureBodyVelocity[2];
            departureState[7] = thisDV; // DV added to file just to be helpful

//            std::cout << "pos1" << departureState[0] << std::endl;

            thisStateMap.insert(std::pair<int, Eigen::VectorXd> (j, departureState));




//                        std::cout << "Saving le trajectory " << std::endl;
            // Define vectors to calculate intermediate points
            std::vector< Eigen::Vector3d > intermediatePositionVector;
            std::vector< double > intermediateTimeVector;
            thisTrajectory.intermediatePoints(0.1*physical_constants::JULIAN_YEAR, intermediatePositionVector, intermediateTimeVector);
////    std::cout << "One of the position vector values: " << intermediatePositionVector[0] << std::endl;
//    std::string saveFilename = "/home/matt/tudatBundle/tudatApplications/EDT_thesis/tudat_apps/main_app/THIS.dat";
//
////    std::string saveFilename = "/this.dat";
            std::string outputPathBase = "/home/matt/tudatBundle/tudatApplications/EDT_thesis/tudat_apps/main_app/SimulationOutput/" + outputSubFolder;
            std::string saveFilename = outputPathBase + "/THIS.dat";
//            std::cout << "output subfolder: " << outputSubFolder_ << std::endl;
//            std::cout << "save file name: " << saveFilename << std::endl;
            tudat::transfer_trajectories::writeTrajectoryToFile( intermediatePositionVector, intermediateTimeVector, saveFilename );




        }
        std::string saveFilename = "/home/matt/tudatBundle/tudatApplications/EDT_thesis/tudat_apps/main_app/THIS.dat";
        std::cout << "Saving " << i <<  " to file" << std::endl;

        std::string outputPath = "/home/matt/tudatBundle/tudatApplications/EDT_thesis/tudat_apps/main_app/SimulationOutput/" + outputSubFolder;
        std::string filename = "initialState_" + fileSuffix + std::to_string(i) + ".dat";
//        std::cout << outputPath << std::endl;
        // Write dependent variables to file
        input_output::writeDataMapToTextFile(thisStateMap,
                                             filename,
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");











    }









//    //////////////////////////////////////// Do the non-perturbed patched-conic case /////////////////////////////////////
//
//    std::cout << "\n Beginning the second stage part \n" << std::endl;
//
//    // Create body map.
//    std::vector< double > gravitationalParametersTransferBodies;
//    gravitationalParametersTransferBodies.push_back( simulationVariables["constants"]["gravitationalParameters"]["EarthMoon"] );
//    gravitationalParametersTransferBodies.push_back( simulationVariables["constants"]["gravitationalParameters"]["Jupiter"] );
//
//    std::vector< ephemerides::EphemerisPointer > ephemerisVectorTransferBodies;
//    ephemerisVectorTransferBodies.push_back( std::make_shared< ephemerides::ApproximatePlanetPositions>(
//            ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter) ) ;
//    ephemerisVectorTransferBodies.push_back( std::make_shared< ephemerides::ApproximatePlanetPositions>(
//            ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter) ) ;
//
//    // Define central body of trajectory, and the body to be propagated TODO: make sure body is right
//    std::vector< std::string > centralBody;
//    centralBody.push_back("Sun");
//    std::string bodyToPropagate = "Vehicle";
//
//    // Create and define the transfer body trajectory, using the previously defined trajectory sequence
//    std::vector< std::string > transferBodyTrajectory = gen::FlybySequence2TransferBodyTrajectory(FlybySequence);
//
//
//    simulation_setup::NamedBodyMap bodyMap = propagators::setupBodyMapFromUserDefinedEphemeridesForPatchedConicsTrajectory(
//            centralBody[0], bodyToPropagate, transferBodyTrajectory,
//            ephemerisVectorTransferBodies, gravitationalParametersTransferBodies, "ECLIPJ2000");
//
//    // Create acceleration map. NOTE: acceleration map corresponds to patched conic assumptions exactly
//    std::vector< basic_astrodynamics::AccelerationMap > patchedConicAccelerationMap = propagators::setupAccelerationMapPatchedConicsTrajectory(
//            transferBodyTrajectory.size(), centralBody[0], bodyToPropagate, bodyMap);
//
//    // Define a bunch of things from the original problem struct
//    std::vector< TransferLegType > legTypeVector = earthPlanetProbStruct.getLegTypeVector();
//    std::vector< double>  semiMajorAxes = gen::EigenVectorXd2StdVectorDoubles(earthPlanetProbStruct.getSemiMajorAxes());
//    std::vector< double> eccentricities = gen::EigenVectorXd2StdVectorDoubles(earthPlanetProbStruct.getEccentricities());
//    std::vector< double> minimumPericenterRadii = gen::EigenVectorXd2StdVectorDoubles(earthPlanetProbStruct.getMinimumPericenterRadii());
//
//
//
//
//    // For loop to run over all the islands
//    for (int j=0; j<islandsList.size(); j++) {
//        island thisIsland = islandsList[j];
//        std::vector<vector_double> theseVariableVectorYears = thisIsland.get_population().get_x();
//
//
//        // Define map to save this island with
//        std::map<double, Eigen::VectorXd> thisIslandDataMap;
//
//        // For loop to run over all of the variable vectors in this island
//        for (int i = 0; i < theseVariableVectorYears.size(); i++) {
//            std::vector<double> variableVectorYears = theseVariableVectorYears[i];
////            std::cout << "VARIABLE VECTOR: " << std::endl;
////
////
////
////            for (std::vector<double>::const_iterator i = variableVectorYears.begin();
////                 i != variableVectorYears.end(); ++i)
////                std::cout << *i << ' ';
//
//            std::vector<double> variableVector = gen::vectorScaling(variableVectorYears,
//                                                                    physical_constants::JULIAN_YEAR); // TODO: maybe rename to secs?
//            variableVector.push_back(1); // Add dummy variable for TOF of capture leg
//
//            // Define integrator settings. TODO: possibly make integrator setting specific to this, currently just uses same as main sim
//            nlohmann::json integratorSettingsJson = simulationVariables["GuidanceConfigs"]["integratorSettings"];
//            std::string RKCoefficients = integratorSettingsJson["integrator"];
//            RungeKuttaCoefficients::CoefficientSets coefficientSet = gen::getIntegratorCoefficientSet(RKCoefficients);
//            double initialEphemerisTime = variableVector[0];
//
//
////            std::cout << "This is before integrator settings" << std::endl;
//
//            // NOTE: These integrator settings don't make much sense for an actual simulation, but this is done only to get the initial conditions so
//            std::shared_ptr<IntegratorSettings<> > integratorSettings =
//                    std::make_shared<RungeKuttaVariableStepSizeSettings<> >(
//                            initialEphemerisTime,
//                            integratorSettingsJson["initialTimeStep"],
//                            coefficientSet,
//                            integratorSettingsJson["minimumStepSize"],
//                            integratorSettingsJson["maximumStepSize"],
//                            0.5,
//                            0.5);
////            std::shared_ptr<IntegratorSettings<> > integratorSettings =
////                    std::make_shared<RungeKuttaVariableStepSizeSettings<> >(
////                            initialEphemerisTime,
////                            integratorSettingsJson["initialTimeStep"],
////                            coefficientSet,
////                            integratorSettingsJson["minimumStepSize"],
////                            integratorSettingsJson["maximumStepSize"],
////                            integratorSettingsJson["relativeErrorTolerance"],
////                            integratorSettingsJson["absoluteErrorTolerance"]);
//
//
//
//            // Calculate patched conic solution and propagation results of full dynamics, for each leg. Is done by creating trajectory maps and filling them with the relevant function
//            std::map<int, std::map<double, Eigen::Vector6d> > patchedConicsTrajectory;
//            std::map<int, std::map<double, Eigen::Vector6d> > fullProblemPatchedConicTrajectory;
//            std::map<int, std::map<double, Eigen::VectorXd> > dependentVariablesResultForEachLeg;
//
//            // Define a list of dependent variables, and create dependent variable save settings
//            std::vector<std::shared_ptr<SingleDependentVariableSaveSettings> > dependentVariablesList;
//
//            dependentVariablesList.push_back(
//                    std::make_shared<SingleDependentVariableSaveSettings>(relative_position_dependent_variable, "Earth",
//                                                                          "Sun")
//            );
//            dependentVariablesList.push_back(
//                    std::make_shared<SingleDependentVariableSaveSettings>(relative_position_dependent_variable,
//                                                                          planetToFlyby, "Sun")
//            );
////    dependentVariablesList.push_back(
////            std::make_shared< SingleDependentVariableSaveSettings >(relative_position_dependent_variable, "Saturn", "Sun")
////    ); TODO: add me back in later, need to add Saturn to acceleration map
//
//            std::shared_ptr<DependentVariableSaveSettings> dependentVariablesToSave =
//                    std::make_shared<DependentVariableSaveSettings>(dependentVariablesList,  false);
//
////            std::cout << "This is before propagation" << std::endl;
//
//
//
//            propagators::fullPropagationPatchedConicsTrajectory(bodyMap,
//                                                                patchedConicAccelerationMap[0],
//                                                                transferBodyTrajectory,
//                                                                centralBody[0],
//                                                                bodyToPropagate,
//                                                                legTypeVector,
//                                                                variableVector,
//                                                                minimumPericenterRadii,
//                                                                semiMajorAxes,
//                                                                eccentricities,
//                                                                integratorSettings,
//                                                                patchedConicsTrajectory,
//                                                                fullProblemPatchedConicTrajectory,
//                                                                dependentVariablesResultForEachLeg,
//                                                                true,
//                                                                dependentVariablesToSave,
//                                                                cowell);
//
////            patchedConicsTrajectory
////            std::cout << "This is before adding to map" << std::endl;
//
//            Eigen::Vector6d thisDepartureVector;
//            double thisDepartureTimeYears;
//
//            for (auto itr : patchedConicsTrajectory) {
//                thisDepartureVector = patchedConicsTrajectory[itr.first].begin()->second;
//                thisDepartureTimeYears = (patchedConicsTrajectory[itr.first].begin()->first) / physical_constants::JULIAN_YEAR;
////                std::cout << "Departure time years: " << thisDepartureTimeYears << std::endl;
////                std::cout << "Departure time V2: " << (patchedConicsTrajectory[0][0]) / physical_constants::JULIAN_YEAR << std::endl;
//
//            }
//
//            thisIslandDataMap.insert(std::pair<double, Eigen::VectorXd>(thisDepartureTimeYears, thisDepartureVector));
//
//
//
//
//
//        // Do some printing during simulation. First does differences, second does raw
////        for (auto itr : patchedConicsTrajectory) {
////            std::cout << "Leg " << itr.first << "\n\n";
////            std::cout << "Departure difference: " << fullProblemPatchedConicTrajectory[itr.first].begin()->second -
////                                                     patchedConicsTrajectory[itr.first].begin()->second << "\n\n";
////            std::cout << "Arrival difference: " << fullProblemPatchedConicTrajectory[itr.first].rbegin()->second -
////                                                   patchedConicsTrajectory[itr.first].rbegin()->second << "\n\n";
////        }
////        for (auto itr : patchedConicsTrajectory) {
//////            std::cout << "Leg " << itr.first << "\n\n";
//////            Eigen::Vector6d departureVector = fullProblemPatchedConicTrajectory[itr.first].begin()->second;
//////            std::cout << "Departure state vector: " << departureVector
//////                      << "\n\n";
////////            std::cout << "Departure state vector: " << fullProblemPatchedConicTrajectory[itr.first].begin()->second
////////                      << "\n\n";
////            double departureTime = fullProblemPatchedConicTrajectory[itr.first].begin()->first;
//////            std::cout << "Departure time: " << departureTime << "\n\n";
////            std::cout << "Departure time years: " << departureTime / 365.25 / 24 / 60 / 60 << std::endl;
////
//////            std::cout << "Arrival state vector: " << fullProblemPatchedConicTrajectory[itr.first].rbegin()->second
//////                      << "\n\n";
////        }
////        for( auto itr2 : fullProblemPatchedConicTrajectory[ itr.first ]){
////            std::cout << itr2.first << std::endl;
////        }
//
//
//
//        }
//        std::cout << "data map size: " << thisIslandDataMap.size() << std::endl;
//        std::string islandDataMapSavename = "initialState_GA_EJ_" + std::to_string(j) + ".dat";
//        std::cout << "Saving " << islandDataMapSavename << std::endl;
//        // Writes data to file for this variable vector
//        input_output::writeDataMapToTextFile(thisIslandDataMap,
//                                             islandDataMapSavename,
//                                             tudat_applications::getOutputPath() + outputSubFolder,
//                                             "",
//                                             std::numeric_limits<double>::digits10,
//                                             std::numeric_limits<double>::digits10,
//                                             ",");
//
//    }
//
//
//
//    // TODO: Reinstate below when / if actually needed
////    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    //////////////// Do the perturbed  case using initial conditions from patched conic   //////////////////////////////
////    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////    ////////////////////////   Get initial conditions from patchedConicsTrajectory    //////////////////////////////////
////    double departureEpoch;
////    Eigen::Vector6d departureStateVector;
////    double rawDepartureEpoch;
////    Eigen::Vector6d rawDepartureStateVector;
////    double fullSimDelay = 14 * 24*60*60;
////    for( auto itr : patchedConicsTrajectory )
////    {
////        if (itr.first == 0){
////            // Get departure epoch and state vector
////            rawDepartureEpoch = patchedConicsTrajectory[ itr.first ].begin( )->first;
////            rawDepartureStateVector = patchedConicsTrajectory[ itr.first ].begin( )->second;
////            std::cout<<"Raw departure epoch, years: " << rawDepartureEpoch /60 /60 /24 /365.25 << std::endl;
////
////            // Implement to delay to starting initial simulation
////            bool departureEpochFound = false;
////            for( auto itr2 : fullProblemPatchedConicTrajectory[ itr.first ]){
//////                std::cout<<itr2.first<<std::endl;
//////                std::cout << itr2.first << std::endl;
//////                std::cout << rawDepartureEpoch + fullSimDelay << "\n\n" << std::endl;
////                if ( (itr2.first > rawDepartureEpoch + fullSimDelay) and (departureEpochFound == false)){
//////                    std::cout << "here" << std::endl;
////                    departureEpoch = itr2.first;
////                    departureStateVector = itr2.second;
////                    departureEpochFound = true;
////                    std::cout<<"Raw departure epoch, years: " << rawDepartureEpoch /60 /60 /24 /365.25 << std::endl;
////                }
////            }
////
////        }
////    }
////
////    ////////////////////////// Create propagation mechanism in the conventional way ////////////////////////////////
////    // NOTE: uses the methods created for the EDT simulations, for consistency, but does not require full functionlaity
////    // therefore some dummy variables etc are created
////
////    // Create EDTConfig for GA
////    std::cout<< " -- Creating Config class -- " << std::endl;
////    nlohmann::json configVariables = dummyVariablesJson["EDTConfigs"];
////    std::string configType = dummyVariablesJson["EDTConfigs"]["configType"];
////    nlohmann::json SRPVariables = dummyVariablesJson["scConfigs"]["SRP"];
////    nlohmann::json materialProperties = dummyVariablesJson["materialProperties"];
////    EDTs::EDTConfig GAConfig = EDTs::EDTConfig(configVariables, configType, SRPVariables, materialProperties);
////
////    // Create EDTEnvironment for GA, using some dummy variables
////    std::cout<< " -- Creating environment class -- " << std::endl;
////    NamedBodyMap GABodyMap;
////    double dummyphi0;
////    double dummyR0;
////    std::vector<double> dummytwoSinePars;
//////    nlohmann::json dummyISMFVariables = dummyVariablesJson["InterstellarMagField"];
////    EDTEnvironment GAEnvironment = EDTEnvironment(dummytwoSinePars, dummyphi0, dummyR0, GABodyMap, simulationVariables);
////
////    // Create EDTGuidance for GA
////    std::cout<< " -- Creating Guidance class -- " << std::endl;
////    std::string placeholderThrustMagnitudeConfig = dummyVariablesJson["GuidanceConfigs"]["thrustMagnitudeConfig"];
////    std::string placeholderThrustDirectionConfig = dummyVariablesJson["GuidanceConfigs"]["thrustDirectionConfig"];
////    EDTGuidance GAGuidance = EDTGuidance(
////            placeholderThrustMagnitudeConfig,
////            placeholderThrustDirectionConfig,
////            GABodyMap,
////            GAEnvironment,
////            GAConfig);
////
////    // Create propBodies class for GA
////    std::cout<< " -- Creating Propbodies class -- " << std::endl;
//////    nlohmann::json jsonBodiesToInclude = simulationVariables["Spice"]["bodiesToInclude"]; // TODO: Check if these should be used
////    univ::propBodies GAPropBodies = univ::propBodies(GAConfig, GAGuidance, GABodyMap, simulationVariables, false);
////
////    // Create propSettings class for GA
////    std::cout<< " -- Creating Propsettings class -- " << std::endl;
////    Eigen::Vector3d departurePosition;
////    departurePosition << departureStateVector[0], departureStateVector[1], departureStateVector[2];
////    Eigen::Vector3d departureVelocity;
////    departureVelocity << departureStateVector[3], departureStateVector[4], departureStateVector[5];
////
////    // Set json info
////    nlohmann::json terminationSettingsJson = simulationVariables["GuidanceConfigs"]["terminationSettings"];
////    nlohmann::json integratorSettingsJsonPerturbed;
////    double simulationTimeUpperLimit = terminationSettingsJson["timeTerminationYears"];
////
////    univ::propSettings GAPropSettings = univ::propSettings(GAPropBodies,
////                                                           departurePosition,
////                                                           departureVelocity,
////                                                            integratorSettingsJson,
////                                                            terminationSettingsJson,
////                                                            departureEpoch);
////
////    std::cout << "Done creating GAPropSettings" << std::endl;
////    // Ensure environment is properly updated
//////    GAGuidance.updateAllEnviro();
////
////    std::cout<< "====================Propagating Perturbed Case==========================" << std::endl;
////    std::cout << "Departure epoch (years): " << departureEpoch /60 /60 /24 /365.25 << std::endl;
////    std::cout << "Sim time upper limit: " << simulationTimeUpperLimit << std::endl;
////    // Create simulation object and propagate dynamics.
////    SingleArcDynamicsSimulator< > dynamicsSimulator(
////            GAPropBodies.bodyMap, GAPropSettings.integratorSettings, GAPropSettings.propagatorSettings);
////    std::map< double, Eigen::VectorXd > integrationResultPerturbedCase = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
////    std::map< double, Eigen::VectorXd > dependentVariableResultPerturbedCase = dynamicsSimulator.getDependentVariableHistory();
////
////
////
////
////    //////////////////////////////////////// Save variables and maps ////////////////////////////////////////////////////
////
////    std::map< double, Eigen::VectorXd > dependentVariableResultPatchedConic = dependentVariablesResultForEachLeg[0];
//////    std::map< double, Eigen::VectorXd > dependentVariableResultPerturbed = dependentVariablesPerturbedCase[0]; TODO: uncomment me
////
////    // Save stuff to file
//////    std::string outputSubFolder = PlanetConfigs["saveDataInfo"]["outputSubFolder"]; //TODO: combine / differentiate this from the one in the first stage
////    // TODO: make file naming convention work for jupiter / saturn
////    for( std::map< int, std::map< double, Eigen::Vector6d > >::iterator itr = patchedConicsTrajectory.begin( );
////         itr != patchedConicsTrajectory.end( ); itr++ ) {
////
////        // Saving for the unperturbed patched conic case
////        input_output::writeDataMapToTextFile( patchedConicsTrajectory[itr->first],
////                                              "unperturbed_patchedConic_leg_" + std::to_string(itr->first) + ".dat",
////                                              tudat_applications::getOutputPath( ) + outputSubFolder,
////                                              "",
////                                              std::numeric_limits< double >::digits10,
////                                              std::numeric_limits< double >::digits10,
////                                              "," );
////
////        input_output::writeDataMapToTextFile( fullProblemPatchedConicTrajectory[itr->first],
////                                              "unperturbed_fullProblem_leg_" + std::to_string(itr->first) + ".dat",
////                                              tudat_applications::getOutputPath( ) + outputSubFolder,
////                                              "",
////                                              std::numeric_limits< double >::digits10,
////                                              std::numeric_limits< double >::digits10,
////                                              "," );
////
////        input_output::writeDataMapToTextFile( dependentVariableResultPatchedConic,
////                                              "unperturbed_depVars_leg_" + std::to_string(itr->first) + ".dat",
////                                              tudat_applications::getOutputPath( ) + outputSubFolder,
////                                              "",
////                                              std::numeric_limits< double >::digits10,
////                                              std::numeric_limits< double >::digits10,
////                                              "," );
////
////        // Saving for the perturbed case
////        input_output::writeDataMapToTextFile( integrationResultPerturbedCase,
////                                              "perturbed_fullProblem_leg_0.dat",
////                                              tudat_applications::getOutputPath( ) + outputSubFolder,
////                                              "",
////                                             std::numeric_limits< double >::digits10,
////                                             std::numeric_limits< double >::digits10,
////                                             "," );
////        input_output::writeDataMapToTextFile( integrationResultPerturbedCase,
////                                              "perturbed_depVars_leg_0.dat",
////                                              tudat_applications::getOutputPath( ) + outputSubFolder,
////                                              "",
////                                              std::numeric_limits< double >::digits10,
////                                              std::numeric_limits< double >::digits10,
////                                              "," );
////    }
////
////
////
////
////
////
////
    return 0;
}