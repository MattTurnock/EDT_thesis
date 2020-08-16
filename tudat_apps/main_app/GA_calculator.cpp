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
#include "general_functions.h"

#include "../../../../tudatExampleApplications/libraryExamples/PaGMOEx/Problems/getAlgorithm.h"
#include "../../../../tudatExampleApplications/libraryExamples/PaGMOEx/Problems/saveOptimizationResults.h"

int main()
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    using namespace tudat;
//    using namespace pagmo;

    // Read json file for general use
    nlohmann::json simulationVariables = gen::readJson("testVariables.json");
    nlohmann::json GAConfigs = simulationVariables["GAConfigs"];
    std::string planetToFlyby = GAConfigs["planetToFlyby"];

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            Do Lambert solver stage              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////// Set up the problem, with bounds //////////////////////////////////////

    // Set seed for reproducible results
    pagmo::random_device::set_seed(123);

    // Uses 3 decision variables (start dates, flight duration, delta V) each with lower and upper bounds (which are
    // implemented with a vector of vectors)
    std::vector< std::vector< double > > JupiterBounds(2, std::vector< double > (2, 0.0));
    std::vector< std::vector< double > > SaturnBounds(2, std::vector< double > (2, 0.0));
    std::vector< double > JupiterDVBounds(2);
    std::vector< double > SaturnDVBounds(2);

    // Define bounds, by using values defined in the json TODO: set reasonable values in json
    JupiterBounds[0][0] = gen::year2MJDSeconds(GAConfigs["JupiterConfigs"]["Bounds"]["StartYearLower"]);
    JupiterBounds[1][0] = gen::year2MJDSeconds( GAConfigs["JupiterConfigs"]["Bounds"]["StartYearUpper"]);
    JupiterBounds[0][1] = gen::years2Seconds(GAConfigs["JupiterConfigs"]["Bounds"]["FlightDurationYearsLower"]);
    JupiterBounds[1][1] = gen::years2Seconds(GAConfigs["JupiterConfigs"]["Bounds"]["FlightDurationYearsUpper"]);
    JupiterDVBounds[0] = gen::km2m(GAConfigs["JupiterConfigs"]["Bounds"]["DeltaVKmsLower"]);
    JupiterDVBounds[1] = gen::km2m(GAConfigs["JupiterConfigs"]["Bounds"]["DeltaVKmsUpper"]);

    SaturnBounds[0][0] = gen::year2MJDSeconds(GAConfigs["SaturnConfigs"]["Bounds"]["StartYearLower"]);
    SaturnBounds[1][0] = gen::year2MJDSeconds(GAConfigs["SaturnConfigs"]["Bounds"]["StartYearUpper"]);
    SaturnBounds[0][1] = gen::years2Seconds(GAConfigs["SaturnConfigs"]["Bounds"]["FlightDurationYearsLower"]);
    SaturnBounds[1][1] = gen::years2Seconds(GAConfigs["SaturnConfigs"]["Bounds"]["FlightDurationYearsUpper"]);
    SaturnDVBounds[0] = gen::km2m(GAConfigs["SaturnConfigs"]["Bounds"]["DeltaVKmsLower"]);
    SaturnDVBounds[1] = gen::km2m(GAConfigs["SaturnConfigs"]["Bounds"]["DeltaVKmsUpper"]);

    // Define the flyby sequences (ie Earth Saturn / Earth Jupiter)
    std::vector< int > JupiterFlybySequence;
    JupiterFlybySequence.push_back(3);
    JupiterFlybySequence.push_back(5);
    std::vector< int > SaturnFlybySequence;
    SaturnFlybySequence.push_back(3);
    SaturnFlybySequence.push_back(6);

    // Create problem object for problem fitness, for Jupiter and Saturn Cases
    pagmo::problem JupiterProb{ EarthPlanetTransfer(JupiterBounds, JupiterDVBounds, JupiterFlybySequence) };
    pagmo::problem SaturnProb{ EarthPlanetTransfer(SaturnBounds, SaturnDVBounds, SaturnFlybySequence) };

    // Select algorithm for problem, based on the one in the json TODO: expand possibilites of algorithsm that can be uesed with a custom function
    // NOTE: indices are: [0,1,2,3,4,5,6,7,8,9,10] = [de, sade, de1220, pso, sea, sga, simulated_annealing, bee_colony,
    // cmaes, ihs, xnes]. Can add more with custom get algorithm function if desired
    std::string algorithmName = GAConfigs["algorithmName"];
    int algorithmIndex;
    if (algorithmName == "nsga2"){
        algorithmIndex = 0;
    }
    else{
        std::cout << "ALGORITHM NOT RECOGNISED: terminating simulation" << std::endl;
        return 1;
    }
    algorithm algo = getMultiObjectiveAlgorithm(algorithmIndex);

    // Create islands with individuals defined by json TODO: make json numbers reasonable
    island JupiterIsl{algo, JupiterProb, GAConfigs["islandSize"]};
    island SaturnIsl{algo, SaturnProb, GAConfigs["islandSize"]};

    // Evolve for number of generations defined by json TODO: make json numbers reasonable
    int noGenerations = GAConfigs["noGenerations"];
    for( int i=0 ; i<noGenerations ; i++){

        // Evolve for Jupiter
        JupiterIsl.evolve();
        while( JupiterIsl.status( ) != pagmo::evolve_status::idle &&
                JupiterIsl.status( ) != pagmo::evolve_status::idle_error ){
            JupiterIsl.wait( );
        }
        JupiterIsl.wait_check( ); // Raises errors

        // Evolve for Saturn
        SaturnIsl.evolve();
        while( SaturnIsl.status( ) != pagmo::evolve_status::idle &&
                SaturnIsl.status( ) != pagmo::evolve_status::idle_error ){
            SaturnIsl.wait( );
        }
        SaturnIsl.wait_check( ); // Raises errors

        std::string outputSubFolder = GAConfigs["saveDataInfo"]["outputSubFolder"];
        std::string filePrefix = ""; //TODO: decide if an actual prefix needed
        // Write current iteration results to file for Jupiter NOTE: outputs as "filePrefix populatin_/fitness_ fileSuffix .dat
        gen::printPopulationToFile(JupiterIsl.get_population( ).get_x( ),filePrefix, "GA_EJ_" + std::to_string( i ), outputSubFolder, false );
        gen::printPopulationToFile(JupiterIsl.get_population( ).get_f( ),filePrefix, "GA_EJ_" + std::to_string( i ), outputSubFolder,true );

        // Write current iteration results to file for Jupiter
        gen::printPopulationToFile(SaturnIsl.get_population( ).get_x( ),filePrefix, "GA_ES_" + std::to_string( i ), outputSubFolder,false );
        gen::printPopulationToFile(SaturnIsl.get_population( ).get_f( ),filePrefix, "GA_ES_" + std::to_string( i ), outputSubFolder,true );

        // Print current generation
        std::cout<<i<<std::endl;


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            Second stage of optimisation              //////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define integrator settings. TODO: change the settings to something more reasonable (ie not fixed step)
        double initialTime = 0.0;
        double fixedStepSize = 100.0;
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
                std::make_shared < numerical_integrators::IntegratorSettings < > > ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize);




    }


    return 0;
}