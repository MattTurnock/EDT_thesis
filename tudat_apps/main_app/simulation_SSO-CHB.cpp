// Main file for SSO-CHB simulation run

///////// Start by defining the SSO simulation ///////////////

#include "simulation_SSO-CHB.h"
#include "EDT_configs.h"
#include "EDTGuidance.h"
#include "universal_settings.h"
#include "applicationOutput.h"
#include "environment_settings.h"




int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::interpolators;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_mathematics;
    using namespace tudat::basic_astrodynamics;


    using namespace univ;
    using namespace EDTs;

    double pi = 3.14159;
    double deg2rad = pi / 180;
    double AU = 1.496e11;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             SIMULATION PREP            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout<< "===============Prepping Sim==================" << std::endl;
    std::cout<< "Loading spice kernels" << std::endl;


    //Load spice kernels (merged).
//    spice_interface::loadStandardSpiceKernels( );


    // Get spice path for merged spice kernels and load them
    std::string currentFilePath(__FILE__);
    boost::filesystem::path spicePath(currentFilePath);
    spicePath = spicePath.parent_path().parent_path().parent_path().parent_path().parent_path();
    spicePath += "/tudat/Tudat/External/SpiceInterface/Kernels/tudat_EDT_spk_kernel.bsp";

    const std::vector<std::string> spicePathVector = {spicePath.string()};
    spice_interface::loadStandardSpiceKernels( spicePathVector );


    // Create base body map to be built on by other classes + give initial numbers for magfield data
    NamedBodyMap baseBodyMap;

    double B0 = 12E-9;
    double phi0 = 55 * deg2rad;
    double R0 = 1*AU;

    // Create EDT Environment class
    std::cout<< "Creating environment class" << std::endl;
    EDTEnvironment CHBEDTEnviro = EDTEnvironment(B0, phi0, R0, baseBodyMap);

    // Create EDT Guidance class
    std::cout<< "Creating Guidance class" << std::endl;
    std::string thrustMagnitudeConfig = "nominal";
    std::string thrustDirectionConfig = "nominalPrograde";

    EDTGuidance CHBEDTGuidance = EDTGuidance(
            thrustMagnitudeConfig,
            thrustDirectionConfig,
            baseBodyMap,
            CHBEDTEnviro);

    // Create EDT config class and set constant thrust in guidance class
    std::cout<< "Creating Config class" << std::endl;
    std::string configType = "CHB";
    EDTs::EDTConfig CHBEDTConfig = EDTs::EDTConfig(CHBEDTGuidance, configType);
    CHBEDTGuidance.setThrustMagnitudeConstant(CHBEDTConfig.getConstantThrust());

    // Get universal class for propagation bodies
    std::cout<< "Creating Propbodies class" << std::endl;
    univ::propBodies SSOPropBodies = univ::propBodies(CHBEDTConfig, baseBodyMap);

    // Get universal class for propagation settings + set vehicle initial state
    std::cout<< "Creating Propsettings class" << std::endl;
    univ::propSettings SSOPropSettings = univ::propSettings(SSOPropBodies,
                                                            {1.496E11, 0, 0},
                                                            {0, 29.78E3, 0},
                                                            1.0E7,
                                                            "nominalTimeTermination",
                                                            10);

    // Ensure environment is properly updated
    CHBEDTEnviro.updateAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "====================Propagating==========================" << std::endl;

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
            SSOPropBodies.bodyMap, SSOPropSettings.integratorSettings, SSOPropSettings.propagatorSettings);
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory();


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             SAVE DATA                  ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "====================Saving Data==========================" << std::endl;

    // Set output folder for all data
    std::string outputSubFolder = "SSO-CHB-Test-custom-2/";

    // Write satellite propagation history to file.
    input_output::writeDataMapToTextFile( integrationResult,
                                          "SSO-CH-Test-out-expo-enviro__E-6_propData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write dependent variables to file
    input_output::writeDataMapToTextFile( dependentVariableResult,
                                          "SSO-CH-Test-out-expo-enviro__E-6_depVarData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write magField history to file.
    input_output::writeDataMapToTextFile( CHBEDTGuidance.getMagFieldMap(),
                                          "SSO-CH-Test-out-expo-enviro__E-6_magData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write ionosphere history to file.
    input_output::writeDataMapToTextFile( CHBEDTGuidance.getIonosphereMap(),
                                          "SSO-CH-Test-out-expo-enviro__E-6_ionoData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write thrust history to file.
    input_output::writeDataMapToTextFile( CHBEDTGuidance.getThrustMap(),
                                          "SSO-CH-Test-out-expo-enviro__E-6_thrustData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write current history to file.
    input_output::writeDataMapToTextFile( CHBEDTGuidance.getCurrentMap(),
                                          "SSO-CH-Test-out-expo-enviro__E-6_currentData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );
//

//    std::cout<< "============= TESTING ===============" << std::endl;
//
//    std::map< std::string, Eigen::VectorXd > testMap;
////    std::pair<std::string, double> testPair;
//    Eigen::VectorXd testPairVector;
//    testPairVector.resize(3);
//
////    testPairVector <<1,2,3;
//    testPairVector[0] =  1;
//    testPairVector[1] =  2;
//    testPairVector[2] =  3;
////    testPair = ("one", 1);
//    testMap.insert( std::pair<std::string, Eigen::VectorXd> ("one", testPairVector) );
//
//
//    testPairVector[0] =  4;
//    testPairVector[1] =  5;
//    testPairVector[2] =  6;
//
////    testPair = ("one", 1);
//    testMap.insert( std::pair<std::string, Eigen::VectorXd> ("two", testPairVector) );
//
////    testPairVector = {4,5,6};
////    testPair = ("two", 2);
////    testMap.insert( testPair  );
//
//
//
//    std::string outputSubFolder = "testDump/";
//    input_output::writeDataMapToTextFile( testMap,
//                                          "testMap.dat",
//                                          tudat_applications::getOutputPath( ) + outputSubFolder,
//                                          "",
//                                          std::numeric_limits< double >::digits10,
//                                          std::numeric_limits< double >::digits10,
//                                          "," );

//    Eigen::Matrix3d LvlhToPlanetocentric = tudat::reference_frames::getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(0,0);
//
//    std::cout << LvlhToPlanetocentric << std::endl;


//
//    std::string filepath_(__FILE__);
//    boost::filesystem::path boostPath(filepath_);
//    boostPath = boostPath.parent_path().parent_path().parent_path().parent_path().parent_path();
//    boostPath += "/tudat/Tudat/External/SpiceInterface/Kernels/tudat_merged_spk_kernel.inp";
//    std::cout<< boostPath <<std::endl;
}