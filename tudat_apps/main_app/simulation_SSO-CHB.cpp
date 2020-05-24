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
    using namespace tudat::orbital_element_conversions;


    using namespace univ;
    using namespace EDTs;
    using namespace gen;



//    std::cout << "TudatBundleRootPathBoost:  " << gen::tudatBundleRootPathBoost << std::endl;
//    std::cout << "TudatBundleRootPathTemp:  " << gen::tudatBundleRootPathTemp << std::endl;
//    std::cout << "TudatBundleRootPath:  " << gen::tudatBundleRootPath << std::endl;
//    std::cout << "TudatRootPath:  " << gen::tudatRootPath << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             SIMULATION PREP            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout<< "===============Loading Variables From File==================" << std::endl;

    nlohmann::json simulationVariables = gen::readJson("testVariables.json");


    std::cout<< "===============Prepping Sim==================" << std::endl;
    std::cout<< "Loading spice kernels" << std::endl;

    // Get spice path for merged spice kernels and load them
    std::string customKernelName = simulationVariables["Spice"]["customKernelName"];
    std::string customSpiceKernelPath = gen::tudatKernelsRootPath + customKernelName;
    const std::vector<std::string> spicePathVector = {customSpiceKernelPath};
    spice_interface::loadStandardSpiceKernels( spicePathVector );

    // Create EDT Environment class
    std::cout<< "Creating environment class" << std::endl;
    // Create base body map to be built on by other classes + give initial numbers for magfield data
    NamedBodyMap baseBodyMap;

    // Load parker variables from json then convert to proper units
    double B0 = simulationVariables["ParkerMagField"]["B0_tesla"]; // Note: B0 not currently used since is calculated using sine relationship
    double phi0_deg = simulationVariables["ParkerMagField"]["phi0_deg"];
    double R0_au = simulationVariables["ParkerMagField"]["R0_au"];
    double phi0 = phi0_deg * deg2rad;
    double R0 = R0_au * AU;

    std::vector<double> twoSinePars;
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["a1"]);
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["b1"]);
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["c1"]);
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["a2"]);
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["b2"]);
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["c2"]);
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["d"]);

    nlohmann::json ISMFVariables = simulationVariables["InterstellarMagField"];

    EDTEnvironment CHBEDTEnviro = EDTEnvironment(twoSinePars, phi0, R0, baseBodyMap, ISMFVariables);

    // Create EDT Guidance class
    std::cout<< "Creating Guidance class" << std::endl;
    std::string thrustMagnitudeConfig = simulationVariables["GuidanceConfigs"]["thrustMagnitudeConfig"];
    std::string thrustDirectionConfig = simulationVariables["GuidanceConfigs"]["thrustDirectionConfig"];

    EDTGuidance CHBEDTGuidance = EDTGuidance(
            thrustMagnitudeConfig,
            thrustDirectionConfig,
            baseBodyMap,
            CHBEDTEnviro);

    // Create EDT config class and set constant thrust in guidance class
    std::cout<< "Creating Config class" << std::endl;
    std::string configType = simulationVariables["EDTConfigs"]["configType"];
    double tetherLength = simulationVariables["EDTConfigs"]["tetherLength"];
    double tetherDiameterInner = simulationVariables["EDTConfigs"]["tetherDiameterInner"];
    double tetherDiameterOuter = simulationVariables["EDTConfigs"]["tetherDiameterOuter"];
    nlohmann::json hoytetherVariables = simulationVariables["EDTConfigs"]["hoytether"];
    nlohmann::json SRPVariables = simulationVariables["scConfigs"]["SRP"];
    EDTs::EDTConfig CHBEDTConfig = EDTs::EDTConfig(CHBEDTGuidance, configType, hoytetherVariables, SRPVariables, tetherLength, tetherDiameterInner, tetherDiameterOuter);
    CHBEDTGuidance.setThrustMagnitudeConstant(CHBEDTConfig.getConstantThrust());

    // Get universal class for propagation bodies
    std::cout<< "Creating Propbodies class" << std::endl;
    univ::propBodies SSOPropBodies = univ::propBodies(CHBEDTConfig, baseBodyMap);

    // Get universal class for propagation settings + set vehicle initial state
    std::cout<< "Creating Propsettings class" << std::endl;

    // Set vehicle initial state using keplerian elements (load from json first and convert to proper value)
    double a_au = simulationVariables["GuidanceConfigs"]["vehicleInitialKeplerian"]["a_au"];
    double e = simulationVariables["GuidanceConfigs"]["vehicleInitialKeplerian"]["e"];
    double i_deg = simulationVariables["GuidanceConfigs"]["vehicleInitialKeplerian"]["i_deg"];
    double aop_deg = simulationVariables["GuidanceConfigs"]["vehicleInitialKeplerian"]["aop_deg"];
    double raan_deg = simulationVariables["GuidanceConfigs"]["vehicleInitialKeplerian"]["raan_deg"];
    double ta_deg = simulationVariables["GuidanceConfigs"]["vehicleInitialKeplerian"]["ta_deg"];

    Eigen::Vector6d vehicleInitialKeplerian;
    vehicleInitialKeplerian << a_au*AU, e, i_deg*deg2rad, aop_deg*deg2rad, raan_deg*deg2rad, ta_deg*deg2rad;

    // Convert to cartesian state, and create position and velocity vectors
    double solarGravPar = SSOPropBodies.bodyMap["Sun"]->getGravityFieldModel()->getGravitationalParameter();
    Eigen::Vector6d vehicleInitialCartesian = convertKeplerianToCartesianElements(
            vehicleInitialKeplerian,
            solarGravPar);

    Eigen::Vector3d vehicleInitialPosition;
    vehicleInitialPosition << vehicleInitialCartesian[0], vehicleInitialCartesian[1], vehicleInitialCartesian[2];
    Eigen::Vector3d vehicleInitialVelocity;
    vehicleInitialVelocity << vehicleInitialCartesian[3], vehicleInitialCartesian[4], vehicleInitialCartesian[5];

//    std::cout << "Solar gravitational parameter: " << solarGravPar << std::endl;
//    std::cout << "Vehicle initial state Keplerian: " << vehicleInitialKeplerian << std::endl;
//    std::cout << "Vehicle initial state Cartesian: " << vehicleInitialCartesian << std::endl;
//    std::cout << "Vehicle initial position: " << vehicleInitialPosition << std::endl;
//    std::cout << "Vehicle initial velocity: " << vehicleInitialVelocity << std::endl;

    // Actually create the propSettings class
    double initialEphemerisTime = simulationVariables["GuidanceConfigs"]["initialEphemerisTime"];
    std::string terminationType = simulationVariables["GuidanceConfigs"]["terminationType"];
    double simulationTimeYears = simulationVariables["GuidanceConfigs"]["simulationTimeYears"];
    univ::propSettings SSOPropSettings = univ::propSettings(SSOPropBodies,
                                                            vehicleInitialPosition,
                                                            vehicleInitialVelocity,
                                                            initialEphemerisTime,
                                                            terminationType,
                                                            simulationTimeYears);
//    univ::propSettings SSOPropSettings = univ::propSettings(SSOPropBodies,
//                                                            {1.496E11, 0, 0},
//                                                            {0, 29.78E3, 0},
//                                                            1.0E7,
//                                                            "nominalTimeTermination",
//                                                            0.1);



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

    // Set output folder and base data title for all data
    std::string outputSubFolder = simulationVariables["saveDataConfigs"]["outputSubFolder"];
    std::string baseFilename = simulationVariables["saveDataConfigs"]["baseFilename"];

    // Write satellite propagation history to file.
    input_output::writeDataMapToTextFile( integrationResult,
                                          baseFilename + "propData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write dependent variables to file
    input_output::writeDataMapToTextFile( dependentVariableResult,
                                          baseFilename + "depVarData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write magField history to file.
    input_output::writeDataMapToTextFile( CHBEDTGuidance.getMagFieldMap(),
                                          baseFilename + "magData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write ionosphere history to file.
    input_output::writeDataMapToTextFile( CHBEDTGuidance.getIonosphereMap(),
                                          baseFilename + "ionoData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write thrust history to file.
    input_output::writeDataMapToTextFile( CHBEDTGuidance.getThrustMap(),
                                          baseFilename + "thrustData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write current history to file.
    input_output::writeDataMapToTextFile( CHBEDTGuidance.getCurrentMap(),
                                          baseFilename + "currentData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write bodyData history to files
    input_output::writeDataMapToTextFile( CHBEDTGuidance.getBodyDataMap(),
                                          baseFilename + "bodyData.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );
//

    std::cout<< "============= TESTING ===============" << std::endl;

    double tehTime = gen::tudatTime2DecimalYear(1E7);
    std::cout << "teh time is: " << tehTime << std::endl;

//    std::string pathToJson = gen::jsonInputsRootPath + "test1.json";
//    std::ifstream people_file(pathToJson);
//    nlohmann::json test1 = nlohmann::json::parse(people_file);
////    test1 << people_file;
//
//    std::cout << test1 << std::endl;
//    std::cout << test1["Anna"] << std::endl;
//    std::cout << test1["Anna"]["Profession"] << std::endl;
//
//    double annaAge = test1["Anna"]["Age"];
//    std::string annaProfession = test1["Anna"]["Profession"];




//    Eigen::Vector3d localMag;
//    Eigen::Vector3d inertialMag;
//    Eigen::Vector3d localMag2_pos;
//    Eigen::Vector3d localMag2_neg;
//    double theta = 0*deg2rad;
//    localMag << -0.06883, 0.098298, 0;
//
//    inertialMag = gen::LvlhToInertial(localMag, theta);
//    localMag2_pos = gen::LvlhToInertial(inertialMag, theta);
//    localMag2_neg = gen::LvlhToInertial(inertialMag, -theta);
//
//    std::cout << "localMag     : " << localMag << std::endl;
//    std::cout << "inertialMag  : " << inertialMag << std::endl;
//    std::cout << "localMag2_pos: " << localMag2_pos << std::endl;
//    std::cout << "localMag2_neg: " << localMag2_neg << std::endl;

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