// Main file for SSO-CHB simulation run

///////// Start by defining the SSO simulation ///////////////
//TODO: anywwhere an uncrecognised thing is done (eg magfield type) put a termination in there
#include "simulation_SSO-CHB.h"
//#include "EDT_configs.h"
#include "EDTGuidance.h"
#include "universal_settings.h"
#include "environment_settings.h"
//#include <experimental/filesystem>

int main(int argc, char *argv[] )
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


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             SIMULATION PREP            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout<< "===============Loading Variables From File==================" << std::endl;

    // Read json file for general use (name of json is passed as first argument, or default set here)
    std::string jsonName;
    if (argc != 1){
        jsonName = argv[1];
    }
    else{
//        jsonName = "VnV/testVariablesParkerVnV.json"; // TODO: change to nominal test variables
        jsonName = "testVariables.json";
    }
    nlohmann::json simulationVariables = gen::readJson(jsonName);


    std::cout<< "===============Prepping Sim==================" << std::endl;
    std::cout<< " -- Loading spice kernels -- " << std::endl;

    // Get spice path for merged spice kernels and load them
    std::string customKernelName = simulationVariables["Spice"]["customKernelName"];
    std::string customSpiceKernelPath = gen::tudatKernelsRootPath + customKernelName;
    const std::vector<std::string> spicePathVector = {customSpiceKernelPath};
    spice_interface::loadStandardSpiceKernels( spicePathVector );

    // Create EDT Environment class
    std::cout<< " -- Creating environment class -- " << std::endl;
    // Create base body map to be built on by other classes + give initial numbers for magfield data
    NamedBodyMap baseBodyMap;

    // Load parker variables from json then convert to proper units
    double phi0_deg = simulationVariables["ParkerMagField"]["phi0_deg"];
    double R0_au = simulationVariables["ParkerMagField"]["R0_au"];
    double phi0 = phi0_deg * deg2rad;
    double R0 = R0_au * AU;
    std::string configType = simulationVariables["EDTConfigs"]["configType"];

    std::vector<double> twoSinePars;
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["a1"]);
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["b1"]);
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["c1"]);
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["a2"]);
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["b2"]);
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["c2"]);
    twoSinePars.push_back(simulationVariables["ParkerMagField"]["twoSinePars"]["d"]);

    nlohmann::json ISMFVariables = simulationVariables["InterstellarMagField"];
    double testYear = 2001.3;
    std::cout << "B0 for year " << testYear << " : " << gen::twoSines(testYear, twoSinePars) << std::endl;
    EDTEnvironment CHBEDTEnviro = EDTEnvironment(twoSinePars, phi0, R0, baseBodyMap, ISMFVariables);

    // Create EDT config class and set constant thrust in guidance class
    std::cout<< " -- Creating Config class -- " << std::endl;
    nlohmann::json configVariables = simulationVariables["EDTConfigs"];
    nlohmann::json SRPVariables = simulationVariables["scConfigs"]["SRP"];
    nlohmann::json materialProperties = simulationVariables["materialProperties"];
    EDTs::EDTConfig CHBEDTConfig = EDTs::EDTConfig(configVariables, configType, SRPVariables, materialProperties);


    // Create EDT Guidance class
    std::cout<< " -- Creating Guidance class -- " << std::endl;
    std::string thrustMagnitudeConfig = simulationVariables["GuidanceConfigs"]["thrustMagnitudeConfig"];
    std::string thrustDirectionConfig = simulationVariables["GuidanceConfigs"]["thrustDirectionConfig"];

    EDTGuidance CHBEDTGuidance = EDTGuidance(
            thrustMagnitudeConfig,
            thrustDirectionConfig,
            baseBodyMap,
            CHBEDTEnviro,
            CHBEDTConfig);
    CHBEDTGuidance.setThrustMagnitudeConstant(CHBEDTConfig.getConstantThrust()); // TODO: Check if this works / is needed


    // Get universal class for propagation bodies
    std::cout<< " -- Creating Propbodies class -- " << std::endl;
//    nlohmann::json jsonBodiesToInclude = simulationVariables["Spice"]["bodiesToInclude"];
    univ::propBodies SSOPropBodies = univ::propBodies(CHBEDTConfig, CHBEDTGuidance, baseBodyMap, simulationVariables, true);

    // Get universal class for propagation settings + set vehicle initial state
    std::cout<< " -- Creating Propsettings class -- " << std::endl;

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

    std::cout << "Vehicle initial cartesian components: "
              << vehicleInitialCartesian[0] << " "
              << vehicleInitialCartesian[1] << " "
              << vehicleInitialCartesian[2] << " "
              << vehicleInitialCartesian[3] << " "
              << vehicleInitialCartesian[4] << " "
              << vehicleInitialCartesian[5]
              << std::endl; //TODO:Remove me

    Eigen::Vector3d vehicleInitialPosition;
    vehicleInitialPosition << vehicleInitialCartesian[0], vehicleInitialCartesian[1], vehicleInitialCartesian[2];
    Eigen::Vector3d vehicleInitialVelocity;
    vehicleInitialVelocity << vehicleInitialCartesian[3], vehicleInitialCartesian[4], vehicleInitialCartesian[5];

    // Set various json variables
    double initialEphemerisYear = simulationVariables["GuidanceConfigs"]["initialEphemerisYear"];
    nlohmann::json  terminationSettingsJson = simulationVariables["GuidanceConfigs"]["terminationSettings"];
    nlohmann::json integratorSettingsJson = simulationVariables["GuidanceConfigs"]["integratorSettings"];

    // Actually create the propSettings class
    univ::propSettings SSOPropSettings = univ::propSettings(SSOPropBodies,
                                                            vehicleInitialPosition,
                                                            vehicleInitialVelocity,
                                                            integratorSettingsJson,
                                                            terminationSettingsJson,
                                                            initialEphemerisYear);



    // Ensure environment is properly updated
    CHBEDTGuidance.updateAllEnviro();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "====================Propagating==========================" << std::endl;

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
            SSOPropBodies.bodyMap, SSOPropSettings.integratorSettings, SSOPropSettings.propagatorSettings);
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////             SAVE DATA                  ////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "====================Saving Data==========================" << std::endl;

    // Set output folder and base data title for all data
    std::string outputSubFolder = simulationVariables["saveDataConfigs"]["outputSubFolder"];
    std::string outputPath = tudat_applications::getOutputPath() + outputSubFolder;
    std::string baseFilename = simulationVariables["saveDataConfigs"]["baseFilename"];

    // Delete all existing files in directory before saving
    boost::filesystem::remove_all(outputPath);

    // Set variables for data types to be saved
    nlohmann::json dataTypesToSave = simulationVariables["saveDataConfigs"]["dataTypesToSave"];
    bool savePropData = dataTypesToSave["propData"];
    bool saveMagneticFieldData = dataTypesToSave["magneticField"];
    bool saveIonosphereData = dataTypesToSave["ionosphere"];
    bool saveThrustData = dataTypesToSave["thrust"];
    bool saveCurrentData = dataTypesToSave["current"];
    bool saveBodyData = dataTypesToSave["bodyData"];
    bool saveDependentVariablesData = dataTypesToSave["dependentVariables"];

    // Check which data types to save, and save them

    if (savePropData) {
        std::cout << "Saving propData" << std::endl;

        // Write satellite propagation history to file.
        input_output::writeDataMapToTextFile(integrationResult,
                                             baseFilename + "propData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        std::cout << "Ignoring propData" << std::endl;
    }

    if (saveMagneticFieldData) {
        std::cout << "Saving magneticFieldData" << std::endl;

        // Write magField history to file.
        input_output::writeDataMapToTextFile(CHBEDTGuidance.getMagFieldMap(),
                                             baseFilename + "magData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        std::cout << "Ignoring magneticFieldData" << std::endl;
    }

    if (saveIonosphereData) {
        std::cout << "Saving ionosphericData" << std::endl;

        // Write ionosphere history to file.
        input_output::writeDataMapToTextFile(CHBEDTGuidance.getIonosphereMap(),
                                             baseFilename + "ionoData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        std::cout << "Ignoring ionosphericData" << std::endl;
    }

    if (saveThrustData) {
        std::cout << "Saving thrustData" << std::endl;

        // Write thrust history to file.
        input_output::writeDataMapToTextFile(CHBEDTGuidance.getThrustMap(),
                                             baseFilename + "thrustData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        std::cout << "Ignoring thrustData" << std::endl;
    }

    if (saveCurrentData) {
        std::cout << "Saving currentData" << std::endl;

        // Write current history to file.
        input_output::writeDataMapToTextFile(CHBEDTGuidance.getCurrentMap(),
                                             baseFilename + "currentData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        std::cout << "Ignoring currentData" << std::endl;
    }

    if (saveBodyData) {
        std::cout << "Saving bodyData" << std::endl;

        // Write bodyData history to files
        input_output::writeDataMapToTextFile(CHBEDTGuidance.getBodyDataMap(),
                                             baseFilename + "bodyData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        std::cout << "Ignoring bodyData" << std::endl;
    }

    if (saveDependentVariablesData) {
        std::cout << "Saving dependentVariableData" << std::endl;

        // Write dependent variables to file
        input_output::writeDataMapToTextFile(dependentVariableResult,
                                             baseFilename + "depVarData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        std::cout << "Ignoring dependentVariableData" << std::endl;
    }




    std::cout << "Program finished " << std::endl;
}