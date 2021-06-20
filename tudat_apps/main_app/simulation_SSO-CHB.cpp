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


    // Read json file for general use (name of json is passed as first argument, or default set here)
    std::string jsonName;
    if (argc != 1){
        jsonName = argv[1];
    }
    else{
        jsonName = "VnV/testVariablesCurrentVnV_1a.json"; // TODO: change to nominal test variables
//        jsonName = "finalSims/SSO.json"; // TODO: change to nominal test variables
//        jsonName = "finalSims/InO.json"; // TODO: change to nominal test variables
//        jsonName = "VnV/testVariablesIntegratorGA_Reference.json"; // TODO: Make python file to run all variants as needed
//        jsonName = "testVariables.json";
//        jsonName = "finalSims/SOKGA_reference.json";
//        jsonName = "testCase_Guidance.json";
//        jsonName = "testCase_NoGuidance.json";
//        jsonName = "testCase_NoGuidanceNeg.json";
//        jsonName = "finalSims/SSO_Config_Sensitivity/SSO_Bare_Base.json";
//        jsonName = "finalSims/SSO_Config_Sensitivity/SSO_Trans_Base.json";
//        jsonName = "GATestVariablesNominal.json";
//        jsonName = "finalSims/SSO_Config_Sensitivity_Nominal.json";
//        jsonName = "finalSims/SSOP/SSOP_Temp.json";
//        jsonName = "finalSims/SOKGA_Stage2/SOKGA_Stage2_Jupiter_112.json";
    }
    nlohmann::json simulationVariables = gen::readJson(jsonName);
    bool verbosity = simulationVariables["saveDataConfigs"]["verbosity"];

    if (verbosity) {
        std::cout<< "===============Loaded Variables From File==================" << std::endl;
        std::cout << "===============Prepping Sim==================" << std::endl;
        std::cout << " -- Loading spice kernels -- " << std::endl;
    }

    // Get spice path for merged spice kernels and load them
    std::string customKernelName = simulationVariables["Spice"]["customKernelName"];
    std::string customSpiceKernelPath = gen::tudatKernelsRootPath + customKernelName;
    const std::vector<std::string> spicePathVector = {customSpiceKernelPath};
    spice_interface::loadStandardSpiceKernels( spicePathVector );

    // Create EDT Environment class
    if (verbosity) {std::cout<< " -- Creating environment class -- " << std::endl;}
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

//    nlohmann::json ISMFVariables = simulationVariables["InterstellarMagField"];
//    double testYear = 2001.3;
//    std::cout << "B0 for year " << testYear << " : " << gen::twoSines(testYear, twoSinePars) << std::endl;
    EDTEnvironment CHBEDTEnviro = EDTEnvironment(twoSinePars, phi0, R0, baseBodyMap, simulationVariables);

    // Create EDT config class and set constant thrust in guidance class
    if (verbosity) {std::cout<< " -- Creating Config class -- " << std::endl;}
    nlohmann::json configVariables = simulationVariables["EDTConfigs"];
    nlohmann::json SRPVariables = simulationVariables["scConfigs"]["SRP"];
    nlohmann::json materialProperties = simulationVariables["materialProperties"];
    EDTs::EDTConfig CHBEDTConfig = EDTs::EDTConfig(simulationVariables, configVariables, configType, SRPVariables, materialProperties);


    // Create EDT Guidance class
    if (verbosity) {std::cout<< " -- Creating Guidance class -- " << std::endl;}
    std::string thrustMagnitudeConfig = simulationVariables["GuidanceConfigs"]["thrustMagnitudeConfig"];
    std::string thrustDirectionConfig = simulationVariables["GuidanceConfigs"]["thrustDirectionConfig"];


    EDTGuidance CHBEDTGuidance = EDTGuidance(
            thrustMagnitudeConfig,
            thrustDirectionConfig,
            baseBodyMap,
            CHBEDTEnviro,
            CHBEDTConfig,
            simulationVariables);
    CHBEDTGuidance.setThrustMagnitudeConstant(CHBEDTConfig.getConstantThrust()); // TODO: Check if this works / is needed


    // Get universal class for propagation bodies
    if (verbosity) {std::cout<< " -- Creating Propbodies class -- " << std::endl;}
//    nlohmann::json jsonBodiesToInclude = simulationVariables["Spice"]["bodiesToInclude"];

    univ::propBodies SSOPropBodies = univ::propBodies(CHBEDTConfig, CHBEDTGuidance, baseBodyMap, simulationVariables);

    // Get universal class for propagation settings + set vehicle initial state
    if (verbosity) {std::cout<< " -- Creating Propsettings class -- " << std::endl;}


    // Set vehicle state either directly from cartesian elements, or via keplerian elements, depending which stateType is specified
    std::string initialStateType = simulationVariables["GuidanceConfigs"]["initialStateType"];
    Eigen::Vector6d vehicleInitialCartesian;

    // Directly set cartesian coordinates
    if (initialStateType == "Cartesian"){
        vehicleInitialCartesian[0] = simulationVariables["GuidanceConfigs"]["vehicleInitialCartesian"]["x1_m"];
        vehicleInitialCartesian[1] = simulationVariables["GuidanceConfigs"]["vehicleInitialCartesian"]["x2_m"];
        vehicleInitialCartesian[2] = simulationVariables["GuidanceConfigs"]["vehicleInitialCartesian"]["x3_m"];
        vehicleInitialCartesian[3] = simulationVariables["GuidanceConfigs"]["vehicleInitialCartesian"]["v1_ms"];
        vehicleInitialCartesian[4] = simulationVariables["GuidanceConfigs"]["vehicleInitialCartesian"]["v2_ms"];
        vehicleInitialCartesian[5] = simulationVariables["GuidanceConfigs"]["vehicleInitialCartesian"]["v3_ms"];
    }
    // Set cartesian coordinates via keplerian ones
    else if (initialStateType=="Keplerian"){
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
        vehicleInitialCartesian = convertKeplerianToCartesianElements(
                vehicleInitialKeplerian,
                solarGravPar);
    }
    else{
        std::cout<< "ERROR: intialStateType [" << initialStateType << "] not recognised!" << std::endl;
    }


//    std::cout << "Vehicle initial cartesian components: "
//              << vehicleInitialCartesian[0] << " "
//              << vehicleInitialCartesian[1] << " "
//              << vehicleInitialCartesian[2] << " "
//              << vehicleInitialCartesian[3] << " "
//              << vehicleInitialCartesian[4] << " "
//              << vehicleInitialCartesian[5]
//              << std::endl; //TODO:Remove me

    Eigen::Vector3d vehicleInitialPosition;
    vehicleInitialPosition << vehicleInitialCartesian[0], vehicleInitialCartesian[1], vehicleInitialCartesian[2];
    Eigen::Vector3d vehicleInitialVelocity;
    vehicleInitialVelocity << vehicleInitialCartesian[3], vehicleInitialCartesian[4], vehicleInitialCartesian[5];

    // Set various json variables
    double initialEphemerisYear = simulationVariables["GuidanceConfigs"]["initialEphemerisYear"];
    nlohmann::json  terminationSettingsJson = simulationVariables["GuidanceConfigs"]["terminationSettings"];
    nlohmann::json integratorSettingsJson = simulationVariables["GuidanceConfigs"]["integratorSettings"];
//    std::cout << ""

    // Actually create the propSettings class
    univ::propSettings SSOPropSettings = univ::propSettings(SSOPropBodies,
                                                            vehicleInitialPosition,
                                                            vehicleInitialVelocity,
                                                            integratorSettingsJson,
                                                            terminationSettingsJson,
                                                            simulationVariables,
                                                            initialEphemerisYear);



    // Ensure environment is properly updated
    CHBEDTGuidance.updateAllEnviro();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (verbosity) {std::cout<< "====================Propagating==========================" << std::endl;}

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
            SSOPropBodies.bodyMap, SSOPropSettings.integratorSettings, SSOPropSettings.propagatorSettings);
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////             SAVE DATA                  ////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (verbosity) {std::cout<< "====================Saving Data==========================" << std::endl;}

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
    bool saveCurrentVNVData = dataTypesToSave["currentVNV"];
    bool saveConfigInfo = dataTypesToSave["configInfo"];

//    std::cout << integrationResult[1294891093.79353] << std::endl;
////    std::cout << integrationResult.find(1294891093.79353)->second << std::endl;
////    // Remove unwanted rows from the custom save maps
//    std::map<double, Eigen::VectorXd> newMap;
//    for (std::map<double, Eigen::VectorXd>::iterator i=integrationResult.begin(); i!=integrationResult.end(); i++ ){
//        double time = i->first;
//        newMap[time] = CHBEDTGuidance.getMagFieldMap()[time];
//    }

    // Check which data types to save, and save them
    if (savePropData) {
        if (verbosity) {std::cout << "Saving propData" << std::endl;}

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
        if (verbosity) {std::cout << "Ignoring propData" << std::endl;}
    }

    if (saveMagneticFieldData) {
        if (verbosity) {std::cout << "Saving magneticFieldData" << std::endl;}

        // Make custom map conform to regular timestep saving
        std::map<double, Eigen::VectorXd> magFieldMap = CHBEDTGuidance.getMagFieldMap();
        std::map<double, Eigen::VectorXd> newMagFieldMap;

        for (std::map<double, Eigen::VectorXd>::iterator i=integrationResult.begin(); i!=integrationResult.end(); i++ ){
            double time = i->first;
            newMagFieldMap[time] = magFieldMap[time];
        }


        // Write magField history to file.
        input_output::writeDataMapToTextFile(newMagFieldMap,
                                             baseFilename + "magData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        if (verbosity) {std::cout << "Ignoring magneticFieldData" << std::endl;}
    }

    if (saveIonosphereData) {
        if (verbosity) {std::cout << "Saving ionosphericData" << std::endl;}

        // Make custom map conform to regular timestep saving
        std::map<double, Eigen::VectorXd> ionosphereMap = CHBEDTGuidance.getIonosphereMap();
        std::map<double, Eigen::VectorXd> newIonosphereMap;
        for (std::map<double, Eigen::VectorXd>::iterator i=integrationResult.begin(); i!=integrationResult.end(); i++ ){
            double time = i->first;
            newIonosphereMap[time] = ionosphereMap[time];
        }

        // Write ionosphere history to file.
        input_output::writeDataMapToTextFile(newIonosphereMap,
                                             baseFilename + "ionoData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        if (verbosity) {std::cout << "Ignoring ionosphericData" << std::endl;}
    }

    if (saveThrustData) {
        if (verbosity) {std::cout << "Saving thrustData" << std::endl;}

        // Make custom map conform to regular timestep saving
        std::map<double, Eigen::VectorXd> thrustMap = CHBEDTGuidance.getThrustMap();
        std::map<double, Eigen::VectorXd> newThrustMap;
        for (std::map<double, Eigen::VectorXd>::iterator i=integrationResult.begin(); i!=integrationResult.end(); i++ ){
            double time = i->first;
            newThrustMap[time] = thrustMap[time];
        }

        // Write thrust history to file.
        input_output::writeDataMapToTextFile(newThrustMap,
                                             baseFilename + "thrustData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        if (verbosity) {std::cout << "Ignoring thrustData" << std::endl;}
    }

    if (saveCurrentData) {
        if (verbosity) {std::cout << "Saving currentData" << std::endl;}

        // Make custom map conform to regular timestep saving
        std::map<double, Eigen::VectorXd> currentMap = CHBEDTGuidance.getCurrentMap();
        std::map<double, Eigen::VectorXd> newCurrentMap;
        for (std::map<double, Eigen::VectorXd>::iterator i=integrationResult.begin(); i!=integrationResult.end(); i++ ){
            double time = i->first;
            newCurrentMap[time] = currentMap[time];
        }

        // Write current history to file.
        input_output::writeDataMapToTextFile(newCurrentMap,
                                             baseFilename + "currentData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        if (verbosity) {std::cout << "Ignoring currentData" << std::endl;}
    }

    if (saveCurrentVNVData) {
        if (verbosity) {std::cout << "Saving currentVNVData" << std::endl;}

        // Make custom map conform to regular timestep saving
        std::map<double, Eigen::VectorXd> currentVNVMap = CHBEDTGuidance.getCurrentVNVMap();
        std::map<double, Eigen::VectorXd> newCurrentVNVMap;
        for (std::map<double, Eigen::VectorXd>::iterator i=integrationResult.begin(); i!=integrationResult.end(); i++ ){
            double time = i->first;
            newCurrentVNVMap[time] = currentVNVMap[time];
        }

        // Write dependent variables to file
        input_output::writeDataMapToTextFile(newCurrentVNVMap,
                                             baseFilename + "currentVNVData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        if (verbosity) {std::cout << "Ignoring currentVNVData" << std::endl;}
    }

    if (saveBodyData) {
        if (verbosity) {std::cout << "Saving bodyData" << std::endl;}

        // Make custom map conform to regular timestep saving
        std::map<double, Eigen::VectorXd> bodyDataMap = CHBEDTGuidance.getBodyDataMap();
        std::map<double, Eigen::VectorXd> newBodyDataMap;
        for (std::map<double, Eigen::VectorXd>::iterator i=integrationResult.begin(); i!=integrationResult.end(); i++ ){
            double time = i->first;
            newBodyDataMap[time] = bodyDataMap[time];
        }

        // Write bodyData history to files
        input_output::writeDataMapToTextFile(newBodyDataMap,
                                             baseFilename + "bodyData.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        if (verbosity) {std::cout << "Ignoring bodyData" << std::endl;}
    }

    if (saveDependentVariablesData) {
        if (verbosity) {std::cout << "Saving dependentVariableData" << std::endl;}

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
        if (verbosity) {std::cout << "Ignoring dependentVariableData" << std::endl;}
    }

    ///// This portion saves all relevant info about the configuration including derived variables such as mass etc ////

    // Initialise save map and vector
    std::map < std::string, Eigen::VectorXd > configInfoMap;
    Eigen::VectorXd configInfoVector;
    configInfoVector.resize(31);
//    std::string TEMPNAME = "Boop";
    configInfoVector[0] = CHBEDTConfig.getResistivity();
    configInfoVector[1] = CHBEDTConfig.getVehicleConductivity();
    configInfoVector[2] = CHBEDTConfig.getEDTResistance();
    configInfoVector[3] = CHBEDTConfig.getVehicleMass();
    configInfoVector[4] = CHBEDTConfig.getTetherMass();
    configInfoVector[5] = CHBEDTConfig.getCompositeMaterialDensity();
    configInfoVector[6] = CHBEDTConfig.getTetherAreaInner();
    configInfoVector[7] = CHBEDTConfig.getTetherAreaOuter();
    configInfoVector[8] = CHBEDTConfig.getTetherAreaInnerSecondary();
    configInfoVector[9] = CHBEDTConfig.getTetherAreaOuterSecondary();
    configInfoVector[10] = CHBEDTConfig.getTetherDiameterInner();
    configInfoVector[11] = CHBEDTConfig.getTetherDiameterOuter();
    configInfoVector[12] = CHBEDTConfig.getTetherDiameterInnerSecondary();
    configInfoVector[13] = CHBEDTConfig.getTetherDiameterOuterSecondary();
    configInfoVector[14] = CHBEDTConfig.getXSecAreaTotal();
    configInfoVector[15] = CHBEDTConfig.getXSecAreaInner();
    configInfoVector[16] = CHBEDTConfig.getXSecAreaOuter();
    configInfoVector[17] = CHBEDTConfig.getXSecAreaConducting();
    configInfoVector[18] = CHBEDTConfig.getNoTetherSegments();
    configInfoVector[19] = CHBEDTConfig.getNoPrimaryLinks();
    configInfoVector[20] = CHBEDTConfig.getNoSecondaryLinks();
    configInfoVector[21] = CHBEDTConfig.getTotalPrimaryLineLength();
    configInfoVector[22] = CHBEDTConfig.getTotalSecondaryLineLength();
    configInfoVector[23] = CHBEDTConfig.getPrimaryLineSeparation();
    configInfoVector[24] = CHBEDTConfig.getPrimaryLineSegmentLength();
    configInfoVector[25] = CHBEDTConfig.getSecondaryLineSegmentLength();
    configInfoVector[26] = CHBEDTConfig.getTotalPrimarySRPArea();
    configInfoVector[27] = CHBEDTConfig.getTotalSecondarySRPArea();
    configInfoVector[28] = CHBEDTConfig.getEffectiveHoytetherSRPArea();
    configInfoVector[29] = CHBEDTConfig.getEffectiveSRPArea();
    configInfoVector[30] = CHBEDTConfig.getEffectiveRadiationPressureCoefficient();



    configInfoMap.insert(std::pair<std::string, Eigen::VectorXd> (jsonName, configInfoVector));

    if (saveConfigInfo) {
        if (verbosity) {std::cout << "Saving configInfo" << std::endl;}

        // Write dependent variables to file
        input_output::writeDataMapToTextFile(configInfoMap,
                                             baseFilename + "configInfo.dat",
                                             outputPath,
                                             "",
                                             std::numeric_limits<double>::digits10,
                                             std::numeric_limits<double>::digits10,
                                             ",");
    }
    else{
        if (verbosity) {std::cout << "Ignoring configInfo" << std::endl;}
    }



//    std::cout << "Program ran for " << integrationResult[0] << std::endl;


    std::cout << "Task failed successfully" << std::endl;
}