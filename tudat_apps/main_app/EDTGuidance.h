//
// Created by matt on 23/03/2020.
//

#ifndef TUDATBUNDLE_EDTGUIDANCE_H
#define TUDATBUNDLE_EDTGUIDANCE_H

#include "EDT_configs.h"
#include "environment_settings.h"

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace EDTs;

class EDTGuidance{
public:
    // Thrust direction and magnitude settings pointers
    std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings;
    std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings;

    EDTGuidance(std::string& thrustMagnitudeConfig,
            std::string& thrustDirectionConfig,
            NamedBodyMap& environmentBodyMap,
            EDTEnvironment& guidanceEnvironment,
            EDTs::EDTConfig& vehicleConfig,
            nlohmann::json simulationVariables,
            const double throttleSetting = 1 ):
            thrustMagnitudeConfig_(thrustMagnitudeConfig),
            thrustDirectionConfig_(thrustDirectionConfig),
            environmentBodyMap_(environmentBodyMap),
            guidanceEnvironment_(guidanceEnvironment),
            vehicleConfig_(vehicleConfig),
            simulationVariables_(simulationVariables),
            throttleSetting_(throttleSetting){

        // Set Sun gravitation parameter
        SunMu_ = simulationVariables_["constants"]["gravitationalParameters"]["Sun"];

        // Set sizes for save vectors
        magFieldSaveVector_.resize(9);
        ionosphereSaveVector_.resize(1);
        thrustSaveVector_.resize(7);
        currentSaveVector_.resize(7);
        currentVNVSaveVector_.resize(17);
        bodyDataSaveVector_.resize(7);

        // Set configType_ from config class
        configType_ = vehicleConfig_.getConfigType();

        // Do initialisation update of guidance settings
        std::cout << "Updating guidance settings" << std::endl;
        updateGuidanceSettings();

        // Perihelion limiter
        minimumPerihelionAU_ = simulationVariables_["GuidanceConfigs"]["minPeAU"];
        minimumPerihelion_ = minimumPerihelionAU_ * AU;
    };

    /////////////////////// Update functions based on config settings for custom thrust ///////////////////////////////
    // TODO: Consider updating environment differently (currently done in both magnitude and direction functions)
    // Function to update thrust magnitude based on time and the config used
    void updateThrustMagnitude(const double simulationTime){
        updateAllEnviro(simulationTime_);
        if (forceZeroCurrentMagnitude_ or (thrustDirectionConfig_ == "disabled")){
            thrustMagnitude_ = 0;
        }
        else{
            if (thrustMagnitudeConfig_ == "constant"){
                thrustMagnitude_ = thrustMagnitudeConstant_;
            }
            else if (thrustMagnitudeConfig_ == "nominal") {

                thrustVector_ = (guidanceEnvironment_.getCurrent()).cross(guidanceEnvironment_.getMagFieldInertial());
                thrustMagnitude_ = thrustVector_.norm();
            }
        }
//        std::cout << "Thrust magnitude: " << thrustMagnitude_ <<  std::endl; //TODO: remove me
//        std::cout << "Thrust vector: " << thrustVector_ <<  std::endl;\\TODO: remove me
//        std::cout << "Thrust vector local: " << thrustVectorLocal_ <<  std::endl;\\TODO: remove me
        // Save simulation time and relevant data to save map
        simulationTime_ = simulationTime;
        saveParametersToMap();
    }

    void updateThrustDirection(const double simulationTime){

//        ///////// TODO: Below is old version of thrust direction logic, use the new one!!! /////////////////////
//        // Set current vector for prograde and retrograde, normalised vector pointing upwards and downwards
//        if ((thrustDirectionConfig_ == "nominalPrograde") or (thrustDirectionConfig_ == "nominalRetrograde") or (thrustDirectionConfig_ == "currentArbitraryVNV")){
//            if (thrustDirectionConfig_ == "nominalPrograde") {
//
//                // For prograde use current oriented in positive Z-direction
//                currentDirectionVector_ << 0, 0, 1; //TODO: Fix to use EDT length unit vector, instead of hardcode
//            }
//            else if (thrustDirectionConfig_ == "nominalRetrograde") {
//
//                // For retrograde use current oriented in negative Z-direction
//                currentDirectionVector_ << 0, 0, -1; //TODO: Fix to use EDT length unit vector, instead of hardcode
//            }
//            else if (thrustDirectionConfig_ == "currentArbitraryVNV"){
//
//                // Manually use the value used for the arbitrary VNV case
//                currentDirectionVector_ << 0.17583754,  0.87544648, -0.45019399;
//            }
////            std::cout << "CurrentDirectionVector at set: " << currentDirectionVector_ << std::endl;
//            // Set current in guidance class and update it
////            guidanceEnvironment_.setCurrent(currentDirectionVector_); // TODO: check if this is needed or no - ie move to update script
//            updateAllEnviro(simulationTime_);
//
//            // Get thrust vector and normalise to use as guidance vector
//            thrustVector_ = (guidanceEnvironment_.getCurrent()).cross(guidanceEnvironment_.getMagFieldInertial());
//            thrustDirection_ = thrustVector_.normalized();
//
//            // Convert inertial thrust vector to local vector for saving
//            thrustVectorLocal_ = gen::InertialToLvlh(thrustVector_, guidanceEnvironment_.getTheta());
//        }
//
//        // Save simulation time and relevant data to save map
//        simulationTime_ = simulationTime;
//        saveParametersToMap();



        ////////////// THE NEW ONE ///////////////////////////////

        // Define positive and negative current direction vectors
        Eigen::Vector3d currentPositiveDirectionVector;
        currentPositiveDirectionVector << 0, 0, 1;
        Eigen::Vector3d currentNegativeDirectionVector;
        currentNegativeDirectionVector << 0, 0, -1;
        Eigen::Vector3d currentDirectionVectorNoThrust;
        currentDirectionVectorNoThrust << 0, 0, -1;
        Eigen::Vector3d currentDirectionArbitraryVnV;
        currentDirectionArbitraryVnV <<  0.17583754,  0.87544648, -0.45019399;

        // Initialise boolean forcing current magnitude to 0, if thrusting should be disabled
        forceZeroCurrentMagnitude_ = false;

        // Grab spacecraft state at the current time and convert to keplerian, then calculate Ap and Pe for T0
        Eigen::Vector6d spacecraftStateCartesianT0 = guidanceEnvironment_.getVehicleState();
        Eigen::Vector6d spacecraftStateKeplerianT0 = tudat::orbital_element_conversions::convertCartesianToKeplerianElements(spacecraftStateCartesianT0, SunMu_);
        double ApT0 =  spacecraftStateKeplerianT0[0] * (1 + spacecraftStateKeplerianT0[1]);
        double PeT0 = spacecraftStateKeplerianT0[0] * (1 - spacecraftStateKeplerianT0[1]);

        // Grab thrust direction at T0, for both positive and negative current cases
        Eigen::Vector3d thrustVectorPositive = currentPositiveDirectionVector.cross(guidanceEnvironment_.getMagFieldInertial());
        Eigen::Vector3d thrustDirectionPositive = thrustVectorPositive.normalized();

        Eigen::Vector3d thrustVectorNegative = currentNegativeDirectionVector.cross(guidanceEnvironment_.getMagFieldInertial());
        Eigen::Vector3d thrustDirectionNegative = thrustVectorNegative.normalized();

        // Use normalised thrust direction directly as a 1m/s impulse on the original state velocity, to get new Ap and Pe, for both positive and negative cases
        Eigen::Vector6d spacecraftStateDeltaPositive;
        spacecraftStateDeltaPositive << 0, 0, 0, thrustDirectionPositive[0], thrustDirectionPositive[1], thrustDirectionPositive[2];
        Eigen::Vector6d spacecraftStateCartesianT1Positive = spacecraftStateCartesianT0 + spacecraftStateDeltaPositive;
        Eigen::Vector6d spacecraftStateKeplerianT1Positive = tudat::orbital_element_conversions::convertCartesianToKeplerianElements(spacecraftStateCartesianT1Positive, SunMu_);
        double ApT1Positive = spacecraftStateKeplerianT1Positive[0] * (1 + spacecraftStateKeplerianT1Positive[1]);
        double PeT1Positive = spacecraftStateKeplerianT1Positive[0] * (1 - spacecraftStateKeplerianT1Positive[1]);

        Eigen::Vector6d spacecraftStateDeltaNegative;
        spacecraftStateDeltaNegative << 0, 0, 0, thrustDirectionNegative[0], thrustDirectionNegative[1], thrustDirectionNegative[2];
        Eigen::Vector6d spacecraftStateCartesianT1Negative = spacecraftStateCartesianT0 + spacecraftStateDeltaNegative;
        Eigen::Vector6d spacecraftStateKeplerianT1Negative = tudat::orbital_element_conversions::convertCartesianToKeplerianElements(spacecraftStateCartesianT1Negative, SunMu_);
        double ApT1Negative = spacecraftStateKeplerianT1Negative[0] * (1 + spacecraftStateKeplerianT1Negative[1]);
        double PeT1Negative = spacecraftStateKeplerianT1Negative[0] * (1 - spacecraftStateKeplerianT1Negative[1]);

        // Thrust logic only used at lower altitudes, ie when altitude less than 10 AU
        bool useThrustLogic;
        double altitude = guidanceEnvironment_.getVehiclePosition().norm();
        if (altitude > 10*AU){
            useThrustLogic = false;
        }
        else{
            useThrustLogic = true;
        }

        // For whether using nominalPrograde or nominalRetrograde, apply the relevant logic to achieve a final thrust direction
        if ((thrustDirectionConfig_ == "nominalPrograde") and (useThrustLogic) ){

            if ((ApT1Positive > ApT0) and (PeT1Positive > PeT0)){
                currentDirectionVector_ = currentPositiveDirectionVector;
            }
            else if ((ApT1Negative > ApT0) and (PeT1Negative > PeT0)){
                currentDirectionVector_ = currentNegativeDirectionVector;
            }
            else if ((ApT1Positive > ApT0) and (PeT1Positive < PeT0) and (PeT1Positive > minimumPerihelion_)){
                currentDirectionVector_ = currentPositiveDirectionVector;
            }
            else if ((ApT1Negative > ApT0) and (PeT1Negative < PeT0) and (PeT1Negative > minimumPerihelion_)){
                currentDirectionVector_ = currentNegativeDirectionVector;
            }
            else{
                currentDirectionVector_ = currentDirectionVectorNoThrust;
                forceZeroCurrentMagnitude_ = true;
            }
        }

        else if ((thrustDirectionConfig_ == "nominalRetrograde") and (useThrustLogic)){

            if ((ApT1Positive > ApT0) and (PeT1Positive < PeT0) and (PeT1Positive > minimumPerihelion_)){
                currentDirectionVector_ = currentPositiveDirectionVector;
            }
            else if ((ApT1Negative > ApT0) and (PeT1Negative < PeT0) and (PeT1Negative > minimumPerihelion_)){
                currentDirectionVector_ = currentNegativeDirectionVector;
            }
            else if ((ApT1Positive < ApT0) and (PeT1Positive < PeT0) and (PeT1Positive > minimumPerihelion_)){
                currentDirectionVector_ = currentPositiveDirectionVector;
            }
            else if ((ApT1Negative < ApT0) and (PeT1Negative < PeT0) and (PeT1Negative > minimumPerihelion_)){
                currentDirectionVector_ = currentNegativeDirectionVector;
            }
            else{
                currentDirectionVector_ = currentDirectionVectorNoThrust;
                forceZeroCurrentMagnitude_ = true;
            }
        }

        else if (thrustDirectionConfig_ == "currentArbitraryVNV"){
            currentDirectionVector_ = currentDirectionArbitraryVnV;
        }

//        std::cout << "Force current bool: " << forceZeroCurrentMagnitude_ << std::endl;
        updateAllEnviro(simulationTime_, forceZeroCurrentMagnitude_);

//        std::cout << "Current is: " << guidanceEnvironment_.getCurrent() << std::endl; // TODO: remove me
        // Get thrust vector and normalise to use as guidance vector
        thrustVector_ = (guidanceEnvironment_.getCurrent()).cross(guidanceEnvironment_.getMagFieldInertial());
        thrustDirection_ = thrustVector_.normalized();

        // Convert inertial thrust vector to local vector for saving
        thrustVectorLocal_ = gen::InertialToLvlh(thrustVector_, guidanceEnvironment_.getTheta());

        // Save simulation time and relevant data to save map
        simulationTime_ = simulationTime;
        saveParametersToMap();

    }

    void updateGuidanceSettings(){
        // Bind functions properly
        thrustMagnitudeFunction_ = std::bind( &EDTGuidance::getThrustMagnitude, this, std::placeholders::_1 );
        thrustDirectionFunction_ = std::bind( &EDTGuidance::getThrustDirection, this, std::placeholders::_1);

        // Set thrust direction and magnitude using custom settings
        thrustMagnitudeSettings =
                std::make_shared< FromFunctionThrustMagnitudeSettings >(
                        thrustMagnitudeFunction_,
                        [ = ]( const double ){ return 999; });

        thrustDirectionGuidanceSettings =
                std::make_shared<CustomThrustDirectionSettings>(
                        thrustDirectionFunction_);
    }

    void saveParametersToMap(){
        // Save magfield parameters
        magFieldSaveVector_[0] = guidanceEnvironment_.getMagFieldMagnitude();
        magFieldSaveVector_[1] = guidanceEnvironment_.getTheta();
        magFieldSaveVector_[2] = guidanceEnvironment_.getMagFieldLocal()[0];
        magFieldSaveVector_[3] = guidanceEnvironment_.getMagFieldLocal()[1];
        magFieldSaveVector_[4] = guidanceEnvironment_.getMagFieldLocal()[2];
        magFieldSaveVector_[5] = guidanceEnvironment_.getMagFieldInertial()[0];
        magFieldSaveVector_[6] = guidanceEnvironment_.getMagFieldInertial()[1];
        magFieldSaveVector_[7] = guidanceEnvironment_.getMagFieldInertial()[2];
        magFieldSaveVector_[8] = guidanceEnvironment_.getB0();

        magFieldMap_.insert(std::pair<double, Eigen::VectorXd> (simulationTime_, magFieldSaveVector_));

        // Save Ionosphere parameters TODO: TBF
        ionosphereSaveVector_[0] = 12.3;
        ionosphereMap_.insert(std::pair<double, Eigen::VectorXd> (simulationTime_, ionosphereSaveVector_));

        // Save Thrust parameters
        thrustSaveVector_[0] = thrustMagnitude_;
        thrustSaveVector_[1] = thrustVector_[0];
        thrustSaveVector_[2] = thrustVector_[1];
        thrustSaveVector_[3] = thrustVector_[2];
        thrustSaveVector_[4] = thrustVectorLocal_[0];
        thrustSaveVector_[5] = thrustVectorLocal_[1];
        thrustSaveVector_[6] = thrustVectorLocal_[2];
        thrustMap_.insert(std::pair<double, Eigen::VectorXd> (simulationTime_, thrustSaveVector_));


        // Save Current parameters TODO: TBF
        currentSaveVector_[0] = guidanceEnvironment_.getCurrentMagnitude();
        currentSaveVector_[1] = guidanceEnvironment_.getCurrent()[0];
        currentSaveVector_[2] = guidanceEnvironment_.getCurrent()[1];
        currentSaveVector_[3] = guidanceEnvironment_.getCurrent()[2];
        currentSaveVector_[4] = velocityWrtMagField_[0];
        currentSaveVector_[5] = velocityWrtMagField_[1];
        currentSaveVector_[6] = velocityWrtMagField_[2];
        currentMap_.insert(std::pair<double, Eigen::VectorXd> (simulationTime_, currentSaveVector_));

        // Save current intermediary VNV data
        currentVNVSaveVector_[0] = motionalEMF_;
        currentVNVSaveVector_[1] = unitCurrent_;
        currentVNVSaveVector_[2] = dimensionlessCurrentC_;
        currentVNVSaveVector_[3] = dimensionlessVoltageA_;
        currentVNVSaveVector_[4] = avgDimensionlessCurrent_;
        currentVNVSaveVector_[5] = avgTrueCurrent_;
        currentVNVSaveVector_[6] = guidanceEnvironment_.getCurrent()[0];
        currentVNVSaveVector_[7] = guidanceEnvironment_.getCurrent()[1];
        currentVNVSaveVector_[8] = guidanceEnvironment_.getCurrent()[2];
        currentVNVSaveVector_[9] = trueCurrentC_;

        currentVNVSaveVector_[10] = guidanceEnvironment_.getCurrentDirectionVector()[0];
        currentVNVSaveVector_[11] = guidanceEnvironment_.getCurrentDirectionVector()[1];
        currentVNVSaveVector_[12] = guidanceEnvironment_.getCurrentDirectionVector()[2];
        currentVNVSaveVector_[13] = vehicleConfig_.getVehicleConductivity();
        currentVNVSaveVector_[14] = vehicleConfig_.getXSecAreaConducting();
        currentVNVSaveVector_[15] = vehicleConfig_.getTetherLength();
        currentVNVSaveVector_[16] = forceZeroCurrentMagnitude_;
        currentVNVMap_.insert(std::pair<double, Eigen::VectorXd> (simulationTime_, currentVNVSaveVector_));
//        std::cout << "CURRENTVNV SAVE map;;;;;;; " << currentVNVSaveMap_ << std::endl; // TODO: remove me

        // Save body data
        bodyDataSaveVector_[0] = guidanceEnvironment_.getR();
        bodyDataSaveVector_[1] = guidanceEnvironment_.getVehicleState()[0];
        bodyDataSaveVector_[2] = guidanceEnvironment_.getVehicleState()[1];
        bodyDataSaveVector_[3] = guidanceEnvironment_.getVehicleState()[2];
        bodyDataSaveVector_[4] = guidanceEnvironment_.getVehicleState()[3];
        bodyDataSaveVector_[5] = guidanceEnvironment_.getVehicleState()[4];
        bodyDataSaveVector_[6] = guidanceEnvironment_.getVehicleState()[5];
//        bodyDataSaveVector_[7] =
        bodyDataMap_.insert(std::pair<double, Eigen::VectorXd> (simulationTime_, bodyDataSaveVector_));

    }

    //////////////////////////////////// Update all function moved from environment_settings ///////////////

    // Function to update current, different depending on EDT config
    void updateCurrentMagnitude(){
//        std::cout << "Config type: " << configType_ <<  std::endl; //TODO: Remove me
        // For spacecraft using a bare tether concept
        if ( (configType_ == "CHB") or (configType_ == "AlMB") ){

            // Set intermediate values, from simulation data (or calculated elsewhere)
            velocityWrtMagField_ = guidanceEnvironment_.getVehicleVelocity(); // TODO: This just uses inertial velocity as relative velocity, make sure solar magfield is not moving (eg with solar wind)
            trueCurrentC_ = vehicleConfig_.getEmitterCurrent(); // TODO: Check emitter current can be the true current - ie check against ionosphere assumptions for intake current

//            std::cout<< "velocityWrtMagField_ : " << velocityWrtMagField_ << std::endl; // TODO: Remove me
//            std::cout<< "trueCurrentC_ : " << trueCurrentC_ << std::endl; // TODO: Remove me
//
//            std::cout<< "magfield: " << guidanceEnvironment_.getMagFieldInertial() << std::endl; // TODO: Remove me
//            std::cout << "currentDirectionVecotr: " << guidanceEnvironment_.getCurrentDirectionVector() << std::endl; // TODO: remove me
//            std::cout << "sigma: " << vehicleConfig_.getVehicleConductivity() << std::endl; // TODO: remove me
//            std::cout << "Area: " << vehicleConfig_.getXSecAreaConducting() << std::endl; // TODO: remove me

            // Calculate unit current, using sigma, Em and A
            motionalEMF_ = gen::getMotionalEMF(velocityWrtMagField_,
                                               guidanceEnvironment_.getMagFieldInertial(),
                                               guidanceEnvironment_.getCurrentDirectionVector()); // Note: EDT unit vector is the same as the current unit vector, therefore can use it directly for the final variable

//            std::cout << "vehicle state: "
//                    << guidanceEnvironment_.getVehicleState()[0] << " "
//                    << guidanceEnvironment_.getVehicleState()[1] << " "
//                    << guidanceEnvironment_.getVehicleState()[2] << " "
//                    << guidanceEnvironment_.getVehicleState()[3] << " "
//            << guidanceEnvironment_.getVehicleState()[4] << " "
//            << guidanceEnvironment_.getVehicleState()[5]
//            << std::endl; //TODO:Remove me
//            std::cout << "velocity wrt magfield components: " << velocityWrtMagField_[0] << " " << velocityWrtMagField_[1] << " " << velocityWrtMagField_[2] << std::endl; //TODO:Remove me
//            std::cout << "magfield inertial components: " << guidanceEnvironment_.getMagFieldInertial()[0] << " " << guidanceEnvironment_.getMagFieldInertial()[1] << " " << guidanceEnvironment_.getMagFieldInertial()[2] << std::endl; //TODO:Remove me
//            std::cout << "Current direction vector: " << guidanceEnvironment_.getCurrentDirectionVector()[0] << " " << guidanceEnvironment_.getCurrentDirectionVector()[1] << " " << guidanceEnvironment_.getCurrentDirectionVector()[2] << std::endl; //TODO:Remove me
//            std::cout<< "motionalEMF_ : " << motionalEMF_ << std::endl; // TODO: Remove me

            unitCurrent_ = gen::getUnitCurrent(vehicleConfig_.getVehicleConductivity(),
                                               motionalEMF_,
                                               vehicleConfig_.getXSecAreaConducting());

//            std::cout<< "unitCUrrent_ : " << unitCurrent_ << std::endl; // TODO: Remove me

            // Calculate iavg, from ic, lambda_A and L
            dimensionlessCurrentC_ = gen::trueCurrentDimensionlessConvert(trueCurrentC_, unitCurrent_);
            dimensionlessVoltageA_ = gen::getDimensionlessVoltageA(dimensionlessCurrentC_);
            avgDimensionlessCurrent_ = gen::getAvgDimensionlessCurrent(vehicleConfig_.getTetherLength(), dimensionlessVoltageA_);

//            std::cout<< "dimensionlessCurrentC_ : " << dimensionlessCurrentC_ << std::endl; // TODO: Remove me
//            std::cout<< "dimensionlessVoltageA_ : " << dimensionlessVoltageA_ << std::endl; // TODO: Remove me
//            std::cout<< "avgDimensionlessCurrent_ : " << avgDimensionlessCurrent_ << std::endl; // TODO: Remove me

            // Calculate final average true current, and set current magnitude to it
            avgTrueCurrent_ = gen::trueCurrentDimensionlessConvert(avgDimensionlessCurrent_, unitCurrent_, false);
            currentMagnitude_ = std::abs(avgTrueCurrent_);


//            std::cout<< "avgTrueCurrent_ : " << avgTrueCurrent_ << std::endl; // TODO: Remove me
//            std::cout<< "currentMagnitude_ : " << currentMagnitude_ << std::endl; // TODO: Remove me
//            std::cout << std::endl; // TODO: remove me

        }



            // For spacecraft using a transient-current concept
        else if ( (configType_ == "CHTr") or (configType_ == "AlMTr") ){
            // TODO: Add an update function for transient-current type spacecraft
            std::cout << "No calculation for transient-current made yet" << std::endl;
        }
    }

    // Run all updaters (mostly from environment settings)
    void updateAllEnviro(double simulationTime=0, bool forceZeroCurrentMagnitude=false){


        // Update body and magfield directly from guidanceEnvironment class
        guidanceEnvironment_.updateBodyParameters(simulationTime);
//        std::cout << "vehicle cartesian, first update: " << guidanceEnvironment_.getVehicleState() << std::endl; // todo:remove me
        guidanceEnvironment_.updateMagField();



        // Update current by updating and setting the magnitude and direction separately (note the direction was set by the thrust direction directly)
        if (forceZeroCurrentMagnitude){
            currentMagnitude_ = 0;
        }
        else{
            updateCurrentMagnitude();
        }

        guidanceEnvironment_.setCurrentMagnitude(currentMagnitude_);
        guidanceEnvironment_.setCurrentDirectionVector(currentDirectionVector_);
        guidanceEnvironment_.setCurrent(currentMagnitude_ * currentDirectionVector_);


    }


    /////////////////////// Get variable functions for custom thrust /////////////////////////////////////////////

    double getThrustMagnitude(const double simulationTime){
        updateThrustMagnitude(simulationTime);
        return thrustMagnitude_;
    }

    Eigen::Vector3d getThrustDirection(const double simulationTime){
        updateThrustDirection(simulationTime);
        return thrustDirection_;
    }

    /////////////////////// Other misc functions /////////////////////////////////////////////

    // get and set functions for configs
    std::string getThrustMagnitudeConfig(){
        return thrustMagnitudeConfig_;
    }

    std::string getThrustDirectionConfig(){
        return thrustDirectionConfig_;
    }

    void setThrustMagnitudeConfig(std::string newConfig){
        thrustMagnitudeConfig_ = newConfig;
    }

    void setThrustDirectionConfig(std::string newConfig){
        thrustDirectionConfig_ = newConfig;
    }

    // Get and set for thrust magnitude and throttle

    double getThrustMagnitudeConstant(){
        return thrustMagnitudeConstant_;
    }

    void setThrustMagnitudeConstant(double newConstantThrust){
        thrustMagnitudeConstant_ = newConstantThrust;
    }

    double getThrottleSetting(){
        return throttleSetting_;
    }

    void setThrottleSetting(double newThrottleSetting){
        throttleSetting_ = newThrottleSetting;
    }

    std::map < double, Eigen::VectorXd > getMagFieldMap(){
        return magFieldMap_;
    }

    std::map < double, Eigen::VectorXd > getIonosphereMap(){
        return ionosphereMap_;
    }

    std::map < double, Eigen::VectorXd > getThrustMap(){
        return thrustMap_;
    }

    std::map < double, Eigen::VectorXd > getCurrentMap(){
        return currentMap_;
    }

    std::map < double, Eigen::VectorXd > getCurrentVNVMap(){
        return currentVNVMap_;
    }

    std::map<double, Eigen::VectorXd> getBodyDataMap(){
        return bodyDataMap_;
    }




protected:
    /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////
    // EDT config variables
    std::string& thrustMagnitudeConfig_;
    std::string& thrustDirectionConfig_;

    // environment body map
    NamedBodyMap& environmentBodyMap_;

    // Class for EDT environment
    EDTEnvironment& guidanceEnvironment_;

    // Vehicle config class
    EDTs::EDTConfig& vehicleConfig_;

    /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////
    // Thrust config variables (can be given default values)
    double thrustMagnitudeConstant_;

    // Factor used for throttling in the nominal mode (from 0-1)
    double throttleSetting_;

    /////////////////////////// Other set parameters ///////////////////////////////
    // Thrust vector, magnitude and direction variables
    Eigen::Vector3d thrustVector_;
    Eigen::Vector3d thrustVectorLocal_;
    double thrustMagnitude_;
    Eigen::Vector3d thrustDirection_;
    double minimumPerihelionAU_;
    double minimumPerihelion_;

    // Current vector used for thrust direction guidance
    Eigen::Vector3d currentDirectionVector_;

    // Simulation time
    double simulationTime_;

    // Save maps and vectors to send to pair
    // NOTE: to change sizes of saveVectors, see class constructor
    std::map < double, Eigen::VectorXd > magFieldMap_;
    Eigen::VectorXd magFieldSaveVector_;

    std::map < double, Eigen::VectorXd > ionosphereMap_;
    Eigen::VectorXd ionosphereSaveVector_;

    std::map < double, Eigen::VectorXd > thrustMap_;
    Eigen::VectorXd thrustSaveVector_;

    std::map < double, Eigen::VectorXd > currentMap_;
    Eigen::VectorXd currentSaveVector_;

    std::map < double, Eigen::VectorXd > currentVNVMap_;
    Eigen::VectorXd currentVNVSaveVector_;

    std::map < double, Eigen::VectorXd > bodyDataMap_;
    Eigen::VectorXd bodyDataSaveVector_;

    // Custom Thrust standard functions
    std::function< double(const double) > thrustMagnitudeFunction_;
    std::function< Eigen::Vector3d(const double) > thrustDirectionFunction_;

    // Variables for general current calcs
    Eigen::Vector3d current_; // TODO: Decide if this vector needed, or can just use magnitude and direction vector
    double currentMagnitude_;
    std::string configType_; // Same as config type found in EDT_configs.h, sourced from jsons
    // Variables specific to bare tether concept
    double avgTrueCurrent_;
    double avgDimensionlessCurrent_;
    double dimensionlessCurrentC_;
    double trueCurrentC_;
    double unitCurrent_;
//    double EDTConductivity_;
    double motionalEMF_;
    Eigen::Vector3d velocityWrtMagField_; // TODO: Consider turning this into a common variable
    double dimensionlessVoltageA_;
//    double tetherLength_; // TODO: Make a common variable from jsons etc


//    // Current factor, TODO: model properly!
//    double currentFactor_ = 1*10E4; //

    nlohmann::json simulationVariables_;
    double SunMu_;

    bool forceZeroCurrentMagnitude_;

};



#endif //TUDATBUNDLE_EDTGUIDANCE_H
