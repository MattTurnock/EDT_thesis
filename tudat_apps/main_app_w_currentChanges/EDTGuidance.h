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
            const double throttleSetting = 1 ):
            thrustMagnitudeConfig_(thrustMagnitudeConfig),
            thrustDirectionConfig_(thrustDirectionConfig),
            environmentBodyMap_(environmentBodyMap),
            guidanceEnvironment_(guidanceEnvironment),
            vehicleConfig_(vehicleConfig),
            throttleSetting_(throttleSetting){
        // Set sizes for save vectors
        magFieldSaveVector_.resize(9);
        ionosphereSaveVector_.resize(1);
        thrustSaveVector_.resize(7);
        currentSaveVector_.resize(4);
        bodyDataSaveVector_.resize(7);

        // Set configType_ from config class
        configType_ = vehicleConfig_.getConfigType();

        // Do initialisation update of guidance settings
        std::cout << "Updating guidance settings" << std::endl;
        updateGuidanceSettings();
    };

    /////////////////////// Update functions based on config settings for custom thrust ///////////////////////////////
    // TODO: Consider updating environment differently (currently done in both magnitude and direction functions)
    // Function to update thrust magnitude based on time and the config used
    void updateThrustMagnitude(const double simulationTime){
        if (thrustMagnitudeConfig_ == "constant"){
            thrustMagnitude_ = thrustMagnitudeConstant_;
        }
        else if (thrustMagnitudeConfig_ == "nominal") {
            updateAllEnviro(simulationTime_);
            thrustVector_ = (guidanceEnvironment_.getCurrent()).cross(guidanceEnvironment_.getMagFieldInertial());
            thrustMagnitude_ = thrustVector_.norm();
        }

        // Save simulation time and relevant data to save map
        simulationTime_ = simulationTime;
        saveParametersToMap();
    }

    void updateThrustDirection(const double simulationTime){

        // Set current vector for prograde and retrograde, normalised vector pointing upwards and downwards
        if ((thrustDirectionConfig_ == "nominalPrograde") or (thrustDirectionConfig_ == "nominalRetrograde") ){
            if (thrustDirectionConfig_ == "nominalPrograde") {

                // For prograde use current oriented in positive Z-direction
                currentDirectionVector_ << 0, 0, 1; //TODO: Fix to use EDT length unit vector, instead of hardcode
            }
            else if (thrustDirectionConfig_ == "nominalRetrograde") {

                // For retrograde use current oriented in negative Z-direction
                currentDirectionVector_ << 0, 0, -1; //TODO: Fix to use EDT length unit vector, instead of hardcode
            }
            // Set current in guidance class and update it
//            guidanceEnvironment_.setCurrent(currentDirectionVector_); // TODO: check if this is needed or no - ie move to update script
            updateAllEnviro(simulationTime_);

            // Get thrust vector and normalise to use as guidance vector
            thrustVector_ = (guidanceEnvironment_.getCurrent()).cross(guidanceEnvironment_.getMagFieldInertial());
            thrustDirection_ = thrustVector_.normalized();

            // Convert inertial thrust vector to local vector for saving
            thrustVectorLocal_ = gen::InertialToLvlh(thrustVector_, guidanceEnvironment_.getTheta());
        }

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

//        std::cout << ":: current save vector components: " << currentSaveVector_[0]
//        << currentSaveVector_[1]
//        << currentSaveVector_[2]
//        << currentSaveVector_[3]
//        << currentSaveVector_[4]
//        << currentSaveVector_[5]
//        << currentSaveVector_[6] << std::endl;

        // Save body data
        bodyDataSaveVector_[0] = guidanceEnvironment_.getR();
        bodyDataSaveVector_[1] = guidanceEnvironment_.getVehicleState()[0];
        bodyDataSaveVector_[2] = guidanceEnvironment_.getVehicleState()[1];
        bodyDataSaveVector_[3] = guidanceEnvironment_.getVehicleState()[2];
        bodyDataSaveVector_[4] = guidanceEnvironment_.getVehicleState()[3];
        bodyDataSaveVector_[5] = guidanceEnvironment_.getVehicleState()[4];
        bodyDataSaveVector_[6] = guidanceEnvironment_.getVehicleState()[5];
        bodyDataMap_.insert(std::pair<double, Eigen::VectorXd> (simulationTime_, bodyDataSaveVector_));

        std::cout << "vehicle state: "
                  << guidanceEnvironment_.getVehicleState()[0] << " "
                  << guidanceEnvironment_.getVehicleState()[1] << " "
                  << guidanceEnvironment_.getVehicleState()[2] << " "
                  << guidanceEnvironment_.getVehicleState()[3] << " "
                  << guidanceEnvironment_.getVehicleState()[4] << " "
                  << guidanceEnvironment_.getVehicleState()[5]
                  << std::endl; //TODO:Remove me


    }

    //////////////////////////////////// Update all function moved from environment_settings ///////////////

    // Function to update current, different depending on EDT config
    void updateCurrentMagnitude(){
        std::cout << "Config type: " << configType_ <<  std::endl; //TODO: Remove me
        // For spacecraft using a bare tether concept
        if ( (configType_ == "CHB") or (configType_ == "AlMB") ){

            // Set intermediate values, from simulation data (or calculated elsewhere)
            velocityWrtMagField_ = guidanceEnvironment_.getVehicleVelocity(); // TODO: This just uses inertial velocity as relative velocity, make sure solar magfield is not moving (eg with solar wind)
            trueCurrentC_ = vehicleConfig_.getEmitterCurrent(); // TODO: Check emitter current can be the true current - ie check against ionosphere assumptions for intake current

            std::cout<< "trueCurrentC_ : " << trueCurrentC_ << std::endl; // TODO: Remove me

            // Calculate unit current, using sigma, Em and A
            motionalEMF_ = gen::getMotionalEMF(velocityWrtMagField_,
                                               guidanceEnvironment_.getMagFieldInertial(),
                                               guidanceEnvironment_.getCurrentDirectionVector()); // Note: EDT unit vector is the same as the current unit vector, therefore can use it directly for the final variable

            std::cout << "vehicle state: "
                    << guidanceEnvironment_.getVehicleState()[0] << " "
                    << guidanceEnvironment_.getVehicleState()[1] << " "
                    << guidanceEnvironment_.getVehicleState()[2] << " "
                    << guidanceEnvironment_.getVehicleState()[3] << " "
            << guidanceEnvironment_.getVehicleState()[4] << " "
            << guidanceEnvironment_.getVehicleState()[5]
            << std::endl; //TODO:Remove me
            std::cout << "velocity wrt magfield components: " << velocityWrtMagField_[0] << " " << velocityWrtMagField_[1] << " " << velocityWrtMagField_[2] << std::endl; //TODO:Remove me
            std::cout << "magfield inertial components: " << guidanceEnvironment_.getMagFieldInertial()[0] << " " << guidanceEnvironment_.getMagFieldInertial()[1] << " " << guidanceEnvironment_.getMagFieldInertial()[2] << std::endl; //TODO:Remove me
            std::cout << "Current direction vector: " << guidanceEnvironment_.getCurrentDirectionVector()[0] << " " << guidanceEnvironment_.getCurrentDirectionVector()[1] << " " << guidanceEnvironment_.getCurrentDirectionVector()[2] << std::endl; //TODO:Remove me
            std::cout<< "motionalEMF_ : " << motionalEMF_ << std::endl; // TODO: Remove me

            unitCurrent_ = gen::getUnitCurrent(vehicleConfig_.getVehicleConductivity(),
                                               motionalEMF_,
                                               vehicleConfig_.getXSecAreaConducting());

            std::cout<< "unitCUrrent_ : " << unitCurrent_ << std::endl; // TODO: Remove me

            // Calculate iavg, from ic, lambda_A and L
            dimensionlessCurrentC_ = gen::trueCurrentDimensionlessConvert(trueCurrentC_, unitCurrent_);
            dimensionlessVoltageA_ = gen::getDimensionlessVoltageA(dimensionlessCurrentC_);
            avgDimensionlessCurrent_ = gen::getAvgDimensionlessCurrent(tetherLength_, dimensionlessVoltageA_);

            // Calculate final average true current, and set current magnitude to it
            avgTrueCurrent_ = gen::trueCurrentDimensionlessConvert(avgDimensionlessCurrent_, unitCurrent_);
            currentMagnitude_ = avgTrueCurrent_;

        }

            // For spacecraft using a transient-current concept
        else if ( (configType_ == "CHTr") or (configType_ == "AlMTr") ){
            // TODO: Add an update function for transient-current type spacecraft
            std::cout << "No calculation for transient-current made yet" << std::endl;
        }
    }

    // Run all updaters (mostly from environment settings)
    void updateAllEnviro(double simulationTime=0){


        // Update body and magfield directly from guidanceEnvironment class
        guidanceEnvironment_.updateBodyParameters(simulationTime);
        std::cout << "vehicle cartesian, first update: " << guidanceEnvironment_.getVehicleState() << std::endl; // todo:remove me
        guidanceEnvironment_.updateMagField();



        // Update current by updating and setting the magnitude and direction separately (note the direction was set by the thrust direction directly)
        updateCurrentMagnitude();
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
    double tetherLength_; // TODO: Make a common variable from jsons etc


//    // Current factor, TODO: model properly!
//    double currentFactor_ = 1*10E4; //

};



#endif //TUDATBUNDLE_EDTGUIDANCE_H
