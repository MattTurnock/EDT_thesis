//
// Created by matt on 23/03/2020.
//

#ifndef TUDATBUNDLE_EDTGUIDANCE_H
#define TUDATBUNDLE_EDTGUIDANCE_H


#include "environment_settings.h"

using namespace tudat;
using namespace tudat::simulation_setup;

class EDTGuidance{
public:
    EDTGuidance(std::string& thrustMagnitudeConfig,
            std::string& thrustDirectionConfig,
            NamedBodyMap& environmentBodyMap,
            EDTEnvironment& guidanceEnvironment,
            const double throttleSetting = 1 ):
            thrustMagnitudeConfig_(thrustMagnitudeConfig),
            thrustDirectionConfig_(thrustDirectionConfig),
            environmentBodyMap_(environmentBodyMap),
            guidanceEnvironment_(guidanceEnvironment),
            throttleSetting_(throttleSetting){
        // Set sizes for save vectors
        magFieldSaveVector_.resize(8);
        ionosphereSaveVector_.resize(1);
        thrustSaveVector_.resize(7);
        currentSaveVector_.resize(4);
        bodyDataSaveVector_.resize(7);
    };

    /////////////////////// Update functions based on config settings for custom thrust ///////////////////////////////
    // TODO: Consider updating environment differently (currently done in both magnitude and direction functions)
    // Function to update thrust magnitude based on time and the config used
    void updateThrustMagnitude(const double simulationTime){
        if (thrustMagnitudeConfig_ == "constant"){
            thrustMagnitude_ = thrustMagnitudeConstant_;
        }
        else if (thrustMagnitudeConfig_ == "nominal") {
            guidanceEnvironment_.updateAll();
            thrustVector_ = (guidanceEnvironment_.getCurrent()).cross(guidanceEnvironment_.getMagFieldInertial());
            thrustMagnitude_ = thrustVector_.norm();
//            std::cout << std::to_string(thrustMagnitude_) << std::endl ; // TODO: remove me (for printing acceleration)
        }

        // Save simulation time and relevant data to save map
        simulationTime_ = simulationTime;
        saveParametersToMap();
    }

    void updateThrustDirection(const double simulationTime){

        // Set current vector for prograde and retrograde, normalised vector pointing upwards and downwards
        if ((thrustDirectionConfig_ == "nominalPrograde") or (thrustDirectionConfig_ == "nominalRetrograde") ){
            if (thrustDirectionConfig_ == "nominalPrograde") {

                currentToSet_ << 0, 0, 1*currentFactor_; //TODO: implement properly
            }
            else if (thrustDirectionConfig_ == "nominalRetrograde") {
                // Set current vector for retrograde, normalised vector pointing downwards
                currentToSet_ << 0, 0, -1*currentFactor_; //TODO: implement properly
            }
            // Set current in guidance class and update it
            guidanceEnvironment_.setCurrent(currentToSet_);
            guidanceEnvironment_.updateAll();
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
        currentSaveVector_[0] = currentFactor_;
        currentSaveVector_[1] = currentToSet_[0];
        currentSaveVector_[2] = currentToSet_[1];
        currentSaveVector_[3] = currentToSet_[2];
        currentMap_.insert(std::pair<double, Eigen::VectorXd> (simulationTime_, currentSaveVector_));

        // Save body data
        bodyDataSaveVector_[0] = guidanceEnvironment_.getR();
        bodyDataSaveVector_[1] = guidanceEnvironment_.getVehicleState()[0];
        bodyDataSaveVector_[2] = guidanceEnvironment_.getVehicleState()[1];
        bodyDataSaveVector_[3] = guidanceEnvironment_.getVehicleState()[2];
        bodyDataSaveVector_[4] = guidanceEnvironment_.getVehicleState()[3];
        bodyDataSaveVector_[5] = guidanceEnvironment_.getVehicleState()[4];
        bodyDataSaveVector_[6] = guidanceEnvironment_.getVehicleState()[5];
        bodyDataMap_.insert(std::pair<double, Eigen::VectorXd> (simulationTime_, bodyDataSaveVector_));


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
    Eigen::Vector3d currentToSet_;

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


    // Current factor, TODO: model properly!
    double currentFactor_ = 0* 10E6; //











};



#endif //TUDATBUNDLE_EDTGUIDANCE_H
