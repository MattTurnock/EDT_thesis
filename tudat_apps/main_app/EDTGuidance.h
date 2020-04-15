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
    EDTGuidance(std::string thrustMagnitudeConfig,
            std::string thrustDirectionConfig,
            NamedBodyMap& environmentBodyMap,
            EDTEnvironment& guidanceEnvironment,
            const double throttleSetting = 1 ):
            thrustMagnitudeConfig_(thrustMagnitudeConfig),
            thrustDirectionConfig_(thrustDirectionConfig),
            environmentBodyMap_(environmentBodyMap),
            guidanceEnvironment_(guidanceEnvironment),
            throttleSetting_(throttleSetting){};

    /////////////////////// Update functions based on config settings for custom thrust ///////////////////////////////
    // TODO: Consider updating environment differently (currently done in both magnitude and direction functions)
    // Function to update thrust magnitude based on time and the config used
    void updateThrustMagnitude(const double simulationTime){
        if (thrustMagnitudeConfig_ == "constant"){
            thrustMagnitude_ = thrustMagnitudeConstant_;
        }
        else if (thrustMagnitudeConfig_ == "nominal") {
            guidanceEnvironment_.updateAll();
            thrustVector_ = (guidanceEnvironment_.getCurrent()).cross(guidanceEnvironment_.getMagField());
            thrustMagnitude_ = thrustVector_.norm();
            std::cout << std::to_string(thrustMagnitude_) << std::endl ;
        }
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
            thrustVector_ = (guidanceEnvironment_.getCurrent()).cross(guidanceEnvironment_.getMagField());
            thrustDirection_ = thrustVector_.normalized();
        }
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




protected:
    /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////
    // EDT config variables
    std::string thrustMagnitudeConfig_;
    std::string thrustDirectionConfig_;

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
    double thrustMagnitude_;
    Eigen::Vector3d thrustDirection_;

    // Current vector used for thrust direction guidance
    Eigen::Vector3d currentToSet_;

    // Current factor, TODO: model properly!
    double currentFactor_ = 1E-6; //











};



#endif //TUDATBUNDLE_EDTGUIDANCE_H
