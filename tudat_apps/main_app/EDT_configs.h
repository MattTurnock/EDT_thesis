//
// Created by matt on 29/02/2020. Used to create classes describing each EDT config
//

#ifndef TUDATBUNDLE_EDT_CONFIGS_H
#define TUDATBUNDLE_EDT_CONFIGS_H

#include "EDTGuidance.h"

using namespace tudat;
using namespace tudat::simulation_setup;

namespace EDTs {

    class EDTConfig {
    public:
        // Create some variables that will be used in the constructor
        std::shared_ptr<Body> EDTBody = std::make_shared<simulation_setup::Body>();
        std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings;
        std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings;

        // Constructor
        EDTConfig(
                EDTGuidance& guidanceClass,
                double vehicleMass=100,
                double vehicleISP=9999,
                double constantThrustAccel = 0.325E-3,
                Eigen::Vector3d bodyFixedThrustDirection = {1, 0, 0}):
                guidanceClass_(guidanceClass),
                vehicleMass_(vehicleMass),
                vehicleISP_(vehicleISP),
                constantThrustAccel_(constantThrustAccel),
                bodyFixedThrustDirection_(bodyFixedThrustDirection){

            // Set constant thrust based on mass and accel, and set in guidance class
            constantThrust_ = vehicleMass_ * constantThrustAccel_;
            guidanceClass_.setThrustMagnitudeConstant(constantThrust_);

            // Normalize body fixed thrust direction
            bodyFixedThrustDirection_.normalize();

            // Set EDT constant bodyMass
            EDTBody->setConstantBodyMass(vehicleMass);

            // Do initialisation update of guidance settings
            updateGuidanceSettings();
        }

        void updateGuidanceSettings(){
            // Bind functions properly
            thrustMagnitudeFunction_ = std::bind( &EDTGuidance::getThrustMagnitude, guidanceClass_, std::placeholders::_1 );
            thrustDirectionFunction_ = std::bind( &EDTGuidance::getThrustDirection, guidanceClass_, std::placeholders::_1);

            // Set thrust direction and magnitude using custom settings
            thrustMagnitudeSettings =
                    std::make_shared< FromFunctionThrustMagnitudeSettings >(
                            thrustMagnitudeFunction_,
                            [ = ]( const double ){ return vehicleISP_; });

            thrustDirectionGuidanceSettings =
                    std::make_shared<CustomThrustDirectionSettings>(
                            thrustDirectionFunction_);
        }

        Eigen::Vector3d getBodyFixedThrustDirection( )
        {
            return bodyFixedThrustDirection_;
        }

        double getConstantThrust(){
            return constantThrust_;
        }

    protected:

        /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////

        // Guidance class
        EDTGuidance& guidanceClass_;

        /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////
        double vehicleMass_;
        double vehicleISP_;
        double constantThrustAccel_;


        /////////////////////////// Other set parameters ///////////////////////////////
        // Initialise constant thrust value
        double constantThrust_;

        // Initialised body fixed thrust direction
        Eigen::Vector3d bodyFixedThrustDirection_;

        // Custom Thrust standard functions
        std::function< double(const double) > thrustMagnitudeFunction_;
        std::function< Eigen::Vector3d(const double) > thrustDirectionFunction_;

    };

}
#endif //TUDATBUNDLE_EDT_CONFIGS_H
