//
// Created by matt on 29/02/2020. Used to create classes describing each EDT config
//

#ifndef TUDATBUNDLE_EDT_CONFIGS_H
#define TUDATBUNDLE_EDT_CONFIGS_H

using namespace tudat;
using namespace tudat::simulation_setup;

namespace EDTs {

//    class EDTConfigTEST {
//    public:
//        std::shared_ptr<Body> EDTBody = std::make_shared<simulation_setup::Body>();
//        // Constructor
//        EDTConfigTEST(
//                double vehicleMass,
//                double vehicleAccel) :
//                vehicleMass_(vehicleMass) {
//            // Constructor function
//
//            // Set EDT vehicle mass
//            EDTBody->setConstantBodyMass(vehicleMass);
//
//            // Set EDT constant thrust based on desired accel. Also specific impulse is maximum
//            double vehicleThrust = vehicleMass * vehicleAccel;
//            double specificImpulse = 999999;
//
//            // Create some settings variables to be changed at class creation
//            std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
//                    std::make_shared< ThrustDirectionFromStateGuidanceSettings >( "Sun", true, false );
//            std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings =
//                    std::make_shared< ConstantThrustMagnitudeSettings >( vehicleThrust, specificImpulse );
//        }

//    private:
//        // Set values of inputs as private
//        double vehicleMass_;
//    };

    class EDTConfig {
    public:
        // Create some variables that will be used in the constructor
        std::shared_ptr<Body> EDTBody = std::make_shared<simulation_setup::Body>();
        std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings;
        std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings;

        // Constructor
        EDTConfig(
                double vehicleMass,
                bool testmode=false,
                double testmodeVehicleAccel = 0.325E-3) : //TODO: change this to something reasonable Some settings: 0 for free flight, 0.325E-3 for B0 (ie at Earth)
                vehicleMass_(vehicleMass) {

            if (testmode){

                // Set EDT constant thrust based on desired accel. Also specific impulse is maximum
                double vehicleThrust = vehicleMass * testmodeVehicleAccel;
                double specificImpulse = 999999;

                // Set thrust direction to prograde, and constant thrust magnitude
                thrustDirectionGuidanceSettings =
                        std::make_shared< ThrustDirectionFromStateGuidanceSettings >( "Sun", true, true );
                thrustMagnitudeSettings =
                        std::make_shared< ConstantThrustMagnitudeSettings >( vehicleThrust, specificImpulse );
            }
            else{
                // Do regular stuff TODO: determine what this is
            }

            // Set EDT constant bodyMass
            EDTBody->setConstantBodyMass(vehicleMass);

        }

    private:
        // Set values of inputs as private
        double vehicleMass_;
    };



}
#endif //TUDATBUNDLE_EDT_CONFIGS_H
