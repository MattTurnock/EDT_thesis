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
                vehicleMass_(vehicleMass),
                testmode_(testmode),
                testmodeVehicleAccel_(testmodeVehicleAccel) {

            vehicleThrust_ = vehicleMass_ * testmodeVehicleAccel_;
            vehicleISP_ = 99999;
            bodyFixedThrustDirection_ << 1.0, 0.0, 0.0;
            bodyFixedThrustDirection_.normalize();

            if (testmode) {
                // Retrieve relevant test functions
                thrustMagnitudeFunction_ =
                        std::bind( &EDTConfig::getThrustMagnitudeConstant, this, std::placeholders::_1 );
            }
            else {
                // Do regular stuff TODO: determine what this is

                // Make a little tester for exponentially decreasing thrust TODO: Change this in future tho

                thrustMagnitudeFunction_ =
                        std::bind( &EDTConfig::getThrustMagnitude, this, std::placeholders::_1 );

            }

            // Set thrust direction and magnitude using custom settings
            thrustMagnitudeSettings =
                    std::make_shared< FromFunctionThrustMagnitudeSettings >(
                            thrustMagnitudeFunction_,
                            [ = ]( const double ){ return vehicleISP_; } );

            thrustDirectionGuidanceSettings =
                    std::make_shared<ThrustDirectionFromStateGuidanceSettings>("Sun", true, false);


            // Set EDT constant bodyMass
            EDTBody->setConstantBodyMass(vehicleMass);
        }

        // Thrust and direction return functions
        double getThrustMagnitudeConstant(const double currentTime){

            return vehicleThrust_;
        }

        double getThrustMagnitude(const double currentTime){
            return (1.0E8/currentTime)*vehicleThrust_;
        }

        Eigen::Vector3d getBodyFixedThrustDirection( )
        {
            return bodyFixedThrustDirection_;
        }

        



    protected:
        // General values
        double vehicleMass_;
        bool testmode_;
        double testmodeVehicleAccel_;
        double vehicleThrust_;
        double vehicleISP_;
        Eigen::Vector3d bodyFixedThrustDirection_;

        // Binded functions
        std::function< void( const double ) > updateFunction_;
        std::function< Eigen::Vector3d( ) > thrustDirectionFunction_;
        std::function< double(const double) > thrustMagnitudeFunction_;

    };

    


}
#endif //TUDATBUNDLE_EDT_CONFIGS_H
