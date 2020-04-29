//
// Created by matt on 29/02/2020. Used to create classes describing each EDT config
//

#ifndef TUDATBUNDLE_EDT_CONFIGS_H
#define TUDATBUNDLE_EDT_CONFIGS_H

#include "EDTGuidance.h"
#include "general_functions.h"

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::mathematical_constants;
using namespace gen;

namespace EDTs {

    class EDTConfig {
    public:
        // Create EDT body and thrust guidance settings variables
        std::shared_ptr<Body> EDTBody = std::make_shared<simulation_setup::Body>();
        std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings;
        std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings;



        // Constructor
        EDTConfig(
                EDTGuidance& guidanceClass,
                std::string& configType,
                double length=1,
                double diameterInner=0.1,
                double diameterOuter=0.1,
                double vehicleMass=100,
                double vehicleISP=9999,
                double constantThrustAccel = 0.325E-3,
                Eigen::Vector3d bodyFixedThrustDirection = {1, 0, 0}):
                guidanceClass_(guidanceClass),
                configType_(configType),
                length_(length),
                diameterInner_(diameterInner),
                diameterOuter_(diameterOuter),
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


            // Calculate derived EDT variables
            areaTotal_ = gen::getCircleArea(diameterOuter_); //TODO: make as a hoytether config, with many cables
            areaInner_ = gen::getCircleArea(diameterInner_);
            areaDonut_ = gen::getDonutArea(diameterInner_, diameterOuter_);

            // Set some parameters based on config type
            if (configType_ == "CHB"){
                // Calculate conductivity, assuming both inner and outer used for it
            }
            else{
                //TODO: fill in other config types as if-else
                std::cout << "EDT Config type " << configType_ << " does not exist" << std::endl;
            }


        }

        void updateGuidanceSettings(){
            // Bind functions properly
            thrustMagnitudeFunction_ = std::bind( &EDTGuidance::getThrustMagnitude, &guidanceClass_, std::placeholders::_1 );
            thrustDirectionFunction_ = std::bind( &EDTGuidance::getThrustDirection, &guidanceClass_, std::placeholders::_1);

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

        EDTGuidance getGuidanceClass(){
            return guidanceClass_;
        }

    protected:

        /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////

        // Guidance class and config type
        EDTGuidance& guidanceClass_;
        std::string& configType_;

        // Create initial variables for EDT properties
        double length_;
        double diameterInner_;
        double diameterOuter_;

        /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////
        double vehicleMass_;
        double vehicleISP_;
        double constantThrustAccel_;


        /////////////////////////// Other set parameters ///////////////////////////////
        // Initialise constant thrust value
        double constantThrust_;

        // Create derived EDT properties
        double areaTotal_;
        double areaInner_;
        double areaDonut_;
        double conductivity_;

        // Initialised body fixed thrust direction
        Eigen::Vector3d bodyFixedThrustDirection_;

        // Custom Thrust standard functions
        std::function< double(const double) > thrustMagnitudeFunction_;
        std::function< Eigen::Vector3d(const double) > thrustDirectionFunction_;

        // Conductivity values for various materials

    };

}
#endif //TUDATBUNDLE_EDT_CONFIGS_H
