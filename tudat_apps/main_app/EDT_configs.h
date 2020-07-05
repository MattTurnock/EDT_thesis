//
// Created by matt on 29/02/2020. Used to create classes describing each EDT config
//

#ifndef TUDATBUNDLE_EDT_CONFIGS_H
#define TUDATBUNDLE_EDT_CONFIGS_H

//#include "EDTGuidance.h"
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

        // Constructor
        EDTConfig(
//                EDTGuidance& guidanceClass,
                std::string& configType,
                nlohmann::json hoytetherVariables,
                nlohmann::json SRPVariables,
                double length=1,
                double diameterInner=0.1,
                double diameterOuter=0.1,
                double vehicleMass=100,
                double vehicleISP=9999,
                double constantThrustAccel = 0.325E-3,
                Eigen::Vector3d bodyFixedThrustDirection = {1, 0, 0}):
//                guidanceClass_(guidanceClass),
                configType_(configType),
                hoytetherVariables_(hoytetherVariables),
                SRPVariables_(SRPVariables),
                length_(length),
                diameterInner_(diameterInner),
                diameterOuter_(diameterOuter),
                vehicleMass_(vehicleMass),
                vehicleISP_(vehicleISP),
                constantThrustAccel_(constantThrustAccel),
                bodyFixedThrustDirection_(bodyFixedThrustDirection){

            // Normalize body fixed thrust direction
            bodyFixedThrustDirection_.normalize();

            // Set EDT constant bodyMass
            EDTBody->setConstantBodyMass(vehicleMass);

            // Set and calculate hoytether variables from the hoytether variables json
            tetherDiameterInnerSecondary_ = hoytetherVariables_["tetherDiameterInnerSecondary"];
            tetherDiameterOuterSecondary_ = hoytetherVariables_["tetherDiameterOuterSecondary"];
            noPrimaryLines_ = hoytetherVariables_["noPrimaryLines"];
            primaryLineSegmentLength_ = hoytetherVariables_["primaryLineSegmentLength"];
            slackCoefficient_ = hoytetherVariables_["slackCoefficient"];
            primaryLineSeparation_ = hoytetherVariables_["primaryLineSeparation"];
            occultationCoefficient_ = hoytetherVariables_["occultationCoefficient"];

            initialiseHoytetherProperties();

            // Set values for SRP variables from json
            endmassArea1_ = SRPVariables_["endmassArea1"];
            endmassArea2_ = SRPVariables_["endmassArea2"];
            endmassRadiationCoefficient_ = SRPVariables_["endmassRadiationCoefficient"];
            tetherRadiationCoefficient_ = SRPVariables_["tetherRadiationCoefficient"];
            rotationCoefficient_ = SRPVariables_["rotationFactor"];
            /////////////////////// Calculate derived EDT variables /////////////////////////

            // Tether cross sectional areas
            xSecAreaTotal_ = gen::getCircleArea(diameterOuter_); //TODO: make a function that calculates cross sectional area, depending on config type
            xSecAreaInner_ = gen::getCircleArea(diameterInner_);
            xSecAreaDonut_ = gen::getDonutArea(diameterInner_, diameterOuter_);

            // Tether

            // Set some parameters based on config type
            if (configType_ == "CHB"){ // TODO: add other ors for other hoytether types
                // Calculate some SRP related properties, for hoytethers
                setSRPForHoytether();


                // TODO: Calculate conductivity, assuming both inner and outer used for it
            }
            else{
                // SRP CALCULATIONS FOR SINGLELINE DESIGNS
                setSRPForSingleLine();

                //TODO: fill in other config types as if-else
                std::cout << "EDT Config type " << configType_ << " does not exist" << std::endl;
            }


        }

        void initialiseHoytetherProperties(){
            // Calculate number of tether segments h
            noTetherSegments_ = length_/primaryLineSegmentLength_;

            // Calculate number of secondary tether links m
            noSecondaryLinks_ = 2 * noPrimaryLines_ * noTetherSegments_;

            // Calculate secondary line segment length ls
            secondaryLineSegmentLength_ = slackCoefficient_ * sqrt( pow(primaryLineSeparation_, 2) + pow(primaryLineSegmentLength_, 2));

            // Calculate number of primary links nl
            noPrimaryLinks_ = noPrimaryLines_ * noTetherSegments_;

            // Calculate total primary and secondary line length Lp and Ls
            totalPrimaryLineLength_ = length_ * noPrimaryLines_;
            totalSecondaryLineLength_ = noSecondaryLinks_ * secondaryLineSegmentLength_;

            // Calculate total primary and secondary SRP areas Ap and As, as well as effective SRP area Aeff
            totalPrimarySRPArea_ = diameterOuter_ * totalPrimaryLineLength_;
            totalSecondarySRPArea_ = tetherDiameterOuterSecondary_ * totalSecondaryLineLength_;
            effectiveHoytetherSRPArea_ = rotationCoefficient_ * occultationCoefficient_ * (totalPrimarySRPArea_ + totalSecondarySRPArea_);
        }

        void setSRPForHoytether(){
            // Set total effective SRP Area
            totalEffectiveSRPArea_ = endmassArea1_ + endmassArea2_ + effectiveHoytetherSRPArea_;

            // Set total effective radiation pressure coefficient, using a weighted sum
            totalEffectiveRadiationPressureCoefficient_ =
                    ( endmassRadiationCoefficient_*(endmassArea1_ + endmassArea2_) + tetherRadiationCoefficient_ * effectiveHoytetherSRPArea_ ) /
                    (endmassArea1_ + endmassArea2_ + effectiveHoytetherSRPArea_);
        }

        void setSRPForSingleLine(){
            // Set total effective SRP area
            effectiveSingleLineSRPArea_ = rotationCoefficient_ * length_ * diameterOuter_;
            totalEffectiveSRPArea_ = endmassArea1_ + endmassArea2_ + effectiveSingleLineSRPArea_;

            // Set total effective radiation pressure coefficient, using a weighted sum
            totalEffectiveRadiationPressureCoefficient_ =
                    ( endmassRadiationCoefficient_*(endmassArea1_ + endmassArea2_) + tetherRadiationCoefficient_ * effectiveSingleLineSRPArea_ ) /
                    (endmassArea1_ + endmassArea2_ + effectiveSingleLineSRPArea_);
        }

        // Get-set functions for protected variables
        Eigen::Vector3d getBodyFixedThrustDirection( )
        {
            return bodyFixedThrustDirection_;
        }

        double getConstantThrust(){
            return constantThrust_;
        }

        double getEffectiveSRPArea(){
            return totalEffectiveSRPArea_;
        }

        double getEffectiveSRPCoefficient(){
            return totalEffectiveRadiationPressureCoefficient_;
        }

        double getXSecAreaConducting(){
            return xSecAreaConducting_;
        }

    protected:

        /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////

        // Guidance class and config type
//        EDTGuidance& guidanceClass_;
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

        // Create derived EDT properties - cross sectional areas
        double xSecAreaTotal_;
        double xSecAreaInner_;
        double xSecAreaDonut_;
        double xSecAreaConducting_; // Refers to the cross sectional area used for carrying current TODO: Add functions to calculate this

        // Create json hoytether EDT properties
        nlohmann::json hoytetherVariables_;
        double tetherDiameterInnerSecondary_;
        double tetherDiameterOuterSecondary_;
        double noPrimaryLines_;
        double primaryLineSegmentLength_;
        double slackCoefficient_;
        double primaryLineSeparation_;
        double occultationCoefficient_;

        // Create derived EDT properties - hoytether stuff
        double noSecondaryLinks_;
        double noTetherSegments_;
        double secondaryLineSegmentLength_;
        double noPrimaryLinks_;
        double totalPrimaryLineLength_;
        double totalSecondaryLineLength_;
        double totalPrimarySRPArea_;
        double totalSecondarySRPArea_;
        double effectiveHoytetherSRPArea_;


        // Create derived EDT properties - other
        double conductivity_;

        // Create json object and regular variables for SRP properties
        nlohmann::json SRPVariables_;
        double endmassArea1_;
        double endmassArea2_;
        double endmassRadiationCoefficient_;
        double tetherRadiationCoefficient_;
        double rotationCoefficient_;

        // Create derived SRP properties, for direct use in cannonball radiation pressure
        double totalEffectiveSRPArea_;
        double totalEffectiveRadiationPressureCoefficient_;
        double effectiveSingleLineSRPArea_;

        // Initialised body fixed thrust direction
        Eigen::Vector3d bodyFixedThrustDirection_;

    };

}
#endif //TUDATBUNDLE_EDT_CONFIGS_H
