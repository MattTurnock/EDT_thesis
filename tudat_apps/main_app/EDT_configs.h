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
                nlohmann::json configVariables,
                std::string& configType,
                nlohmann::json SRPVariables,
                nlohmann::json materialProperties,
                double vehicleMass=100,
                double vehicleISP=9999,
                double constantThrustAccel = 0.325E-3,
                Eigen::Vector3d bodyFixedThrustDirection = {1, 0, 0}):
                configVariables_(configVariables),
                configType_(configType),
                SRPVariables_(SRPVariables),
                materialProperties_(materialProperties),
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

            /////////////////////// Set SRP parameters (based on config type) /////////////////////////////////////
            if ( (configType_ == "CHB") or (configType_ == "CHTr") ){ // TODO: CHECK DESIGN SPECS - do we use tape tethers or no? IMPORTANT!!!!!!
                // Calculate some SRP related properties, for hoytethers
                setSRPForHoytether();


                // TODO: Calculate conductivity, assuming both inner and outer used for it
            }
            else if ( (configType_ == "AlMB") or (configType_ == "AlMTr") ){
                // SRP CALCULATIONS FOR SINGLELINE DESIGNS
                setSRPForSingleLine();

                //TODO: fill in other config types as if-else
                std::cout << "EDT Config type " << configType_ << " does not exist" << std::endl;
            }
            else{
                std::cout << "Config type: " << configType_ << " Not recognised" << std::endl;
            }

            ///////////////////////// Set current-based parameters (based on config type) ///////////////////////
            if ( (configType_ == "CHB") or (configType_ == "CHTr") ) {
                resistivity_ = getCompositeResistivity(length_, resistivityAl_, resistivityCu_, xSecAreaDonut_, xSecAreaInner_);
            }
            else if ( (configType_ == "AlMB") or (configType_ == "AlMTr") ) {
                resistivity_ = resistivityAl_;
            }
            conductivity_ = getConductivity(resistivity_);
            EDTResistance_ = getResistance(resistivity_, length_, xSecAreaConducting_); // TODO: Check length correct



        }

        ////////////////////////////////////////////// SRP and hoytether properties //////////////////////////////////////

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


        ////////////////////// Conductivity and other current-related EDT properties /////////////////////////////

        // Function to calculate conductivity, from length, area, and total resistance
        void setConductivity(double length, double resistance, double area){
            conductivity_ = length / (resistance*area);
        }

        double getResistance(double resistivity, double length, double area){
            return (resistivity*length)/area;
        }

        // Function to calculate total resistance for a composite tether (assumed as 2 resistors in parallel)
        double getCompositeResistance(double resistance1, double resistance2){
            return pow( (1/resistance1 + 1/resistance2) ,-1);
        }

        // Function to get composite resistivity, using resistance as an intermediary
        double getCompositeResistivity(double length, double resistivity1, double resistivity2, double area1, double area2){
            double Atot = area1 + area2;
            double R1 = getResistance(resistivity1, length, area1);
            double R2 = getResistance(resistivity2, length, area2);
            double Rtot = getCompositeResistance(R1, R2);

            return (Rtot * Atot)/length;
        }

        // Function to calculate conductivity, from resistivity (ie the inverse)
        double getConductivity(double resistivity) {
            return 1/resistivity;
        }

        /////////////////////////// Get-set functions /////////////////////////////////////////////////////////////////////
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

        double getVehicleMass(){
            return vehicleMass_;
        }

    protected:

        /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////

        // Config variables from the main json file
        nlohmann::json configVariables_;

        // Guidance class and config type
//        EDTGuidance& guidanceClass_;
        std::string& configType_;

        // Create initial variables for EDT properties
        double length_ = configVariables_["tetherLength"];
        double diameterInner_ = configVariables_["tetherDiameterInner"];
        double diameterOuter_ = configVariables_["tetherDiameterOuter"];

        /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////
        double vehicleMass_;
        double vehicleISP_;
        double constantThrustAccel_;


        /////////////////////////// Other set parameters ///////////////////////////////
        // Initialise constant thrust value
        double constantThrust_;

        // Create derived EDT properties - cross sectional areas. Inner is Cu, outer is Al
        double xSecAreaTotal_;
        double xSecAreaInner_;
        double xSecAreaDonut_;
        double xSecAreaConducting_; // Refers to the cross sectional area used for carrying current TODO: Add functions to calculate this

        // Create json hoytether EDT properties
        nlohmann::json hoytetherVariables_ = configVariables_["hoytether"];
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
        nlohmann::json materialProperties_;
        double resistivityAl_ = materialProperties_["resistivity"]["Al"];
        double resistivityCu_ = materialProperties_["resistivity"]["Cu"];
        double resistivity_;
        double conductivity_;
        double EDTResistance_;
        double emitterCurrent_ = configVariables_["emitterCurrent"]; // Corresponds to Ic, the current from the e- emitter

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
