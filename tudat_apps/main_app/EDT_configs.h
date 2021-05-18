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
                nlohmann::json simulationVariables,
                nlohmann::json configVariables,
                std::string& configType,
                nlohmann::json SRPVariables,
                nlohmann::json materialProperties,
                double constantThrustAccel = 0.325E-3,
                Eigen::Vector3d bodyFixedThrustDirection = {1, 0, 0}):
                simulationVariables_(simulationVariables),
                configVariables_(configVariables),
                configType_(configType),
                SRPVariables_(SRPVariables),
                materialProperties_(materialProperties),
                constantThrustAccel_(constantThrustAccel),
                bodyFixedThrustDirection_(bodyFixedThrustDirection){



            // Normalize body fixed thrust direction
            bodyFixedThrustDirection_.normalize();

            // Set initial EDT properties
            length_ = configVariables_["tetherLength"];
            tetherDiameterOuter_ = configVariables_["tetherDiameter"];
            tetherAreaRatio_ = configVariables["tetherAreaRatio"];
            imposedAreaBool_ = configVariables_["imposedAreaBool"];
            imposedArea_ = configVariables_["imposedArea"];

            tetherDiameterSecondary_ = hoytetherVariables_["tetherDiameterSecondary"];
            tetherAreaRatioSecondary_ = hoytetherVariables_["tetherAreaRatioSecondary"];



            // Set and calculate hoytether variables from the hoytether variables json
            noPrimaryLines_ = hoytetherVariables_["noPrimaryLines"];
            primaryLineSegmentLengthRatio_ = hoytetherVariables_["primaryLineSegmentLengthRatio"];
            slackCoefficient_ = hoytetherVariables_["slackCoefficient"];
            primaryLineSeparationRatio_ = hoytetherVariables_["primaryLineSeparationRatio"];
            SRPOccultationCoefficient_ = hoytetherVariables_["SRPOccultationCoefficient"];
            endmassMass1_ = simulationVariables_["scConfigs"]["endmassMass1"];
            endmassMass2_ = simulationVariables_["scConfigs"]["endmassMass2"];



            // Set values for SRP variables from json
            endmassArea1_ = SRPVariables_["endmassArea1"];
            endmassArea2_ = SRPVariables_["endmassArea2"];
            endmassRadiationCoefficient_ = SRPVariables_["endmassRadiationCoefficient"];
            tetherRadiationCoefficient_ = SRPVariables_["tetherRadiationCoefficient"];
            SRPRotationCoefficient_ = SRPVariables_["SRPRotationCoefficient"];
            useSRP_ = SRPVariables_["useSRP"];

            // Set some material properties
            imposedConductivityBool_ = materialProperties_["imposedConductivityBool"];
            imposedConductivity_ = materialProperties_["imposedConductivity"];
            resistivityAl_ = materialProperties_["resistivity"]["Al"];
            resistivityCu_ = materialProperties_["resistivity"]["Cu"];
            emitterCurrentmA_ = configVariables_["emitterCurrentmA"] ;
            emitterCurrent_ = emitterCurrentmA_ / 1000;




            //////////////////////// Set diameters and Areas /////////////////////
            // Use overall diameter to find inner and outer areas
            double TempArea = gen::calculateCircleArea(tetherDiameterOuter_);
            tetherAreaInner_ = TempArea * (1 - tetherAreaRatio_);
            tetherAreaOuter_ = TempArea * tetherAreaRatio_;

            double TempAreaSecondary = gen::calculateCircleArea(tetherDiameterSecondary_);
            tetherAreaInnerSecondary_ = TempAreaSecondary * (1 - tetherAreaRatioSecondary_);
            tetherAreaOuterSecondary_ = TempAreaSecondary * tetherAreaRatioSecondary_;

            // Use inner and outer areas to calculate relevant diameters
            tetherDiameterOuter_ = gen::calculateCircleDiameter(tetherAreaInner_ + tetherAreaOuter_);
            tetherDiameterInner_ = gen::calculateCircleDiameter(tetherAreaInner_);
            tetherDiameterOuterSecondary_ = gen::calculateCircleDiameter(tetherAreaInnerSecondary_ + tetherAreaOuterSecondary_);
            tetherDiameterInnerSecondary_ = gen::calculateCircleDiameter(tetherAreaInnerSecondary_);

            /////////////////////// Calculate geometric variables based on config /////////////////////////

            if ( (configType_ == "CHB") or (configType_ == "CHTr")){

                // Hoytether based solutions have their properties initialised, and other derived properties calculated
                initialiseHoytetherProperties();
                xSecAreaTotal_ = xSecAreaHoytetherPrimaryTotal_;
                xSecAreaInner_ = xSecAreaHoytetherPrimaryInnerTotal_;
                xSecAreaOuter_ = xSecAreaHoytetherPrimaryOuterTotal_;




            }
            else{
                // TODO: Add some stuff here if we use the other types
            }

            if (imposedAreaBool_){
                xSecAreaTotal_ = imposedArea_;
                xSecAreaOuter_ = imposedArea_;
                xSecAreaInner_ = 0;
            }

            // Tether cross sectional areas
//            xSecAreaTotal_ = gen::getCircleArea(tetherAreaOuter_); //TODO: make a function that calculates cross sectional area, depending on config type
//            xSecAreaInner_ = gen::getCircleArea(tetherAreaInner_);
//            xSecAreaDonut_ = gen::getDonutArea(tetherAreaInner_, tetherAreaOuter_);


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


            ///////////////////////// Set current-based parameters, and other physical parameters (based on config type) ///////////////////////
            if ( (configType_ == "CHB") or (configType_ == "CHTr") ) {
                resistivity_ = calculateCompositeResistivity(length_, resistivityAl_, resistivityCu_, xSecAreaOuter_,
                                                             xSecAreaInner_);
            }
            else if ( (configType_ == "AlMB") or (configType_ == "AlMTr") ) {
                resistivity_ = resistivityAl_;
            }


            if (imposedConductivityBool_){
                resistivity_ = 1.0/imposedConductivity_;
            }

            xSecAreaConducting_ = xSecAreaTotal_;
            conductivity_ = resistivityToConductivity(resistivity_);
            EDTResistance_ = calculateResistance(resistivity_, length_, xSecAreaConducting_); // TODO: Check length correct
//            std::cout << "Vehicle conducting area: " << xSecAreaConducting_ << std::endl;
//            std::cout << "Hoytether area primary total (should be same as above) : " << xSecAreaHoytetherPrimaryTotal_ << std::endl;


            // Calculate and set EDT constant bodyMass
            vehicleMass_ = tetherMass_ + endmassMass1_ + endmassMass2_;
            EDTBody->setConstantBodyMass(vehicleMass_);
//            EDTBody->setConstantBodyMass(100);
//            std::cout << "Vehicle masses (total, endmass1, endmass2, tether mass): " << vehicleMass_ << ", " << endmassMass1_ << ", " << endmassMass2_ << ", " << tetherMass_ << std::endl;

        }

        ////////////////////////////////////////////// SRP and hoytether properties //////////////////////////////////////

        void initialiseHoytetherProperties(){
            // Calculate number of tether segments h
            primaryLineSegmentLength_ = length_ * primaryLineSegmentLengthRatio_;
            noTetherSegments_ = length_/primaryLineSegmentLength_;

            // Calculate number of secondary tether links m
            noSecondaryLinks_ = 2 * noPrimaryLines_ * noTetherSegments_;
//            std::cout << "Number secondary links: " << noSecondaryLinks_ << std::endl;

            // Calculate secondary line segment length ls
            primaryLineSeparation_ = tetherDiameterOuter_ * primaryLineSeparationRatio_;
            secondaryLineSegmentLength_ = slackCoefficient_ * sqrt( pow(primaryLineSeparation_, 2) + pow(primaryLineSegmentLength_, 2));

            // Calculate number of primary links nl
            noPrimaryLinks_ = noPrimaryLines_ * noTetherSegments_;

            // Calculate total primary and secondary line length Lp and Ls
            totalPrimaryLineLength_ = length_ * noPrimaryLines_;
            totalSecondaryLineLength_ = noSecondaryLinks_ * secondaryLineSegmentLength_;

            // Calculate total primary and secondary SRP areas Ap and As, as well as effective SRP area Aeff
            totalPrimarySRPArea_ = tetherDiameterOuter_ * totalPrimaryLineLength_;
            totalSecondarySRPArea_ = tetherDiameterOuterSecondary_ * totalSecondaryLineLength_;
            effectiveHoytetherSRPArea_ = SRPRotationCoefficient_ * SRPOccultationCoefficient_ * (totalPrimarySRPArea_ + totalSecondarySRPArea_);

            // Calculate the total cross sectional areas of primary and secondary lines
            xSecAreaHoytetherPrimaryInnerTotal_ = noPrimaryLines_ * tetherAreaInner_;
            xSecAreaHoytetherPrimaryOuterTotal_  = noPrimaryLines_ * tetherAreaOuter_;
            xSecAreaHoytetherPrimaryTotal_ = xSecAreaHoytetherPrimaryInnerTotal_ + xSecAreaHoytetherPrimaryOuterTotal_;

            double xSecAreaHoytetherSecondaryInnerTotal = noSecondaryLinks_ * tetherAreaInner_;
            double xSecAreaHoytetherSecondaryOuterTotal  = noSecondaryLinks_ * tetherAreaOuter_;
            double xSecAreaHoytetherSecondaryTotal = xSecAreaHoytetherSecondaryInnerTotal + xSecAreaHoytetherSecondaryOuterTotal;


            // Mass calculation
            double AlDensity = materialProperties_["density"]["Al"];
            double CuDensity = materialProperties_["density"]["Cu"];
            compositeMaterialDensity_ = calculateCompositeDensity(AlDensity, CuDensity, tetherAreaOuter_, tetherAreaInner_);

            double volumeHoytetherPrimaryTotal = xSecAreaHoytetherPrimaryTotal_ * length_;
            double volumeHoytetherSecondaryTotal = xSecAreaHoytetherSecondaryTotal * secondaryLineSegmentLength_;
            double totalMaterialVolume = volumeHoytetherPrimaryTotal + volumeHoytetherSecondaryTotal;

//            std::cout << "Volumes: " << volumeHoytetherPrimaryTotal << ",  " << volumeHoytetherSecondaryTotal << std::endl;

            tetherMass_ = totalMaterialVolume * compositeMaterialDensity_;
//            std::cout << "Hoytether area conducting: " << xSecAreaConducting_ << std::endl;

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
            effectiveSingleLineSRPArea_ = SRPRotationCoefficient_ * length_ * tetherAreaOuter_;
            totalEffectiveSRPArea_ = endmassArea1_ + endmassArea2_ + effectiveSingleLineSRPArea_;

            if (useSRP_) {
                // Set total effective radiation pressure coefficient, using a weighted sum
                totalEffectiveRadiationPressureCoefficient_ =
                        (endmassRadiationCoefficient_ * (endmassArea1_ + endmassArea2_) +
                         tetherRadiationCoefficient_ * effectiveSingleLineSRPArea_) /
                        (endmassArea1_ + endmassArea2_ + effectiveSingleLineSRPArea_);
            }
            else {
                totalEffectiveRadiationPressureCoefficient_ = 0;
            }

        }


        ////////////////////// Conductivity and other current-related EDT properties /////////////////////////////

        // Function to calculate conductivity, from length, area, and total resistance
        void calculateConductivity(double length, double resistance, double area){
            conductivity_ = length / (resistance*area);
        }

        double calculateResistance(double resistivity, double length, double area){
            return (resistivity*length)/area;
        }

        // Function to calculate total resistance for a composite tether (assumed as 2 resistors in parallel)
        double calculateCompositeResistance(double resistance1, double resistance2){
            return pow( (1/resistance1 + 1/resistance2) ,-1);
        }

        // Function to get composite resistivity, using resistance as an intermediary
        double calculateCompositeResistivity(double length, double resistivity1, double resistivity2, double area1, double area2){
            double Atot = area1 + area2;
            double R1 = calculateResistance(resistivity1, length, area1);
            double R2 = calculateResistance(resistivity2, length, area2);
            double Rtot = calculateCompositeResistance(R1, R2);

            return (Rtot * Atot)/length;
        }

        // Function to calculate conductivity, from resistivity (ie the inverse)
        double resistivityToConductivity(double resistivity) {
            return 1/resistivity;
        }

        // Function to calculate composite density for a cable
        double calculateCompositeDensity(double material1Density, double material2Density, double material1Area, double material2Area){
            double compositeDensity = (material1Density*material1Area + material2Density*material2Area)/(material1Area + material2Area);

            return compositeDensity;
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

        double getEffectiveRadiationPressureCoefficient(){
            return totalEffectiveRadiationPressureCoefficient_;
        }

        double getXSecAreaConducting(){
            return xSecAreaConducting_;
        }

        double getVehicleMass(){
            return vehicleMass_;
        }

        double getVehicleConductivity(){
            return conductivity_;
        }

        double getEmitterCurrent(){
            return emitterCurrent_;
        }

        std::string getConfigType(){
            return configType_;
        }

        double getTetherLength(){
            return length_;
        }

        /// Newly added ones for config saving ///
        double getResistivity(){
            return resistivity_;
        }

        double getEDTResistance(){
            return EDTResistance_;
        }

        double getTetherMass(){
            return tetherMass_;
        }

        double getCompositeMaterialDensity(){
            return compositeMaterialDensity_;
        }

        double getTetherAreaInner(){
            return tetherAreaInner_;
        }

        double getTetherAreaOuter(){
            return tetherAreaOuter_;
        }

        double getTetherAreaInnerSecondary(){
            return tetherAreaInnerSecondary_;
        }

        double getTetherAreaOuterSecondary(){
            return tetherAreaOuterSecondary_;
        }

        double getTetherDiameterInner(){
            return tetherDiameterInner_;
        }

        double getTetherDiameterOuter(){
            return tetherDiameterOuter_;
        }

        double getTetherDiameterInnerSecondary(){
            return tetherDiameterInnerSecondary_;
        }

        double getTetherDiameterOuterSecondary(){
            return tetherDiameterOuterSecondary_;
        }

        double getXSecAreaTotal(){
            return xSecAreaTotal_;
        }

        double getXSecAreaInner(){
            return xSecAreaInner_;
        }

        double getXSecAreaOuter(){
            return xSecAreaOuter_;
        }

        double getNoTetherSegments(){
            return noTetherSegments_;
        }

        double getNoPrimaryLinks(){
            return noPrimaryLinks_;
        }

        double getNoSecondaryLinks(){
            return noSecondaryLinks_;
        }

        double getTotalPrimaryLineLength(){
            return totalPrimaryLineLength_;
        }

        double getTotalSecondaryLineLength(){
            return totalSecondaryLineLength_;
        }

        double getPrimaryLineSeparation(){
            return primaryLineSeparation_;
        }

        double getPrimaryLineSegmentLength(){
            return primaryLineSegmentLength_;
        }

        double getSecondaryLineSegmentLength(){
            return secondaryLineSegmentLength_;
        }

        double getTotalPrimarySRPArea(){
            return totalPrimarySRPArea_;
        }

        double getTotalSecondarySRPArea(){
            return totalSecondarySRPArea_;
        }

        double getEffectiveHoytetherSRPArea(){
            return effectiveHoytetherSRPArea_;
        }



    protected:

        /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////

        // Config variables from the main json file
        nlohmann::json configVariables_;
        nlohmann::json simulationVariables_;

        // Guidance class and config type
//        EDTGuidance& guidanceClass_;
        std::string& configType_;

        // Create initial variables for EDT properties
        double length_;
//        double tetherDiameter_;
        double tetherAreaRatio_;
        double tetherAreaInner_;
        double tetherAreaOuter_;
        bool imposedAreaBool_;
        double imposedArea_;
        double tetherDiameterInner_;
        double tetherDiameterOuter_; // Is set from the above areas, as is inner


        /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////
        double vehicleMass_;
        double constantThrustAccel_;


        /////////////////////////// Other set parameters ///////////////////////////////
        // Initialise constant thrust value
        double constantThrust_;

        // boolean for using SRP
        bool useSRP_;

        // Create derived EDT properties - cross sectional areas. Inner is Cu, outer is Al
        double xSecAreaTotal_;
        double xSecAreaInner_;
        double xSecAreaOuter_;
        double xSecAreaConducting_; // Refers to the cross sectional area used for carrying current TODO: Add functions to calculate this


        // Create json hoytether EDT properties
        nlohmann::json hoytetherVariables_ = configVariables_["hoytether"];
        double tetherDiameterSecondary_;
        double tetherAreaRatioSecondary_;
        double tetherAreaInnerSecondary_;
        double tetherAreaOuterSecondary_;
        double tetherDiameterInnerSecondary_;
        double tetherDiameterOuterSecondary_;
        double noPrimaryLines_;
        double primaryLineSegmentLengthRatio_;
        double primaryLineSegmentLength_;
        double slackCoefficient_;
        double primaryLineSeparationRatio_;
        double primaryLineSeparation_;
        double SRPOccultationCoefficient_;

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

        double xSecAreaHoytetherPrimaryInnerTotal_;
        double xSecAreaHoytetherPrimaryOuterTotal_;
        double xSecAreaHoytetherPrimaryTotal_;


        // Create derived EDT properties - other
        nlohmann::json materialProperties_;
        bool imposedConductivityBool_;
        double imposedConductivity_;
        double resistivityAl_;
        double resistivityCu_;
        double resistivity_;
        double conductivity_;
        double EDTResistance_;
        double emitterCurrentmA_; // Corresponds to Ic, the current from the e- emitter
        double emitterCurrent_ ;
        double compositeMaterialDensity_;
        double tetherMass_;
        double endmassMass1_;
        double endmassMass2_;


        // Create json object and regular variables for SRP properties
        nlohmann::json SRPVariables_;
        double endmassArea1_;
        double endmassArea2_;
        double endmassRadiationCoefficient_;
        double tetherRadiationCoefficient_;
        double SRPRotationCoefficient_;

        // Create derived SRP properties, for direct use in cannonball radiation pressure
        double totalEffectiveSRPArea_;
        double totalEffectiveRadiationPressureCoefficient_;
        double effectiveSingleLineSRPArea_;

        // Initialised body fixed thrust direction
        Eigen::Vector3d bodyFixedThrustDirection_;

    };

}
#endif //TUDATBUNDLE_EDT_CONFIGS_H
