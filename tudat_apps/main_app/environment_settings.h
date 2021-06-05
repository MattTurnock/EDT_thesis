////
//// Created by matt on 11/03/2020.
////

#ifndef TUDATBUNDLE_ENVIRONMENT_SETTINGS_H
#define TUDATBUNDLE_ENVIRONMENT_SETTINGS_H

#include "general_functions.h"

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::basic_mathematics;
using namespace tudat::mathematical_constants;



class EDTEnvironment{

public:
    EDTEnvironment(std::vector<double>& B0EstimationParameters,
            double& phi0,
            double& R0,
            NamedBodyMap& environmentBodyMap,
            nlohmann::json simulationVariables):
            B0EstimationParameters_(B0EstimationParameters),
            phi0_(phi0),
            R0_(R0),
            environmentBodyMap_(environmentBodyMap),
            simulationVariables_(simulationVariables){
        //////////// CONSTRUCTOR ////////////////

        // Set values from test variables
        ISMFVariables_ = simulationVariables_["InterstellarMagField"];

        // Set values for ISMF values from json data, and calculate the new ones
        BInfMagnitude_ = ISMFVariables_["BInf_nT"];
        longitudeInf_deg_ = ISMFVariables_["longitudeInf_deg"];
        latitudeInf_deg_ = ISMFVariables_["latitudeInf_deg"];
        longitudeInf_ = gen::deg2rad * longitudeInf_deg_;
        latitudeInf_ = gen::deg2rad * latitudeInf_deg_;
        magFieldTransitionType_ = ISMFVariables_["transitionType"];
        sphericalTransitionDistanceAU_ = ISMFVariables_["sphericalTransitionDistanceAU"];
        sphericalTransitionDistance_ = gen::AU * sphericalTransitionDistanceAU_;

        // Set values for imposed magnetic field
        imposingMagField_ = simulationVariables_["ImposedMagField"]["imposingMagField"];
        double imposedB1 = simulationVariables_["ImposedMagField"]["B1_nT"];
        double imposedB2 = simulationVariables_["ImposedMagField"]["B2_nT"];
        double imposedB3 = simulationVariables_["ImposedMagField"]["B3_nT"];
        imposedMagField_ << imposedB1, imposedB2, imposedB3;

        calculateBinfDirection();
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////// Define Update functions ///////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Update magnetic field vector using base variables, and determining if Parker model or other should be used
    void updateMagField(){

        // Checks if magnetic field is directly imposed, and sets it to that if it is
        if (imposingMagField_){
            magField_ = imposedMagField_;
            magFieldMagnitude_ = imposedMagField_.norm();
        }
        else{
            // Check the transition to type to determine where parker vs ismf magfield is used
            if (magFieldTransitionType_ == "spherical"){

                if (R_ <= sphericalTransitionDistance_){
                    magFieldRegion_ = "Parker";
                }
                else{
                    magFieldRegion_ = "Interstellar";
                }

            }
            else if (magFieldTransitionType_ == "elongated"){
                //TODO: Implement elongated transition zone (not just spherical)
            }
            else{
                std::cout << "Magnetic field transition type is unknown -  " << magFieldTransitionType_ << std::endl;
            }

//        std::cout << "Region: " << magFieldRegion_ << std::endl;
//        std::cout << "Mag strength: " << BInfMagnitude_ << std::endl;
//        std::cout << "Radius: " << R_ << std::endl;
//        std::cout << ""  << std::endl;

            if (magFieldRegion_ == "Parker") {
                // Define local magnetic field vector
                double BR = B0_ * cos(phi0_) * pow((R0_ / R_), 2);
                double Bphi = B0_ * sin(phi0_) * pow((R0_ / R_), 1);
                double BzLocal = 0;
                magFieldMaglocal_ << BR, Bphi, BzLocal;

//                // Convert magnetic standard local magnetic field to lvlh
//                magFieldLvlh_ = gen::MaglocalToLvlh(magFieldMaglocal_);
//
//                // Convert local magfield to inertial and calculate magnitude
//                magField_ = gen::LvlhToInertial(magFieldLvlh_, theta_);
//                magFieldMagnitude_ = magFieldMaglocal_.norm();
                magField_ = gen::LvlhToInertial(magFieldMaglocal_, theta_);
                magFieldMagnitude_ = magFieldMaglocal_.norm();



            }

            else if (magFieldRegion_ == "Transitional"){
                //TODO: make transitional magfield info here (NOTE: No longer required)
            }

            else if (magFieldRegion_ == "Interstellar"){

                // set magnetic field magnitude and vector, directly from data, and from the calculated values
                magFieldMagnitude_ = BInfMagnitude_;
                magField_ = magFieldMagnitude_ * BInfDirection_;



            }
            else{
                std::cout << "Magnetic field region is unknown -  " << magFieldRegion_ << std::endl;
            }
        }
    }

    // Update various body parameters from state, or calculation
    void updateBodyParameters(double simulationTime=0){
        vehicleState_ = environmentBodyMap_["Vehicle"]->getState();
        vehiclePosition_ = environmentBodyMap_["Vehicle"]->getPosition();
        vehicleVelocity_ = environmentBodyMap_["Vehicle"]->getVelocity();
        R_ = vehiclePosition_.norm(); // TODO: Maybe find way to get actual altitude / radius directly

        theta_ = atan2(vehiclePosition_[1], vehiclePosition_[0]);
        if (theta_ < 0){
            theta_ += 2*PI;
        }

        // Update B0 using approximate sin relationship (need to convert to proper time, and result is in nT)
        simulationTimeYears_ = gen::tudatTime2DecimalYear(simulationTime);
        B0_ = gen::twoSines(simulationTimeYears_, B0EstimationParameters_);

    }



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////// Define a few misc functions (eg for creating ISMF variables) ///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Function to calculate the ISMF direction vector, to be used in constructor
    void calculateBinfDirection(){

        // Calculate direction vector components (before normalisation)
        double BInfx_unit = cos(latitudeInf_) * cos(longitudeInf_);
        double BInfy_unit = cos(latitudeInf_) * sin(longitudeInf_);
        double BInfz_unit = sin(latitudeInf_);

        // Create, initialise, and normalise (inplace) direction vector from above components
        BInfDirection_ << BInfx_unit, BInfy_unit, BInfz_unit;
        BInfDirection_.normalize();

    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////// Define some set functions (eg for pitch) ///////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void setCurrentMagnitude(double newCurrentMagnitude){
        currentMagnitude_ = newCurrentMagnitude;
//        std::cout << ":: current magnitude: " << currentMagnitude_ << std::endl;
    }

    void setCurrentDirectionVector(Eigen::Vector3d newCurrentDirectionVector){
        currentDirectionVector_ = newCurrentDirectionVector;
//        std::cout << ":: current direction vector components: " << currentDirectionVector_[0] << currentDirectionVector_[1] << currentDirectionVector_[2] << std::endl;
    }

    void setCurrent(Eigen::Vector3d newCurrent){
        current_ = newCurrent;
//        std::cout << ":: current vector components: " << current_[0] << current_[1] << current_[2] << std::endl;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////// Define simple get functions for protected variables ///////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // TODO: Clean these equations + variables (if not in this class then get rid of them)
    Eigen::Vector3d getEDTThrust(){
        return EDTThrust_;
    }

    Eigen::Vector3d getEDTThrustDirection(){
        return EDTThrustDirection_;
    }

    double getEDTThrustMagnitude(){
        return EDTThrustMagnitude_;
    }

    Eigen::Vector3d getCurrent(){
        return current_;
    }

    double getCurrentMagnitude(){
        return currentMagnitude_;
    }

    Eigen::Vector3d getCurrentDirectionVector(){
        return currentDirectionVector_;
    }

    // Return magnetic field values - both nanotesla returns AND tesla returns

    Eigen::Vector3d getMagFieldInertialnT(){
        return magField_;
    }

    Eigen::Vector3d getMagFieldLocalnT(){
        return magFieldMaglocal_;
    }

    double getMagFieldMagnitudenT(){
        return magFieldMagnitude_;
    }

    Eigen::Vector3d getMagFieldInertial(){
        return 1E-9 * magField_;
    }

    Eigen::Vector3d getMagFieldLocal(){
        return 1E-9 * magFieldMaglocal_;
    }

    double getMagFieldMagnitude(){
        return 1E-9 * magFieldMagnitude_;
    }

    double getB0(){
        return B0_;
    }

    double getPhi0(){
        return phi0_;
    }

    double getR0(){
        return R0_;
    }

    Eigen::Vector6d getVehicleState(){
        return vehicleState_;
    }

    Eigen::Vector3d getVehiclePosition(){
        return vehiclePosition_;
    }

    Eigen::Vector3d getVehicleVelocity(){
        return vehicleVelocity_;
    }

    double getR(){
        return R_;
    }

    double getTheta(){
        return theta_;
    }

    double getSphericalTransitionDistance(){
        return sphericalTransitionDistance_;
    }




protected:

    /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////
    // Base variables for Parker model of magnetic field (including vector for estimating B0_)
    double B0_;
    double& phi0_;
    double& R0_;

    std::vector<double>& B0EstimationParameters_;

    // Body map for simulation
    NamedBodyMap& environmentBodyMap_;

    // General test variables
    nlohmann::json simulationVariables_;

    /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////


    /////////////////////////// Other set parameters ///////////////////////////////
    // Vector representing thrust (both direction and magnitude), then just direction unit vector, then just magnitude
    Eigen::Vector3d EDTThrust_;
    Eigen::Vector3d EDTThrustDirection_;
    double EDTThrustMagnitude_;

    // Variables for general current calcs TODO: Check if thes values actually needed for anything here
    Eigen::Vector3d current_;
    Eigen::Vector3d currentDirectionVector_;
    double currentMagnitude_;


    // Vectors for general + parker magnetic field (needs a local orientation magFieldMaglocal_ and real one magField_
    Eigen::Vector3d magFieldMaglocal_;
    Eigen::Vector3d magFieldLvlh_;
    Eigen::Vector3d magField_;
    double magFieldMagnitude_;
    std::string magFieldRegion_;

    bool imposingMagField_;
    Eigen::Vector3d imposedMagField_;

    // Data values for interstellar magnetic field
    Eigen::Vector3d BInfDirection_;
    nlohmann::json ISMFVariables_;
    double BInfMagnitude_;
    double longitudeInf_deg_;
    double latitudeInf_deg_;
    double longitudeInf_;
    double latitudeInf_;
    std::string magFieldTransitionType_;
    double sphericalTransitionDistanceAU_;
    double sphericalTransitionDistance_;

    // Base vehicle variables. State position and velocity vectors, radius from sun, EDT pitch (or effective pitch) (ie bodyParameters)
    Eigen::Vector6d vehicleState_;
    Eigen::Vector3d vehiclePosition_;
    Eigen::Vector3d vehicleVelocity_;
    double R_;
    double theta_; // Angle of rotation in inertial frame

    // Simulation time in years
    double simulationTimeYears_;

};


#endif //TUDATBUNDLE_ENVIRONMENT_SETTINGS_H
