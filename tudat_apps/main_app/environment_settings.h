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
            nlohmann::json ISMFVariables):
            B0EstimationParameters_(B0EstimationParameters),
            phi0_(phi0),
            R0_(R0),
            environmentBodyMap_(environmentBodyMap),
            ISMFVariables_(ISMFVariables){
        //////////// CONSTRUCTOR ////////////////

        // Set values for ISMF values from json data, and calculate the new ones
        BInfMagnitude_ = ISMFVariables_["BInf_tesla"];
        longitudeInf_deg_ = ISMFVariables_["longitudeInf_deg"];
        latitudeInf_deg_ = ISMFVariables_["latitudeInf_deg"];
        longitudeInf_ = gen::deg2rad * longitudeInf_deg_;
        latitudeInf_ = gen::deg2rad * latitudeInf_deg_;

        calculateBinfDirection();
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////// Define Update functions ///////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Update magnetic field vector using base variables, and determining if Parker model or other should be used
    void updateMagField(){
        //TODO: Make some way of determining the magnetic field region. Currently uses one predefined here
        magFieldRegion_ = "Parker";


        if (magFieldRegion_ == "Parker") {
            // Define local magnetic field vector
            double BR = B0_ * cos(phi0_) * pow((R0_ / R_), 2);
            double Bphi = B0_ * sin(phi0_) * pow((R0_ / R_), 1);
            double BzLocal = 0;
            magFieldMaglocal_ << BR, Bphi, BzLocal;

            // Convert magnetic standard local magnetic field to lvlh
            magFieldLvlh_ = gen::MaglocalToLvlh(magFieldMaglocal_);

            // Convert local magfield to inertial and calculate magnitude
            magField_ = gen::LvlhToInertial(magFieldLvlh_, theta_);
            magFieldMagnitude_ = magFieldMaglocal_.norm();
        }

        else if (magFieldRegion_ == "Transitional"){
            //TODO: make transitional magfield info here
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

        // Update B0 using approximate sin relationship (need to convert to proper time, and convert to nanotesla
        simulationTimeYears_ = gen::tudatTime2DecimalYear(simulationTime);
        B0_ = 1E-9 * gen::twoSines(simulationTimeYears_, B0EstimationParameters_);
    }

    // Run all updaters
    void updateAll(double simulationTime=0){
        updateBodyParameters(simulationTime);
        updateMagField();
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

    void setCurrent(Eigen::Vector3d newCurrent){
        current_ = newCurrent;
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

    Eigen::Vector3d getMagFieldInertial(){
        return magField_;
    }

    Eigen::Vector3d getMagFieldLocal(){
        return magFieldMaglocal_;
    }

    double getMagFieldMagnitude(){
        return magFieldMagnitude_;
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




protected:

    /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////
    // Base variables for Parker model of magnetic field (including vector for estimating B0_)
    double B0_;
    double& phi0_;
    double& R0_;

    std::vector<double>& B0EstimationParameters_;

    // Body map for simulation
    NamedBodyMap& environmentBodyMap_;

    /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////


    /////////////////////////// Other set parameters ///////////////////////////////
    // Vector representing thrust (both direction and magnitude), then just direction unit vector, then just magnitude
    Eigen::Vector3d EDTThrust_;
    Eigen::Vector3d EDTThrustDirection_;
    double EDTThrustMagnitude_;

    // Vectors for current
    Eigen::Vector3d current_;
    double currentMagnitude_;

    // Vectors for general + parker magnetic field (needs a local orientation magFieldMaglocal_ and real one magField_
    Eigen::Vector3d magFieldMaglocal_;
    Eigen::Vector3d magFieldLvlh_;
    Eigen::Vector3d magField_;
    double magFieldMagnitude_;
    std::string magFieldRegion_;

    // Data values for interstellar magnetic field
    Eigen::Vector3d BInfDirection_;
    nlohmann::json ISMFVariables_;
    double BInfMagnitude_;
    double longitudeInf_deg_;
    double latitudeInf_deg_;
    double longitudeInf_;
    double latitudeInf_;

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
