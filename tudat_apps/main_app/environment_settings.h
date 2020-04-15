////
//// Created by matt on 11/03/2020.
////

#ifndef TUDATBUNDLE_ENVIRONMENT_SETTINGS_H
#define TUDATBUNDLE_ENVIRONMENT_SETTINGS_H

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::basic_mathematics;


class EDTEnvironment{

public:
    EDTEnvironment(const double B0, const double phi0, const double R0, NamedBodyMap& environmentBodyMap):
    B0_(B0), phi0_(phi0), R0_(R0), environmentBodyMap_(environmentBodyMap){};

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////// Define Update functions ///////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Update magnetic field vector using base variables, and determining if Parker model or other should be used
    void updateMagField(){
        //TODO: Add if-else statement using state for parker model, transitional, or interstellar use. Currently only uses Parker
        double Bx = B0_ * cos(phi0_) * pow((R0_/R_), 2);
        double By = B0_ * sin(phi0_) * pow((R0_/R_), 1);
        double Bz = 0;

        magField_ << Bx, By, Bz;
    }

    void updateBodyParameters(){
//        std::cout<< "TEMP updating parameters" << std::to_string(environmentBodyMap_.size()) << std::endl;
        vehicleState_ = environmentBodyMap_["Vehicle"]->getState();
        vehiclePosition_ = environmentBodyMap_["Vehicle"]->getPosition();
        vehicleVelocity_ = environmentBodyMap_["Vehicle"]->getVelocity();
        R_ = vehicleVelocity_.norm();
    }

    // Run all updaters
    void updateAll(){
        updateBodyParameters();
        updateMagField();
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

    Eigen::Vector3d getMagField(){
        return magField_;
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




protected:

    /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////
    // Base variables for Parker model of magnetic field
    double B0_;
    double phi0_;
    double R0_;

    // Body map for simulation
    NamedBodyMap& environmentBodyMap_;

    /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////


    /////////////////////////// Other set parameters ///////////////////////////////
    // Vector representing thrust (both direction and magnitude), then just direction unit vector, then just magnituden
    Eigen::Vector3d EDTThrust_;
    Eigen::Vector3d EDTThrustDirection_;
    double EDTThrustMagnitude_;

    // Vectors for current and magnetic field
    Eigen::Vector3d current_;
    double currentMagnitude_;
    Eigen::Vector3d magField_;
    double magFieldMagnitude_;

    // Base vehicle variables. State position and velocity vectors, radius from sun, EDT pitch (or effective pitch) (ie bodyParameters)
    Eigen::Vector6d vehicleState_;
    Eigen::Vector3d vehiclePosition_;
    Eigen::Vector3d vehicleVelocity_;
    double R_;

};


#endif //TUDATBUNDLE_ENVIRONMENT_SETTINGS_H
