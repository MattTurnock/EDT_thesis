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
using namespace tudat::mathematical_constants;


class EDTEnvironment{

public:
    EDTEnvironment(double& B0, double& phi0, double& R0, NamedBodyMap& environmentBodyMap):
    B0_(B0), phi0_(phi0), R0_(R0), environmentBodyMap_(environmentBodyMap){};

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////// Define Update functions ///////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Update magnetic field vector using base variables, and determining if Parker model or other should be used
    void updateMagField(){
        //TODO: Add if-else statement using state for parker model, transitional, or interstellar use. Currently only uses Parker

        // Define local magnetic field vector
        double BxLocal = B0_ * cos(phi0_) * pow((R0_/R_), 2);
        double ByLocal = B0_ * sin(phi0_) * pow((R0_/R_), 1);
        double BzLocal = 0;
//        std::cout << "B0_ = " << B0_ << "    " << "phi0_ = " << phi0_ << "    " << "R0_ = " << R0_ << "    " << "R_ = " << R_ << std::endl;
//        std::cout << "cos(phi0) = " << cos(phi0_) << "    sin(phi0) = " << sin(phi0_) << "    R0/R = " << (1.496E11)/(R_) << "    (R0/R)^2 = " << pow((R0_/R_), 2) << std::endl;
        magFieldLocal_ << BxLocal, ByLocal, BzLocal;

        magFieldMagnitude_ = magFieldLocal_.norm();

        // Convert local to inertial magnetic field vector
        double Bx = BxLocal * sin(theta_) + ByLocal * cos(theta_);
        double By = BxLocal * cos(theta_) - ByLocal * sin(theta_);
        double Bz = 0;
        magField_ << Bx, By, Bz;
        //TODO: delete me down there
//        environmentBodyMap_["Vehicle"].rot
//        std::cout << magFieldLocal_ << std::endl;
//        std::cout << "--------------" << std::endl;
//        std::cout << magField_ << std::endl;
//        std::cout << theta_ << std::endl;

    }

    void updateBodyParameters(){
        vehicleState_ = environmentBodyMap_["Vehicle"]->getState();
        vehiclePosition_ = environmentBodyMap_["Vehicle"]->getPosition();
        vehicleVelocity_ = environmentBodyMap_["Vehicle"]->getVelocity();
        R_ = vehiclePosition_.norm();

        theta_ = atan2(vehiclePosition_[1], vehiclePosition_[0]);
        if (theta_ < 0){
            theta_ += 2*PI;
        }
    }

//    void saveMagFieldProperties(){
////        magFieldSaveVector_ =
//    }

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
        return magFieldLocal_;
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
    // Base variables for Parker model of magnetic field
    double& B0_;
    double& phi0_;
    double& R0_;

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

    // Vectors for magnetic field (needs a local orientation magFieldLocal_ and real one magField_
    Eigen::Vector3d magFieldLocal_;
    Eigen::Vector3d magField_;
    double magFieldMagnitude_;

    // Base vehicle variables. State position and velocity vectors, radius from sun, EDT pitch (or effective pitch) (ie bodyParameters)
    Eigen::Vector6d vehicleState_;
    Eigen::Vector3d vehiclePosition_;
    Eigen::Vector3d vehicleVelocity_;
    double R_;
    double theta_; // Angle of rotation in inertial frame

};


#endif //TUDATBUNDLE_ENVIRONMENT_SETTINGS_H
