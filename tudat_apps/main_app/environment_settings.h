//
// Created by matt on 11/03/2020.
//

#ifndef TUDATBUNDLE_ENVIRONMENT_SETTINGS_H
#define TUDATBUNDLE_ENVIRONMENT_SETTINGS_H

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;

namespace EDT_enviro{

    // Class containing properties of the magnetic field and ionosphere
    class MagneticField {
    public:
        MagneticField(double B0, double Bphi0, double R0):
        B0_(B0),
        Bphi_(Bphi0),
        BR_(R0){}

        void updateMagField(double R){

        }


        double getB(){
            return B_;
        }
        double getBphi() {
            return Bphi_;
        }
        double getBR(){
            return BR_;
        }


    protected:
        double B_;
        double Bphi_;
        double BR_;
        double B0_;
        double Bphi0_;
        double BR0_;

    };



//double getMagFieldStrength(double B0, double phi0, double r0, double r){

//    double Br_0 = B0*std::cos(phi0);
//    double Bphi_0 = B0*std::sin(phi0);

//    double BR = Br_0 * std::pow( (r0/r), 2);
//    double Bphi = Bphi_0 * (r0/r);

//    double Bmag = std::sqrt(std::pow(BR, 2) + std::pow(Bphi, 2));

//    return Bmag
//}

//}

#endif //TUDATBUNDLE_ENVIRONMENT_SETTINGS_H
