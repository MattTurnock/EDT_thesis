//
// Created by matt on 11/03/2020.
//

#ifndef TUDATBUNDLE_ENVIRONMENT_SETTINGS_H
#define TUDATBUNDLE_ENVIRONMENT_SETTINGS_H

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;

namespace environment_custom{

    class MagneticField {
    public:
        123;

        MagneticField()

    private:
        456;
    };



double getMagFieldStrength(double B0, double phi0, double r0, double r){

    double Br_0 = B0*std::cos(phi0);
    double Bphi_0 = B0*std::sin(phi0);

    double BR = Br_0 * std::pow( (r0/r), 2);
    double Bphi = Bphi_0 * (r0/r);

    double Bmag = std::sqrt(std::pow(BR, 2) + std::pow(Bphi, 2));

    return Bmag
}

}

#endif //TUDATBUNDLE_ENVIRONMENT_SETTINGS_H
