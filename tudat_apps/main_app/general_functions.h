//
// Created by matt on 19/04/2020.
//

#ifndef TUDATBUNDLE_GENERAL_FUNCTIONS_H
#define TUDATBUNDLE_GENERAL_FUNCTIONS_H

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h>

using namespace tudat::mathematical_constants;

namespace gen {
    double getCircleArea(double diameter){
        double radius = 0.5*diameter;
        double area =  PI * std::pow(radius, 2);

        return area;
    }

    double getDonutArea(double diameterInner, double diameterOuter){
        double areaInner = gen::getCircleArea(diameterInner);
        double areaOuter = gen::getCircleArea(diameterOuter);
        double donutArea = areaOuter - areaInner;

        return donutArea;
    }
}

#endif //TUDATBUNDLE_GENERAL_FUNCTIONS_H
