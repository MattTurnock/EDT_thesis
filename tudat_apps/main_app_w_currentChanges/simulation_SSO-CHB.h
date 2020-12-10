//
// Created by matt on 29/02/2020.
// Header file to include a bunch of only necessary things for the simulation_SSO-CHB.cpp file.
//

#ifndef TUDATBUNDLE_SIMULATION_SSO_CHB_H
#define TUDATBUNDLE_SIMULATION_SSO_CHB_H

/// COPY OF tudatSimulationHeader.h ///

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateRepresentationConversions.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/mapTextFileReader.h"

#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/Filters/createFilter.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Mathematics/Interpolators/createInterpolator.h"

#ifdef USE_CSPICE
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#endif

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"
#include "Tudat/SimulationSetup/EstimationSetup/estimatableParameterSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTerminationSettings.h"
#include "Tudat/SimulationSetup/EstimationSetup/createNumericalSimulator.h"
#include "Tudat/SimulationSetup/PropagationSetup/thrustSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/createMassRateModels.h"

#include "Tudat/JsonInterface/Support/errorHandling.h"

// Add pagmo files
#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/simulated_annealing.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/nsga2.hpp"

#include <fstream>

#endif //TUDATBUNDLE_SIMULATION_SSO_CHB_H
