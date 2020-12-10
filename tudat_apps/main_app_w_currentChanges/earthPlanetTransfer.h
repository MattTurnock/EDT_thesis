//
// Created by matt on 08/08/2020.
//

#ifndef TUDATBUNDLE_EARTH_PLANET_TRANSFER
#define TUDATBUNDLE_EARTH_PLANET_TRANSFER

#include <vector>
#include <utility>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Basics/testMacros.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>

#include <Tudat/InputOutput/basicInputOutput.h>

#include "Tudat/SimulationSetup/PropagationSetup/propagationPatchedConicFullProblem.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationLambertTargeterFullProblem.h"

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"
#include <random>

#include "pagmo/island.hpp"
#include "pagmo/io.hpp"
#include "pagmo/problem.hpp"
#include <pagmo/rng.hpp>

#include <Eigen/Core>

using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::transfer_trajectories;
using namespace tudat;
using namespace pagmo;

/*!
 *  The class defined in this file is to be used in a Pagmo optimization. It defines the objective function for an
 *  Earth-Jupiter or Earth-Saturn two-burn impulsive transfer (although other flyby sequences are supported), using a
 *  Lambert targeter. The independent variables are:
 *
 *  1) Departure Julian day
 *  2) Time-of-flight of Earth-Mars transfer
 *
 *  The problem minimized the TOF, and bounds the problem to an upper DV
 */

//! Function almost identical to MultipleGravityAssist, adjusted to alter fitness etc
struct EarthPlanetTransfer
{
    // Empty constructor for pagmo compatibility
    EarthPlanetTransfer(){}

    EarthPlanetTransfer( std::vector< std::vector< double > > &bounds,
                           std::vector< double > &deltaVBounds,
                           std::vector< int > flybySequence,
                           bool normaliseValues,
                           bool includeDeaprtureDV,
                           bool includeArrivalDV);

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const;

    std::string get_name( ) const;

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    // Sets size of fitness vector (one for TOF, one for DV. NOTE: makes it multi-objective)
    vector_double::size_type get_nobj() const
    {
        return 2u;
    }

    // A bunch of functions to return private variables if required
    std::vector< TransferLegType > getLegTypeVector() const;
    std::vector< std::string > getBodyNamesVector() const;
    std::vector< ephemerides::EphemerisPointer > getEphemerisVector() const;
    Eigen::VectorXd getGravitationalParameterVector() const;
    Eigen::VectorXd getSemiMajorAxes() const;
    Eigen::VectorXd getEccentricities() const;
    Eigen::VectorXd getMinimumPericenterRadii() const;

private:

    const std::vector< std::vector< double > > problemBounds_;

    bool useTripTime_;

    int numberOfLegs_;
    std::vector< TransferLegType > legTypeVector_;
    std::vector< std::string > bodyNamesVector_;
    std::vector< ephemerides::EphemerisPointer > ephemerisVector_;
    Eigen::VectorXd gravitationalParameterVector_;
    Eigen::VectorXd semiMajorAxes_;
    Eigen::VectorXd eccentricities_;
    Eigen::VectorXd minimumPericenterRadii_;

    // CUSTOM ADD-IN FOR DV BOUNDS (and normalisation bool)
    std::vector< double > deltaVBounds_;
    bool normaliseValues_;
    bool includeDepartureDV_;
    bool includeArrivalDV_;
};


#endif // TUDATBUNDLE_EARTH_PLANET_TRANSFER