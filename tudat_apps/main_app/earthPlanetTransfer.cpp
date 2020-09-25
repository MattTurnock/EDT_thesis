//
// Created by matt on 08/08/2020.
//

#include "earthPlanetTransfer.h"

using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::transfer_trajectories;
using namespace tudat;
using namespace pagmo;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TODO: SOURCE ME FROM GENERAL FUNCTIONS - had problems with importing, needed to test quick
namespace gen2 {
// Function to normalise or denormlaise deltaV or TOF, based on given bounds
    double normaliseValue(double initialValue, double lowerBound, double upperBound, bool denormalise = false) {
        double newValue;
        // if statement true, then denormalised
        if (denormalise) {
            double denormalisedValue = initialValue * (upperBound - lowerBound) + lowerBound;
            newValue = denormalisedValue;
        }
            // otherwise normalised
        else {
            double normalisedValue = (initialValue - lowerBound) / (upperBound - lowerBound);
            newValue = normalisedValue;
        }

        return newValue;
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

EarthPlanetTransfer::EarthPlanetTransfer(std::vector< std::vector< double > > &bounds,
                                         std::vector< double > &deltaVBounds,
                                         std::vector< int > flybySequence):
                                         problemBounds_(bounds),
                                         deltaVBounds_(deltaVBounds){

    // Specify number of legs and type of legs
    numberOfLegs_ = flybySequence.size();
    legTypeVector_.resize( numberOfLegs_ );
    legTypeVector_[ 0 ] = mga_Departure;
    legTypeVector_[ numberOfLegs_ - 1 ] = capture;

    // For loop sets intermediate legs (if any) to swingbys
    for(int i = 1; i < numberOfLegs_ - 1; i++){
        legTypeVector_[ i ] = mga_Swingby;
    }

    // Create the ephemeris, gravitational parameter, and minimum pericentre vector.
    ephemerisVector_.resize( numberOfLegs_ );
    gravitationalParameterVector_.resize( numberOfLegs_ );
    minimumPericenterRadii_.resize( numberOfLegs_ );

    // For loop to set gravitational parameters and minimum pericenter radii for each planet. Uses default values
    for(int i = 0; i < numberOfLegs_; i++)
    {
        switch(flybySequence[ i ])
        {
            case( 1 ):
                ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                        ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury );
                gravitationalParameterVector_[ i ] = 2.2032E13;
                minimumPericenterRadii_[ i ] = 2639.7E3;
                break;
            case( 2 ):
                ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                        ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
                gravitationalParameterVector_[ i ] = 3.24859E14;
                minimumPericenterRadii_[ i ] = 6251.8E3;
                break;
            case( 3 ):
                ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                        ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
                gravitationalParameterVector_[ i ] = 3.986004418E14;
                minimumPericenterRadii_[ i ] = 6578.1E3;
                break;
            case( 4 ):
                ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                        ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );
                gravitationalParameterVector_[ i ] = 4.282837E13;
                minimumPericenterRadii_[ i ] = 3596.2E3;
                break;
            case( 5 ):
                ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                        ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter );
                gravitationalParameterVector_[ i ] = 1.26686534E17;
                minimumPericenterRadii_[ i ] = 72000.0E3;
                break;
            case( 6 ):
                ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                        ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn );
                gravitationalParameterVector_[ i ] = 3.7931187E16;
                minimumPericenterRadii_[ i ] = 61000.0E3;
                break;
            case( 7 ):
                ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                        ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::uranus );
                gravitationalParameterVector_[ i ] = 5.793939E15;
                minimumPericenterRadii_[ i ] = 26000.0E3;
                break;
            case( 8 ):
                ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                        ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::neptune );
                gravitationalParameterVector_[ i ] = 6.836529E15;
                minimumPericenterRadii_[ i ] = 25000.0E3;
                break;
            case( 9 ):
                ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                        ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::pluto );
                gravitationalParameterVector_[ i ] = 8.71E11;
                minimumPericenterRadii_[ i ] = 1395.0E3;
                break;
            default:
                std::cerr<<"Planet in flyby sequence is not defined.";
        }
    }

    // Create departure and capture variables. Uses default values - parabolic escape and (basically) parabolic entry
    semiMajorAxes_.resize( 2 );
    eccentricities_.resize( 2 );
    semiMajorAxes_ << std::numeric_limits< double >::infinity( ), 1.0895e8 / 0.02;
    eccentricities_ << 0., 0.98;
}

//! Descriptive name of the problem
std::string EarthPlanetTransfer::get_name() const {
    return "Earth-Jupiter Earth-Saturn trajectories (supports MGA)";
}

//! Get bounds
std::pair<std::vector<double>, std::vector<double> > EarthPlanetTransfer::get_bounds() const {

    return { problemBounds_[0], problemBounds_[1] };
}



//! Implementation of the fitness function. Returns TOF
std::vector<double> EarthPlanetTransfer::fitness( const std::vector<double> &xv ) const {
    // Sun gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;


    // Create variable vector.
    Eigen::VectorXd variableVector ( numberOfLegs_ + 1 );

    double TOF = 0;
    for(int i = 0; i < numberOfLegs_ ; i++){
        variableVector[ i ] = xv[ i ];
        if( i > 0 ){
            TOF += xv[i];
        }
    }
    variableVector[ numberOfLegs_ ] = 1;//dummy
    variableVector *= physical_constants::JULIAN_YEAR; // TOF now given in years, so is required

    // Create the trajectory problem
    Trajectory earthPlanetTraj( numberOfLegs_,
                                legTypeVector_,
                                ephemerisVector_,
                                gravitationalParameterVector_,
                                variableVector,
                                sunGravitationalParameter,
                                minimumPericenterRadii_,
                                semiMajorAxes_,
                                eccentricities_,
                                true,
                                false); // TODO: make sure arrival DV is not included!! ================ <---

    // Create the DV vector and calculate trajectory
    double resultingDeltaV;
    earthPlanetTraj.calculateTrajectory(resultingDeltaV);
    resultingDeltaV = resultingDeltaV / 1000; // TODO: check if this is right
    // Give large TOF penalty if DV exceeds constrained value, or cannot be calculated TODO: check if needs to be readded
    if ( (std::isnan(resultingDeltaV)) or ( resultingDeltaV < deltaVBounds_[0] ) or (resultingDeltaV > deltaVBounds_[1]) ){
        TOF += 1.0E12;
    }

    double normalisedDV = gen2::normaliseValue(resultingDeltaV, deltaVBounds_[0], deltaVBounds_[1], false);
    double normalisedTOF = gen2::normaliseValue(TOF, problemBounds_[0][1], problemBounds_[1][1], false);

    return {normalisedDV, normalisedTOF};
//    return {resultingDeltaV, TOF};
//    return {resultingDeltaV};
}



// Some additional functions for returning important info (see header file)
std::vector< TransferLegType > EarthPlanetTransfer::getLegTypeVector() const {
    return legTypeVector_;
}
std::vector< std::string > EarthPlanetTransfer::getBodyNamesVector() const {
    return bodyNamesVector_;
}
std::vector< ephemerides::EphemerisPointer > EarthPlanetTransfer::getEphemerisVector() const{
    return ephemerisVector_;
}
Eigen::VectorXd EarthPlanetTransfer::getGravitationalParameterVector() const{
    return gravitationalParameterVector_;
}
Eigen::VectorXd EarthPlanetTransfer::getSemiMajorAxes() const{
    return semiMajorAxes_;
}
Eigen::VectorXd EarthPlanetTransfer::getEccentricities() const{
    return eccentricities_;
}
Eigen::VectorXd EarthPlanetTransfer::getMinimumPericenterRadii() const{
    return minimumPericenterRadii_;
}