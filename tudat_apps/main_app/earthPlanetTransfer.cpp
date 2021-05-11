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

    //! Get path for output directory.
    static inline std::string getOutputPath2( const std::string& extraDirectory = "" )
    {
        // Declare file path string assigned to filePath.
        // __FILE__ only gives the absolute path of the header file!
        std::string filePath_( __FILE__ );

        // Strip filename from temporary string and return root-path string.
        std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                                       std::string( "applicationOutput.h" ).length( ) );
        std::string outputPath = reducedPath + "SimulationOutput/";
        if( extraDirectory != "" )
        {
            outputPath += extraDirectory;
        }

        if( outputPath.at( outputPath.size( ) - 1 ) != '/' )
        {
            outputPath += "/";
        }

        return outputPath;
    }

}

namespace matt_trajectories{

//    // Just the sun's gravitational parameter
//    const double sunGravitationalParameter = 1.32712428e20;

    // Function to return number of legs in a trajectory
    int getNumberOfLegs(std::vector< int > flybySequence){
        return flybySequence.size();
    }

    // Function to get the leg type vector - says whether each leg is departure, intermediate, or end. Uses number of legs to do it
    std::vector< TransferLegType > getLegTypeVector(int numberOfLegs){

        std::vector< TransferLegType > legTypeVector;
        legTypeVector.resize( numberOfLegs );
        legTypeVector[ 0 ] = mga_Departure;
        legTypeVector[ numberOfLegs - 1 ] = capture;

        // For loop sets intermediate legs (if any) to swingbys
        for(int i = 1; i < numberOfLegs - 1; i++){
            legTypeVector[ i ] = mga_Swingby;
        }

        return legTypeVector;
    }

    // Function to return the target planet index - used for the type of semi-major axis / eccentricity target
    int getTargetPlanetIndex(std::vector< int > flybySequence, int numberOfLegs){
        return flybySequence[numberOfLegs - 1];
    }

    // Function to modify ephemeris vector, gravitation parameter vector, and minimum pericenter radii vectors in place using pointers
    void set_Ephemerides_GravitationParameters_MinPeRadii(std::vector< ephemerides::EphemerisPointer >& ephemerisVector,
                                                           Eigen::VectorXd& gravitationalParameterVector,
                                                           Eigen::VectorXd& minimumPericenterRadii,
                                                           std::vector< int > flybySequence,
                                                           int numberOfLegs){
        // Create the ephemeris, gravitational parameter, and minimum pericentre vector.
        ephemerisVector.resize( numberOfLegs );
        gravitationalParameterVector.resize( numberOfLegs );
        minimumPericenterRadii.resize( numberOfLegs );

        // For loop to set gravitational parameters and minimum pericenter radii for each planet. Uses default values
        for(int i = 0; i < numberOfLegs; i++)
        {
            switch(flybySequence[ i ])
            {
                case( 1 ):
                    ephemerisVector[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                            ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury );
                    gravitationalParameterVector[ i ] = 2.2032E13;
                    minimumPericenterRadii[ i ] = 2639.7E3;
                    break;
                case( 2 ):
                    ephemerisVector[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                            ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
                    gravitationalParameterVector[ i ] = 3.24859E14;
                    minimumPericenterRadii[ i ] = 6251.8E3;
                    break;
                case( 3 ):
                    ephemerisVector[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                            ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
                    gravitationalParameterVector[ i ] = 3.986004418E14;
                    minimumPericenterRadii[ i ] = 6578.1E3;
                    break;
                case( 4 ):
                    ephemerisVector[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                            ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );
                    gravitationalParameterVector[ i ] = 4.282837E13;
                    minimumPericenterRadii[ i ] = 3596.2E3;
                    break;
                case( 5 ):
                    ephemerisVector[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                            ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter );
                    gravitationalParameterVector[ i ] = 1.26686534E17;
                    minimumPericenterRadii[ i ] = 72000.0E3;
                    break;
                case( 6 ):
                    ephemerisVector[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                            ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn );
                    gravitationalParameterVector[ i ] = 3.7931187E16;
                    minimumPericenterRadii[ i ] = 61000.0E3;
                    break;
                case( 7 ):
                    ephemerisVector[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                            ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::uranus );
                    gravitationalParameterVector[ i ] = 5.793939E15;
                    minimumPericenterRadii[ i ] = 26000.0E3;
                    break;
                case( 8 ):
                    ephemerisVector[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                            ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::neptune );
                    gravitationalParameterVector[ i ] = 6.836529E15;
                    minimumPericenterRadii[ i ] = 25000.0E3;
                    break;
                case( 9 ):
                    ephemerisVector[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                            ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::pluto );
                    gravitationalParameterVector[ i ] = 8.71E11;
                    minimumPericenterRadii[ i ] = 1395.0E3;
                    break;
                default:
                    std::cerr<<"Planet in flyby sequence is not defined.";
            }
        }
    }

    // Function to set semiMajorAxes and eccentricities vectors by modifying in-place
    void set_SMAs_Eccentricities(Eigen::VectorXd& semiMajorAxes,
                                 Eigen::VectorXd& eccentricities,
                                 int targetPlanetIndex){

        // Does nominal case for all values except mars, which has a parabolic arrival orbit
        if (targetPlanetIndex == 4){
            // Create Mars departure and capture variables. Parabolic-parabolic
            semiMajorAxes.resize( 2 );
            eccentricities.resize( 2 );
            semiMajorAxes << std::numeric_limits< double >::infinity( ), std::numeric_limits< double >::infinity( );
            eccentricities << 0., 0.;
        }
            // Nominal values
        else{
            // Create departure and capture variables. Uses default values - parabolic escape and (basically) parabolic entry
            semiMajorAxes.resize( 2 );
            eccentricities.resize( 2 );
            semiMajorAxes << std::numeric_limits< double >::infinity( ), 1.0895e8 / 0.02;
            eccentricities << 0., 0.98;
        }

    }



};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

EarthPlanetTransfer::EarthPlanetTransfer(std::vector< std::vector< double > > &bounds,
                                         std::vector< double > &deltaVBounds,
                                         std::vector< int > flybySequence,
                                         bool normaliseValues,
                                         bool includeDepartureDV,
                                         bool includeArrivalDV,
                                         std::string outputSubFolder):
                                         problemBounds_(bounds),
                                         deltaVBounds_(deltaVBounds),
                                         normaliseValues_(normaliseValues),
                                         includeDepartureDV_(includeDepartureDV),
                                         includeArrivalDV_(includeArrivalDV),
                                         outputSubFolder_(outputSubFolder){

    // Specify number of legs and type of legs
    numberOfLegs_ = matt_trajectories::getNumberOfLegs(flybySequence);
    legTypeVector_ = matt_trajectories::getLegTypeVector(numberOfLegs_);
    int targetPlanetIndex = matt_trajectories::getTargetPlanetIndex(flybySequence, numberOfLegs_);

    // Set ephemeris vector, gravitational parameter vector, and minimum pericenter radii vector
    matt_trajectories::set_Ephemerides_GravitationParameters_MinPeRadii(ephemerisVector_,
                                                                        gravitationalParameterVector_,
                                                                        minimumPericenterRadii_,
                                                                        flybySequence,
                                                                        numberOfLegs_);



    // Set semiMajorAxes and eccentircities vectors
    matt_trajectories::set_SMAs_Eccentricities(semiMajorAxes_, eccentricities_, targetPlanetIndex);
    
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
    const double sunGravitationalParameter = matt_trajectories::sunGravitationalParameter;


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
                                includeDepartureDV_,
                                includeArrivalDV_);

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


    /////////////// TEST SECTION for saving trajectory to file TODO: REMOVE ME IF NOT NEEDED /////////////////////////

//    std::cout << "Saving le trajectory " << std::endl;
//    // Define vectors to calculate intermediate points
//    std::vector< Eigen::Vector3d > intermediatePositionVector;
//    std::vector< double > intermediateTimeVector;
//    earthPlanetTraj.intermediatePoints(1000, intermediatePositionVector, intermediateTimeVector);
//////    std::cout << "One of the position vector values: " << intermediatePositionVector[0] << std::endl;
////    std::string saveFilename = "/home/matt/tudatBundle/tudatApplications/EDT_thesis/tudat_apps/main_app/THIS.dat";
////
//////    std::string saveFilename = "/this.dat";
//    std::string outputPathBase = "/home/matt/tudatBundle/tudatApplications/EDT_thesis/tudat_apps/main_app/SimulationOutput/" + outputSubFolder_;
//    std::string saveFilename = outputPathBase + "/THIS.dat";
//    std::cout << "output subfolder: " << outputSubFolder_ << std::endl;
//    std::cout << "save file name: " << saveFilename << std::endl;
//    tudat::transfer_trajectories::writeTrajectoryToFile( intermediatePositionVector, intermediateTimeVector, saveFilename );

    /////////////////////// TEST SECITON OVER ///////////////////////////////////////////////////////////////


    // If statement for including normalised values
    if (normaliseValues_) {
        return {normalisedDV, normalisedTOF};
    }
    else{
        return {resultingDeltaV, TOF};
    }



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

