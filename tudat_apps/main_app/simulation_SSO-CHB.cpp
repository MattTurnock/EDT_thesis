// Main file for SSO-CHB simulation run

///////// Start by defining the SSO simulation ///////////////

#include "simulation_SSO-CHB.h"
#include "EDT_configs.h"
#include "universal_settings.h"
#include "applicationOutput.h"
#include "environment_settings.h"


int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::interpolators;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_mathematics;
    using namespace tudat::basic_astrodynamics;

    using namespace univ;
    using namespace EDTs;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             SIMULATION PREP            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout<< "Prepping Sim" << std::endl;
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create EDT class
    EDTs::EDTConfig CHBEDTConfig = EDTs::EDTConfig(12.2, false); //TODO: Change from testmode when ready

    // Get universal class for propagation bodies
    univ::propBodies SSOPropBodies = univ::propBodies(CHBEDTConfig);

    // Get universal class for propagation settings + set vehicle initial state
    Eigen::Vector6d vehicleInitialState = Eigen::Vector6d::Zero( ); // NOTE: this is not zeros if default used
    univ::propSettings SSOPropSettings = univ::propSettings(SSOPropBodies, vehicleInitialState);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "Propagating" << std::endl;
    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
            SSOPropBodies.bodyMap, SSOPropSettings.integratorSettings, SSOPropSettings.propagatorSettings);

    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    std::string outputSubFolder = "SSO-CHB-Test-custom/";

    // Write satellite propagation history to file.
    input_output::writeDataMapToTextFile( integrationResult,
                                          "SSO-CH-Test-out-expo-0-000325.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );


    std::cout<< "===========================\nPerforming interpolator test" << std::endl; // TODO: remove this and below sections as necessary

    // Load data.
    std::map< double, Eigen::Vector6d > stateMap;
//    stateMap = ....

    // empty map container
    std::map<double, double> gquiz1;

    // insert elements in random order
    gquiz1.insert(std::pair<int, int>(1, 40));
//    gquiz1.insert(std::pair<int, int>(2, 30));
//    gquiz1.insert(std::pair<int, int>(3, 60));
//    gquiz1.insert(std::pair<int, int>(4, 20));
//    gquiz1.insert(std::pair<int, int>(5, 50));
//    gquiz1.insert(std::pair<int, int>(6, 50));
//    gquiz1.insert(std::pair<int, int>(7, 10));


// Create interpolator
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            std::make_shared< interpolators::InterpolatorSettings >( linear_interpolator );
    std::shared_ptr< OneDimensionalInterpolator< double, double > > interpolator =
            interpolators::createOneDimensionalInterpolator(
                    gquiz1, interpolatorSettings );

// Interpolate
    double interpolatedResult = interpolator->interpolate( 1.2 );

    std::cout << "Interpolated result for 1.2: " << std::to_string(interpolatedResult) << std::endl;


    std::cout<< "Done!" << std::endl;

}