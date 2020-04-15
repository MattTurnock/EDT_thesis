// Main file for SSO-CHB simulation run

///////// Start by defining the SSO simulation ///////////////

#include "simulation_SSO-CHB.h"
#include "EDT_configs.h"
#include "EDTGuidance.h"
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

    double pi = 3.14159;
    double deg2rad = pi / 180;
    double AU = 1.496e11;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             SIMULATION PREP            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout<< "===============Prepping Sim==================" << std::endl;
    std::cout<< "Loading spice kernels" << std::endl;
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create base body map to be built on by other classes + give initial numbers for magfield data
    NamedBodyMap baseBodyMap;

    double B0 = 12E-9;
    double phi0 = 55 * deg2rad;
    double R0 = 1*AU;

    // Create EDT Environment class
    std::cout<< "Creating environment class" << std::endl;
    EDTEnvironment CHBEDTEnviro = EDTEnvironment(B0, phi0, R0, baseBodyMap);

    // Create EDT Guidance class
    std::cout<< "Creating Guidance class" << std::endl;
    EDTGuidance CHBEDTGuidance = EDTGuidance(
            "nominal",
            "nominalPrograde",
            baseBodyMap,
            CHBEDTEnviro);

    // Create EDT config class and set constant thrust in guidance class
    std::cout<< "Creating Config class" << std::endl;
    EDTs::EDTConfig CHBEDTConfig = EDTs::EDTConfig(CHBEDTGuidance);
    CHBEDTGuidance.setThrustMagnitudeConstant(CHBEDTConfig.getConstantThrust());

    // Get universal class for propagation bodies
    std::cout<< "Creating Propbodies class" << std::endl;
    univ::propBodies SSOPropBodies = univ::propBodies(CHBEDTConfig, baseBodyMap);

    // Get universal class for propagation settings + set vehicle initial state
    std::cout<< "Creating Propsettings class" << std::endl;
    univ::propSettings SSOPropSettings = univ::propSettings(SSOPropBodies);

    // Ensure environment is properly updated
    CHBEDTEnviro.updateAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "====================Propagating==========================" << std::endl;

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
            SSOPropBodies.bodyMap, SSOPropSettings.integratorSettings, SSOPropSettings.propagatorSettings);
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    // Write satellite propagation history to file.
    std::string outputSubFolder = "SSO-CHB-Test-custom/";
    input_output::writeDataMapToTextFile( integrationResult,
                                          "SSO-CH-Test-out-expo-enviro__E-6.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

}