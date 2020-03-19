//
// Created by matt on 29/02/2020.
//

#ifndef TUDATBUNDLE_UNIVERSAL_SETTINGS_H
#define TUDATBUNDLE_UNIVERSAL_SETTINGS_H


#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h>

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace EDTs;

namespace univ {

    // Define Class that contains celestial bodies, and creates acceleration map etc from vehicle body map
    class propBodies {
    public:
        // Make relevant parameters accessible to extern
        std::vector<std::string> bodiesToCreate;
        std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;
        NamedBodyMap bodyMap;
        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        basic_astrodynamics::AccelerationMap accelerationModelMap;

        // Constructor
        propBodies(
                EDTs::EDTConfig vehicleConfig) :
                vehicleConfig_(vehicleConfig) {
            // Constructor function

            /////////// Define body settings for simulation and create bodymap./////////////////////////
            bodiesToCreate.push_back("Sun");
            bodiesToCreate.push_back("Earth");
            bodiesToCreate.push_back("Jupiter");
            bodiesToCreate.push_back("Saturn");
//            bodiesToCreate.push_back("Uranus"); TODO: encountered SPICE error -- fix this!
//            bodiesToCreate.push_back("Neptune");

            // Create body objects.
            bodySettings =
                    getDefaultBodySettings( bodiesToCreate );

            bodyMap = createBodies( bodySettings );

            // Add vehicle to bodyMap
            bodyMap["Vehicle"] = vehicleConfig.EDTBody;

            // Finalise body creation
            setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

            ///////////////// Create bodies acceleration map///////////////////////


            // Define acceleration model settings for vehicle
            accelerationsOfVehicle[ "Vehicle" ].push_back(
                    std::make_shared< ThrustAccelerationSettings >( vehicleConfig.thrustDirectionGuidanceSettings, vehicleConfig.thrustMagnitudeSettings ) );

            // Define acceleration model settings for celestials using for loop
            for(std::size_t i=0; i<bodiesToCreate.size(); ++i){
                accelerationsOfVehicle[bodiesToCreate[i]].push_back(std::make_shared< AccelerationSettings >( central_gravity ) );
            }

            // Create acceleration map
            accelerationMap["Vehicle"] = accelerationsOfVehicle;

            // Add bodies to bodies to propagate and central bodies
            bodiesToPropagate.push_back("Vehicle");
            centralBodies.push_back("Sun");

            // Create acceleration models and propagation settings
            accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


        }

    private:
        // Set values of inputs as private
        EDTs::EDTConfig vehicleConfig_;
    };



    // Define Class that creates propagation and integration settings
    class propSettings {
    public:
        // Make relevant parameters accessible to extern
        std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
        std::shared_ptr< PropagationTerminationSettings > terminationSettings;
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings;
        std::shared_ptr< IntegratorSettings< > > integratorSettings;

        // Constructor
        propSettings(
                propBodies simulationPropBodies,
                Eigen::Vector6d vehicleInitialState,
                double initialEphemerisTime = 1.0E7, //TODO: set default ephem time to for eg 2025
                bool useDefaultInitState=true,
                bool useDefaultTimeTermination=true) :
                vehicleInitialState_(vehicleInitialState) {
            // Constructor function

            /////////////// Propagation Settings ///////////////////////

            // Set initial state according to defaults or the ones set by class creator
            if(useDefaultInitState){
                vehicleInitialState = Eigen::Vector6d::Zero( ); //TODO: Make sure current set values (Earth orbit around Sun) are correct
                vehicleInitialState[0] = 1.496e11;
                vehicleInitialState[1] = 0.0;
                vehicleInitialState[2] = 0.0;
                vehicleInitialState[3] = 0.0;
                vehicleInitialState[4] = 29.78E3;
                vehicleInitialState[5] = 0.0;
            }

            // Define propagation termination conditions (varies, so use hybrid termination settings, see lunarascentgroup for an implementation)
            if (useDefaultTimeTermination){
                // Add default time termination to termination settings list
                double defaultFinalEphemTime = initialEphemerisTime + 10*365*24*60*60; //TODO: Make this a reasonable number
                terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                        defaultFinalEphemTime ) );
            }
            else{
                //TODO: work out what happens here
            }

            // Using above if-else really make the termination settings
            terminationSettings = std::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );

            // Define settings for propagation of translational dynamics.
            propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >(
                            simulationPropBodies.centralBodies, simulationPropBodies.accelerationModelMap,
                            simulationPropBodies.bodiesToPropagate, vehicleInitialState, terminationSettings );


            ///////////////////// Integration Settings //////////////////////
            // TODO: Check general settings of integrator, and adjust as necessary (especially tolerances)
            integratorSettings =
                    std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                            initialEphemerisTime,
                            1.0,
                            numerical_integrators::RungeKuttaCoefficients::rungeKutta87DormandPrince, //TODO: change this integrator to FILG if possible (ask dominic?)
                            1E-5,
                            365*24*60*60,
                            1e-14,
                            1e-14);
        }

    private:
        // Set values of inputs as private
        Eigen::Vector6d vehicleInitialState_;
    };




//// getBodies: create NamedBodyMap of standard celestial bodies to simulate
//    NamedBodyMap getCBodyMap() {
//
//        // Define body settings for simulation.
//        std::vector<std::string> bodiesToCreate;
//        bodiesToCreate.push_back("Sun");
//        bodiesToCreate.push_back("Earth");
//        bodiesToCreate.push_back("Jupiter");
//        bodiesToCreate.push_back("Saturn");
//        bodiesToCreate.push_back("Uranus");
//        bodiesToCreate.push_back("Neptune");
//
//        // Create body objects.
//        std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
//                getDefaultBodySettings( bodiesToCreate );
//
//        NamedBodyMap bodyMap = createBodies( bodySettings );
//
//        return bodyMap;
//    }
//
//// Define acceleration model settings based on bodyMap
//    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > getAccelerationSettings() {
//
//    }

    void finBodyCreation(NamedBodyMap bodyMap, std::string frameOrigin = "SSB", std::string frameOrientation = "ECLIPJ2000"){
        setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );
    }


    void printIntVector(std::vector<int> const &input) {
        for (int i = 0; i < input.size(); i++) {
            std::cout << input.at(i) << ' ';
        }
    }

    void printStrVector(std::vector<std::string> const &input) {
        for (int i = 0; i < input.size(); i++) {
            std::cout << input.at(i) << ' ';
        }
    }

}

#endif //TUDATBUNDLE_UNIVERSAL_SETTINGS_H
