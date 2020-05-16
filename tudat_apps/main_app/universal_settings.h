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

//    void IDFunction(std::string bodyName){
//        int bodyID = spice_interface::convertBodyNameToNaifId(bodyName);
//        std::cout << "Body ID for " << bodyName << ": " << std::to_string(bodyID) << std::endl;
//    } // TODO: Delete me



    // Define Class that contains celestial bodies, and creates acceleration map etc from vehicle body map
    class propBodies {
    public:
        // Public simulation variables
        NamedBodyMap& bodyMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;
        basic_astrodynamics::AccelerationMap accelerationModelMap;

        // Constructor
        propBodies(
                EDTs::EDTConfig& vehicleConfig,
                NamedBodyMap& baseBodyMap) :
                vehicleConfig_(vehicleConfig),
                bodyMap(baseBodyMap)
                {
            // Constructor

            /////////// Define body settings for simulation and create bodymap./////////////////////////
            bodiesToCreate_.push_back("Sun"); // TODO: uncomment following bodies:
//            bodiesToCreate_.push_back("Earth");
//            bodiesToCreate_.push_back("Jupiter");
//            bodiesToCreate_.push_back("Saturn");
//            bodiesToCreate_.push_back("Uranus");
//            bodiesToCreate_.push_back("Neptune");

            // Create body objects.
            bodySettings_ = getDefaultBodySettings( bodiesToCreate_ );
            bodyMap = createBodies( bodySettings_ );

            // Add vehicle to bodyMap
            bodyMap["Vehicle"] = vehicleConfig.EDTBody;

            // Finalise body creation
            setGlobalFrameBodyEphemerides( bodyMap, "Sun", "ECLIPJ2000" );

            ///////////////// Create bodies acceleration map (includes perturbations) ///////////////////////

            // Add acceleration model settings for vehicle - ie thrust guidance settings
            accelerationsOfVehicle_[ "Vehicle" ].push_back(
                    std::make_shared< ThrustAccelerationSettings >( vehicleConfig.thrustDirectionGuidanceSettings, vehicleConfig.thrustMagnitudeSettings ) );

            // Add acceleration model settings for celestials using for loop
            for(std::size_t i=0; i<bodiesToCreate_.size(); ++i){
                accelerationsOfVehicle_[bodiesToCreate_[i]].push_back(std::make_shared< AccelerationSettings >( central_gravity ) );
            }

            // Add acceleration model for SRP TODO: make sure "Vehicle" sizing includes a way for canonball pressure to work, https://tudat.tudelft.nl/tutorials/applicationWalkthroughs/perturbedEarthOrbitingSatellite.html?highlight=acceleration
            // TODO: also consider using another pressure model, such as panelled (see link) https://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/availableSettings.html#radiationpressuremodeloptions
//            accelerationsOfVehicle_["Sun"].push_back(std::make_shared< AccelerationSettings >(basic_astrodynamics::cannon_ball_radiation_pressure));

            // Create acceleration map
            accelerationMap_["Vehicle"] = accelerationsOfVehicle_;

            // Add bodies to bodies to propagate and central bodies
            bodiesToPropagate.push_back("Vehicle");
            centralBodies.push_back("Sun");

            // Create acceleration models and propagation settings
            accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap_, bodiesToPropagate, centralBodies );

        }

        NamedBodyMap getBodyMap(){
            return bodyMap;
        }

        EDTs::EDTConfig getVehicleConfig(){
            return vehicleConfig_;
        }

    protected:

        /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////
        EDTs::EDTConfig& vehicleConfig_;


        /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////

        /////////////////////////// Other set parameters ///////////////////////////////
        // Body variables
        std::vector<std::string> bodiesToCreate_;
        std::map< std::string, std::shared_ptr< BodySettings > > bodySettings_;

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap_;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle_;


    };



    // Define Class that creates propagation and integration settings
    class propSettings {
    public:
        // Make relevant parameters accessible to extern
        std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
        std::shared_ptr< PropagationTerminationSettings > terminationSettings;
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings;
        std::shared_ptr< IntegratorSettings< > > integratorSettings;

        // Define list of dependent variables to save.
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;

        // Constructor
        propSettings(
                propBodies& simulationPropBodies,
                Eigen::Vector3d vehicleInitialPosition = {1.496E11, 0, 0},
                Eigen::Vector3d vehicleInitialVelocity = {0, 29.78E3, 0},
                double initialEphemerisTime = 1.0E7,
                std::string terminationMethod = "nominalTimeTermination",
                double simulationLength = 10) :
                simulationPropBodies_(simulationPropBodies),
                vehicleInitialPosition_(vehicleInitialPosition),
                vehicleInitialVelocity_(vehicleInitialVelocity),
                initialEphemerisTime_(initialEphemerisTime),
                terminationMethod_(terminationMethod),
                simulationLength_(simulationLength){

            /////////////// Propagation Settings ///////////////////////
            // Create vehicle initial state from position and velocities
            vehicleInitialState_ << vehicleInitialPosition_, vehicleInitialVelocity_;


            // Define propagation termination conditions
            if (terminationMethod_ == "nominalTimeTermination"){
                // Add default time termination to termination settings list
                double defaultFinalEphemTime = initialEphemerisTime_ + simulationLength_*365*24*60*60;
                terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                        defaultFinalEphemTime ) );
            }
            else{
                //TODO: work out what happens here
            }

            // Using above if-else really make the termination settings
            terminationSettings = std::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );

            // Add dependent variables to save to list (kepler + altitude)
            dependentVariablesList.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >( keplerian_state_dependent_variable,
                            "Vehicle" ,
                            "Sun" ) );

            dependentVariablesList.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(relative_distance_dependent_variable,
                            "Vehicle", "Sun"));


            dependentVariablesToSave_ =
                    std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

            // Define settings for propagation of translational dynamics.
            propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >(
                            simulationPropBodies_.centralBodies, simulationPropBodies_.accelerationModelMap,
                            simulationPropBodies_.bodiesToPropagate, vehicleInitialState_, terminationSettings,
                            cowell, dependentVariablesToSave_ );


            ///////////////////// Integration Settings //////////////////////
            // TODO: Check general settings of integrator, and adjust as necessary (especially tolerances)
            integratorSettings =
                    std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                            initialEphemerisTime_,
                            1.0,
                            numerical_integrators::RungeKuttaCoefficients::rungeKutta87DormandPrince, //TODO: change this integrator to FILG if possible (ask dominic?)
                            1E-5,
                            365*24*60*60,
                            1e-14,
                            1e-14);
        }

        propBodies getSimulationPropBodies(){
            return simulationPropBodies_;
        }

    protected:

        /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////
        propBodies& simulationPropBodies_;
        Eigen::Vector3d vehicleInitialPosition_;
        Eigen::Vector3d vehicleInitialVelocity_;


        /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////
        double initialEphemerisTime_;
        std::string terminationMethod_;
        double simulationLength_;


        /////////////////////////// Other set parameters ///////////////////////////////
        Eigen::Vector6d vehicleInitialState_;
        std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave_;
    };


//    void printIntVector(std::vector<int> const &input) {
//        for (int i = 0; i < input.size(); i++) {
//            std::cout << input.at(i) << ' ';
//        }
//    }
//
//    void printStrVector(std::vector<std::string> const &input) {
//        for (int i = 0; i < input.size(); i++) {
//            std::cout << input.at(i) << ' ';
//        }
//    }

}

#endif //TUDATBUNDLE_UNIVERSAL_SETTINGS_H
