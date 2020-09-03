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
        // Public simulation variables
        NamedBodyMap& bodyMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;
        basic_astrodynamics::AccelerationMap accelerationModelMap;

        // Constructor
        propBodies(
                EDTs::EDTConfig& vehicleConfig,
                EDTGuidance& guidanceClass,
                NamedBodyMap& baseBodyMap,
                nlohmann::json jsonBodiesToInclude,
                bool useThrust = true) :
                vehicleConfig_(vehicleConfig),
                guidanceClass_(guidanceClass),
                bodyMap(baseBodyMap),
                jsonBodiesToInclude_(jsonBodiesToInclude),
                useThrust_(useThrust)
                {

            // Constructor
            /////////// Initialise all bodies list with all possible bodies that can be used //////////
            allBodiesList_.push_back("Sun");
            allBodiesList_.push_back("Mercury");
            allBodiesList_.push_back("Venus");
            allBodiesList_.push_back("Earth");
            allBodiesList_.push_back("Earth");
            allBodiesList_.push_back("Mars");
            allBodiesList_.push_back("Jupiter");
            allBodiesList_.push_back("Saturn");
            allBodiesList_.push_back("Uranus");
            allBodiesList_.push_back("Neptune");


            /////////// Define body settings for simulation and create bodymap./////////////////////////
            // Creat bodies list from json data
            for(std::size_t i=0; i<allBodiesList_.size(); ++i){
                intPlaceholder_ = jsonBodiesToInclude_[allBodiesList_[i]];
                if (jsonBodiesToInclude_[allBodiesList_[i]] == 1){
                    bodiesToCreate_.push_back(allBodiesList_[i]);
                }
            }

            std::cout << "Using bodies: " << std::endl;
            for(std::size_t i=0; i<bodiesToCreate_.size(); ++i){
                std::cout << bodiesToCreate_[i] << std::endl;
            }

            // Create body objects.
            bodySettings_ = getDefaultBodySettings( bodiesToCreate_ );
            bodyMap = createBodies( bodySettings_ );

            // Add vehicle to bodyMap
            bodyMap["Vehicle"] = vehicleConfig.EDTBody;

            // Finalise body creation
            setGlobalFrameBodyEphemerides( bodyMap, "Sun", "ECLIPJ2000" );

            ///////////////// Create bodies acceleration map (includes perturbations) ///////////////////////
            // Add acceleration model settings for vehicle - ie thrust guidance settings. Uses boolean for creating acceleration map without vehicle thrust if wanted
            if (useThrust_) {
                accelerationsOfVehicle_["Vehicle"].push_back(
                        std::make_shared<ThrustAccelerationSettings>(guidanceClass.thrustDirectionGuidanceSettings,
                                                                     guidanceClass.thrustMagnitudeSettings));
            }
            // Add acceleration model settings for celestials using for loop
            for(std::size_t i=0; i<bodiesToCreate_.size(); ++i){
                accelerationsOfVehicle_[bodiesToCreate_[i]].push_back(std::make_shared< AccelerationSettings >( central_gravity ) );
            }

            // Add acceleration model for cannonball SRP
            SRPReferenceArea_ = vehicleConfig_.getEffectiveSRPArea();
            SRPCoefficient_ = vehicleConfig_.getEffectiveSRPCoefficient();

            vehicleRadiationPressureSettings_ = std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", SRPReferenceArea_, SRPCoefficient_);

            SRPInterface_ = createRadiationPressureInterface(
                    vehicleRadiationPressureSettings_, "Vehicle", bodyMap );

            bodyMap["Vehicle"]->setRadiationPressureInterface("Sun", SRPInterface_ );

            accelerationsOfVehicle_["Sun"].push_back(std::make_shared< AccelerationSettings >(basic_astrodynamics::cannon_ball_radiation_pressure));

            // Create acceleration map
            accelerationMap_["Vehicle"] = accelerationsOfVehicle_;

            // Add bodies to bodies to propagate and central bodies
            bodiesToPropagate.push_back("Vehicle");
            centralBodies.push_back("Sun");

            // Create acceleration models and propagation settings
            accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap_, bodiesToPropagate, centralBodies );

        }

        //////////////////////////////// Get-set functions //////////////////////////////////////////////////////
        NamedBodyMap getBodyMap(){
            return bodyMap;
        }

        EDTs::EDTConfig getVehicleConfig(){
            return vehicleConfig_;
        }

        std::shared_ptr< electro_magnetism::RadiationPressureInterface > getSRPInterface(){
            return SRPInterface_;
        }

        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > getAccelerationsOfVehicle(){
            return accelerationsOfVehicle_;
        }

    protected:

        /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////
        EDTs::EDTConfig& vehicleConfig_;
        EDTGuidance& guidanceClass_;
        nlohmann::json jsonBodiesToInclude_;

        /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////
        bool useThrust_;
        /////////////////////////// Other set parameters ///////////////////////////////
        // Body variables
        std::vector<std::string> bodiesToCreate_;
        std::map< std::string, std::shared_ptr< BodySettings > > bodySettings_;
        std::vector<std::string> allBodiesList_;

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap_;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle_;
        double SRPReferenceArea_;
        double SRPCoefficient_;
        std::shared_ptr< RadiationPressureInterfaceSettings > vehicleRadiationPressureSettings_;

        // misc
        int intPlaceholder_;
        std::shared_ptr< electro_magnetism::RadiationPressureInterface > SRPInterface_;


    };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
                Eigen::Vector3d vehicleInitialPosition,
                Eigen::Vector3d vehicleInitialVelocity,
                nlohmann::json integratorSettingsJson,
                double initialEphemerisTime = 1.0E7,
                std::string terminationMethod = "nominalTimeTermination",
                double simulationLength = 10,
                std::string proximityTerminationBody1 = "Vehicle",
                std::string proximityTerminationBody2 = "Jupiter",
                double proximityCutoffAU = 1,
                double maxCPUTimeSecs = 60) :
                simulationPropBodies_(simulationPropBodies),
                vehicleInitialPosition_(vehicleInitialPosition),
                vehicleInitialVelocity_(vehicleInitialVelocity),
                integratorSettingsJson_(integratorSettingsJson),
                initialEphemerisTime_(initialEphemerisTime),
                terminationMethod_(terminationMethod),
                simulationLength_(simulationLength),
                proximityTerminationBody1_(proximityTerminationBody1),
                proximityTerminationBody2_(proximityTerminationBody2),
                proximityCutoffAU_(proximityCutoffAU),
                maxCPUTimeSecs_(maxCPUTimeSecs){
            std::cout << "beginning constructor" << std::endl;
            /////////////// Propagation Settings ///////////////////////
            // Create vehicle initial state from position and velocities
            vehicleInitialState_ << vehicleInitialPosition_, vehicleInitialVelocity_;

            // Set maximum CPU termination time for all cases
            terminationSettingsList.push_back( std::make_shared< PropagationCPUTimeTerminationSettings >(
                    maxCPUTimeSecs_ ) );

            // Set maximum simulation date to Jan 01 2100, since this is when ephemerides run out
            terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                    gen::year2MJDSeconds(2100) ) );

            // Define propagation termination conditions
            // Time termination - ends after a certain number of years defined from simulationLength
            if (terminationMethod_ == "nominalTimeTermination"){
                // Add default time termination to termination settings list
                double defaultFinalEphemTime = initialEphemerisTime_ + simulationLength_*365*24*60*60;
                terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                        defaultFinalEphemTime ) );
            }
            // Proximity termination - define a body and the proximity to it when termination happens. Also uses simulationTime as an upper limit of simulation
            else if (terminationMethod_ == "proximityTermination"){
                // Add time termination
                double maximumEphemTime = initialEphemerisTime_ + simulationLength_*365*24*60*60;
                terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                        maximumEphemTime ) );

                // Add proximity termination
                proximityTerminationDepVar_ =
                        std::make_shared< SingleDependentVariableSaveSettings >(
                                relative_distance_dependent_variable,
                                proximityTerminationBody1_,
                                proximityTerminationBody2_);

                terminationSettingsList.push_back(std::make_shared< PropagationDependentVariableTerminationSettings >(
                        proximityTerminationDepVar_,
                        proximityCutoffAU_*gen::AU,
                        true,
                        false));
            }
            else{
                std::cout << "ERROR: invalid termination setting" << std::endl;
                exit(1);
            }

            // Using above if-else really make the termination settings
            terminationSettings = std::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );
            std::cout << "Set termination settings" << std::endl;
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
            std::cout << "Set dependent variables to save" << std::endl;
            // Define settings for propagation of translational dynamics.
            propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >(
                            simulationPropBodies_.centralBodies, simulationPropBodies_.accelerationModelMap,
                            simulationPropBodies_.bodiesToPropagate, vehicleInitialState_, terminationSettings,
                            cowell, dependentVariablesToSave_ );

            std::cout << "Set proagator settings" << std::endl;

            ///////////////////// Integration Settings //////////////////////
            // TODO: Check general settings of integrator, and adjust as necessary (especially tolerances).
            // TODO: part 2, add some integrator settings to the json file
            std::string RKCoefficients = integratorSettingsJson["integrator"];
            RungeKuttaCoefficients::CoefficientSets coefficientSet = gen::getIntegratorCoefficientSet(RKCoefficients);
            
            integratorSettings =
                    std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                            initialEphemerisTime_,
                            integratorSettingsJson_["initialTimeStep"],
                            coefficientSet, //TODO: change this integrator to FILG if possible (ask dominic?)
                            integratorSettingsJson_["minimumStepSize"],
                            integratorSettingsJson_["maximumStepSize"],
                            integratorSettingsJson_["relativeErrorTolerance"],
                            integratorSettingsJson_["absoluteErrorTolerance"]);
            std::cout << "Set integrator settings" << std::endl;
        }

        propBodies getSimulationPropBodies(){
            return simulationPropBodies_;
        }

    protected:

        /////////////////////////// (Mandatory) Initialisation parameters ///////////////////////////////
        propBodies& simulationPropBodies_;
        Eigen::Vector3d vehicleInitialPosition_;
        Eigen::Vector3d vehicleInitialVelocity_;
        nlohmann::json integratorSettingsJson_;


        /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////
        double initialEphemerisTime_;
        std::string terminationMethod_;
        double simulationLength_;
        std::string proximityTerminationBody1_;
        std::string proximityTerminationBody2_;
        double proximityCutoffAU_;
        double maxCPUTimeSecs_;


        /////////////////////////// Other set parameters ///////////////////////////////
        Eigen::Vector6d vehicleInitialState_;
        std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave_;
        std::shared_ptr< SingleDependentVariableSaveSettings > proximityTerminationDepVar_;

    };

}

#endif //TUDATBUNDLE_UNIVERSAL_SETTINGS_H
