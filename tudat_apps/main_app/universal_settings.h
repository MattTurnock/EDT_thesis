//
// Created by matt on 29/02/2020.
//

#ifndef TUDATBUNDLE_UNIVERSAL_SETTINGS_H
#define TUDATBUNDLE_UNIVERSAL_SETTINGS_H


#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h>
#include <Tudat/External/SpiceInterface/spiceEphemeris.h>

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::gravitation;
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
                nlohmann::json jsonSimulationVariables) :
                vehicleConfig_(vehicleConfig),
                guidanceClass_(guidanceClass),
                bodyMap(baseBodyMap),
                jsonSimulationVariables_(jsonSimulationVariables)
                {

            // Constructor
            // Set various json sets from the simulation variables
            jsonBodiesToInclude_ = jsonSimulationVariables["Spice"]["bodiesToInclude"];
            gravitationalParameters_ = jsonSimulationVariables["constants"]["gravitationalParameters"];


            /////////// Initialise all bodies list with all possible bodies that can be used //////////
            allBodiesList_.push_back("Sun");
            allBodiesList_.push_back("Mercury");
            allBodiesList_.push_back("Venus");
            allBodiesList_.push_back("Earth");
            allBodiesList_.push_back("Mars");
            allBodiesList_.push_back("Jupiter");
            allBodiesList_.push_back("Saturn");
            allBodiesList_.push_back("Uranus");
            allBodiesList_.push_back("Neptune");


            /////////// Define body settings for simulation and create bodymap./////////////////////////
            // Create bodies list from json data
            for(std::size_t i=0; i<allBodiesList_.size(); ++i){
                intPlaceholder_ = jsonBodiesToInclude_[allBodiesList_[i]];
                if (jsonBodiesToInclude_[allBodiesList_[i]] == 1){
                    bodiesToCreate_.push_back(allBodiesList_[i]);
                }
            }

            std::cout << "Using bodies: " << std::endl;
            for(std::size_t i=0; i<bodiesToCreate_.size(); ++i){
                // Set string for body being used
                std::string bodyToCreate = bodiesToCreate_[i];
                std::cout << bodyToCreate << ", ";

                // Create variable for target body, and set observer body. Also set ref frame name
                std::string targetBodyName;
                std::string observerBodyName = "SSB";
                std::string referenceFrameName = "ECLIPJ2000";
                double gravitationalParameter = gravitationalParameters_[bodyToCreate];


                // Set body parameters - Uses planetary system gravitational data (eg earth) but barycenter coordinates
                if (bodyToCreate == "Sun") {
                    targetBodyName = "Sun";
                }
                else if (bodyToCreate == "Mercury"){
                    targetBodyName = "MERCURY_BARYCENTER";
                }
                else if (bodyToCreate == "Venus"){
                    targetBodyName = "VENUS_BARYCENTER";
                }
                else if (bodyToCreate == "Earth"){
                    targetBodyName = "EARTH_BARYCENTER";
                }
                else if (bodyToCreate == "Mars"){
                    targetBodyName = "MARS_BARYCENTER";
                }
                else if (bodyToCreate == "Jupiter"){
                    targetBodyName = "JUPITER_BARYCENTER";
                }
                else if (bodyToCreate == "Saturn"){
                    targetBodyName = "SATURN_BARYCENTER";
                }
                else if (bodyToCreate == "Uranus"){
                    targetBodyName = "URANUS_BARYCENTER";
                }
                else if (bodyToCreate == "Neptune"){
                    targetBodyName = "NEPTUNE_BARYCENTER";
                }

                // Create body in bodyMap
                bodyMap[bodyToCreate] = std::make_shared< Body >( );

                // Set ephemeris for new body
                bodyMap[ bodyToCreate ]->setEphemeris(
                        std::make_shared< SpiceEphemeris >( targetBodyName, observerBodyName,
                                                            true, true,
                                                            false, referenceFrameName ) );

                // Set gravitational parameter for new body
                bodyMap[ bodyToCreate ]->setGravityFieldModel( std::make_shared< GravityFieldModel >( gravitationalParameter ) );
            }
            std::cout << std::endl;

            // Add vehicle to bodyMap
            bodyMap["Vehicle"] = vehicleConfig.EDTBody;

            // Finalise body creation
            setGlobalFrameBodyEphemerides( bodyMap, "Sun", "ECLIPJ2000" );



            ///////////////// Create bodies acceleration map (includes perturbations) ///////////////////////
            // Add acceleration model settings for vehicle - ie thrust guidance settings. Uses boolean for creating acceleration map without vehicle thrust if wanted
            accelerationsOfVehicle_["Vehicle"].push_back(
                    std::make_shared<ThrustAccelerationSettings>(guidanceClass.thrustDirectionGuidanceSettings,
                                                                 guidanceClass.thrustMagnitudeSettings));
//                std::cout << "thrust acceleration added to vehicle" << std::endl; //TODO: remove me

            // Add acceleration model settings for celestials using for loop
            for(std::size_t i=0; i<bodiesToCreate_.size(); ++i){
                accelerationsOfVehicle_[bodiesToCreate_[i]].push_back(std::make_shared< AccelerationSettings >( central_gravity ) );
            }

            // Add acceleration model for cannonball SRP
            SRPReferenceArea_ = vehicleConfig_.getEffectiveSRPArea();
            SRPCoefficient_ = vehicleConfig_.getEffectiveRadiationPressureCoefficient();

//            std::cout << "Effective SRP radiation coefficient of vehicle: " << SRPCoefficient_ << std::endl; // TODO: remove me
//            std::cout << "Effective SRP reference area of vehicle       : " << SRPReferenceArea_ << std::endl; // TODO: remove me
//            std::cout << "Spacecraft Mass                               : " << vehicleConfig_.getVehicleMass() << std::endl; // TODO: remove me

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

        /////////////////////////// Other set parameters ///////////////////////////////
        // simulation variables (from json)
        nlohmann::json jsonSimulationVariables_;
        

        // Body variables
        std::vector<std::string> bodiesToCreate_;
        std::map< std::string, std::shared_ptr< BodySettings > > bodySettings_;
        std::vector<std::string> allBodiesList_;
        nlohmann::json gravitationalParameters_;

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

        // Constructor.
        propSettings(
                propBodies& simulationPropBodies,
                Eigen::Vector3d vehicleInitialPosition,
                Eigen::Vector3d vehicleInitialVelocity,
                nlohmann::json integratorSettingsJson,
                nlohmann::json terminationSettingsJson,
                double initialEphemerisYear = 2000):
                simulationPropBodies_(simulationPropBodies),
                vehicleInitialPosition_(vehicleInitialPosition),
                vehicleInitialVelocity_(vehicleInitialVelocity),
                integratorSettingsJson_(integratorSettingsJson),
                terminationSettingsJson_(terminationSettingsJson),
                initialEphemerisYear_(initialEphemerisYear){

            ////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////// Propagation Settings //////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            /////////////// Misc Initialisations ///////////////////////
            // Set termination Parameters from json data
            std::string terminationType = terminationSettingsJson_["terminationType"];
            double timeTerminationYears = terminationSettingsJson_["timeTerminationYears"];
            std::string proximityTerminationBody1 = terminationSettingsJson_["proximityTerminationBody1"];
            std::string proximityTerminationBody2 = terminationSettingsJson_["proximityTerminationBody2"];
            double proximityTerminationCutoffAU = terminationSettingsJson_["proximityTerminationCutoffAU"];
            double distanceTerminationAU = terminationSettingsJson_["distanceTerminationAU"];
            double maxCPUTimeTerminationSecs = terminationSettingsJson_["maxCPUTimeTerminationSecs"];
            double absoluteTerminationTimeYear = terminationSettingsJson_["absoluteTimeTerminationYear"];

            // Create vehicle initial state from position and velocities (also time)
            initialEphemerisTime_ = gen::year2MJDSeconds(initialEphemerisYear_);
//            std::cout << "INITIAL EPHEMERIS YEAR: " << initialEphemerisYear_ << std::endl;
//            std::cout << "INITIAL EPHEMERIS TIME: " << initialEphemerisTime_ << std::endl;
            vehicleInitialState_ << vehicleInitialPosition_, vehicleInitialVelocity_;

            /////////////// Termination Setttings ///////////////////////

            // Set maximum CPU termination time for all cases
            terminationSettingsList.push_back( std::make_shared< PropagationCPUTimeTerminationSettings >(
                    maxCPUTimeTerminationSecs ) );


            // Define propagation termination conditions
            // Time termination - ends after a certain number of years defined from simulationLength
            if (terminationType == "nominalTimeTermination"){
                // Add default time termination to termination settings list
                double defaultFinalEphemTime = initialEphemerisTime_ + timeTerminationYears * 365 * 24 * 60 * 60;
                terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                        defaultFinalEphemTime ) );
            }
            // Proximity termination - define a body and the proximity to it when termination happens. Also uses simulationTime as an upper limit of simulation
            // Distance termination - like proximity termination, but for when vehicle is further than distance x
            else if (terminationType == "proximityTermination" or terminationType=="distanceTermination") {

                // Set variables for differences between proximity termination (when vehicle is closer than x) vs distance termination (when vehicle is further than x)
                bool useAsLowerLimit;
                double distanceToUse;
                if (terminationType == "proximityTermination") {
                    useAsLowerLimit = true;
                    distanceToUse = proximityTerminationCutoffAU;
                } else if (terminationType == "distanceTermination") {
                    useAsLowerLimit = false;
                    distanceToUse = distanceTerminationAU;
                }


                // Add proximity termination
                proximityTerminationDepVar_ =
                        std::make_shared<SingleDependentVariableSaveSettings>(
                                relative_distance_dependent_variable,
                                proximityTerminationBody1,
                                proximityTerminationBody2);

                terminationSettingsList.push_back(std::make_shared<PropagationDependentVariableTerminationSettings>(
                        proximityTerminationDepVar_,
                        distanceToUse * gen::AU,
                        useAsLowerLimit,
                        false));
            }
            else{
                std::cout << "ERROR: invalid termination setting" << std::endl;
                exit(1);
            }

            // Add in absolute time terminator regardless of setting (avoids epehemeris errors)
            terminationSettingsList.push_back(std::make_shared<PropagationCustomTerminationSettings>(
                    std::bind(&gen::absoluteTimeTermination, std::placeholders::_1, absoluteTerminationTimeYear)
                    ));

            // Using above if-else really make the termination settings
            terminationSettings = std::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );

            /////////////// Set dependent variables to save ///////////////////////

            // Add Keplerian state
            dependentVariablesList.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >( keplerian_state_dependent_variable,
                            "Vehicle" ,
                            "Sun" ) );

            // Add Relative distance
            dependentVariablesList.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(relative_distance_dependent_variable,
                            "Vehicle", "Sun"));

            // Add Single acceleration variable for SRP
            dependentVariablesList.push_back(
                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            cannon_ball_radiation_pressure,
                            "Vehicle",
                            "Sun",
                            false)
                    );

            // Finalise dependent variables to save
            dependentVariablesToSave_ =
                    std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );


            /////////////// Finalise propagation settings ///////////////////////

            // Define settings for propagation of translational dynamics.
            propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >(
                            simulationPropBodies_.centralBodies, simulationPropBodies_.accelerationModelMap,
                            simulationPropBodies_.bodiesToPropagate, vehicleInitialState_, terminationSettings,
                            cowell, dependentVariablesToSave_ );


            ////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////// Integration Settings ////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            // TODO: Check general settings of integrator, and adjust as necessary (especially tolerances).
            std::string RKCoefficients = integratorSettingsJson["integrator"];

            if (RKCoefficients == "RK4"){
                integratorSettings = std::make_shared<IntegratorSettings< > >(
                        rungeKutta4,
                        initialEphemerisTime_,
                        integratorSettingsJson_["initialTimeStep"]
                        );
            }

            else {
                RungeKuttaCoefficients::CoefficientSets coefficientSet = gen::getIntegratorCoefficientSet(
                        RKCoefficients);

                integratorSettings =
                        std::make_shared<RungeKuttaVariableStepSizeSettings<> >(
                                initialEphemerisTime_,
                                integratorSettingsJson_["initialTimeStep"],
                                coefficientSet, //TODO: change this integrator to FILG if possible (ask dominic?)
                                integratorSettingsJson_["minimumStepSize"],
                                integratorSettingsJson_["maximumStepSize"],
                                integratorSettingsJson_["relativeErrorTolerance"],
                                integratorSettingsJson_["absoluteErrorTolerance"]);
            }
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
        nlohmann::json terminationSettingsJson_;


        /////////////////////////// (Optional) Initialisation parameters ///////////////////////////////
        double initialEphemerisYear_;


        /////////////////////////// Other set parameters ///////////////////////////////
        Eigen::Vector6d vehicleInitialState_;
        std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave_;
        std::shared_ptr< SingleDependentVariableSaveSettings > proximityTerminationDepVar_;
        double initialEphemerisTime_;




    };

}

#endif //TUDATBUNDLE_UNIVERSAL_SETTINGS_H
