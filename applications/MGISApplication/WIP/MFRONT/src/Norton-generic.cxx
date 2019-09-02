/*!
* \file   Norton-generic.cxx
* \brief  This file implements the umat interface for the Norton behaviour law
* \author Helfer Thomas
* \date   23 / 11 / 06
*/

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif /* NOMINMAX */
#include <windows.h>
#ifdef small
#undef small
#endif /* small */
#endif /* _WIN32 */

#ifndef MFRONT_SHAREDOBJ
#define MFRONT_SHAREDOBJ TFEL_VISIBILITY_EXPORT
#endif /* MFRONT_SHAREDOBJ */

#include<iostream>
#include<cstdlib>
#include"TFEL/Material/OutOfBoundsPolicy.hxx"
#include"TFEL/Material/Norton.hxx"
#include"MFront/GenericBehaviour/Integrate.hxx"

#include"MFront/GenericBehaviour/Norton-generic.hxx"

static tfel::material::OutOfBoundsPolicy&
Norton_getOutOfBoundsPolicy(){
using namespace tfel::material;
static OutOfBoundsPolicy policy = None;
return policy;
}

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */

MFRONT_SHAREDOBJ const char* 
Norton_mfront_ept = "Norton";

MFRONT_SHAREDOBJ const char* 
Norton_tfel_version = "3.3.0-dev";

MFRONT_SHAREDOBJ unsigned short Norton_mfront_mkt = 1u;

MFRONT_SHAREDOBJ const char *
Norton_mfront_interface = "Generic";

MFRONT_SHAREDOBJ const char *
Norton_src = "Norton.mfront";

MFRONT_SHAREDOBJ unsigned short Norton_nModellingHypotheses = 5u;

MFRONT_SHAREDOBJ const char * 
Norton_ModellingHypotheses[5u] = {"AxisymmetricalGeneralisedPlaneStrain",
"Axisymmetrical",
"PlaneStrain",
"GeneralisedPlaneStrain",
"Tridimensional"};

MFRONT_SHAREDOBJ unsigned short Norton_nMainVariables = 1;
MFRONT_SHAREDOBJ unsigned short Norton_nGradients = 1;

MFRONT_SHAREDOBJ int Norton_GradientsTypes[1] = {1};
MFRONT_SHAREDOBJ const char * Norton_Gradients[1] = {"Strain"};
MFRONT_SHAREDOBJ unsigned short Norton_nThermodynamicForces = 1;

MFRONT_SHAREDOBJ int Norton_ThermodynamicForcesTypes[1] = {1};
MFRONT_SHAREDOBJ const char * Norton_ThermodynamicForces[1] = {"Stress"};
MFRONT_SHAREDOBJ unsigned short Norton_nTangentOperatorBlocks = 2;

MFRONT_SHAREDOBJ const char * Norton_TangentOperatorBlocks[2] = {"Stress",
"Strain"};
MFRONT_SHAREDOBJ unsigned short Norton_BehaviourType = 1u;

MFRONT_SHAREDOBJ unsigned short Norton_BehaviourKinematic = 1u;

MFRONT_SHAREDOBJ unsigned short Norton_SymmetryType = 0u;

MFRONT_SHAREDOBJ unsigned short Norton_ElasticSymmetryType = 0u;

MFRONT_SHAREDOBJ unsigned short Norton_UsableInPurelyImplicitResolution = 1;

MFRONT_SHAREDOBJ unsigned short Norton_nMaterialProperties = 4u;

MFRONT_SHAREDOBJ const char *Norton_MaterialProperties[4u] = {"NortonCoefficient",
"NortonExponent",
"YoungModulus",
"PoissonRatio"};

MFRONT_SHAREDOBJ unsigned short Norton_nInternalStateVariables = 2;
MFRONT_SHAREDOBJ const char * Norton_InternalStateVariables[2] = {"ElasticStrain",
"EquivalentViscoplasticStrain"};
MFRONT_SHAREDOBJ int Norton_InternalStateVariablesTypes [] = {1,0};

MFRONT_SHAREDOBJ unsigned short Norton_nExternalStateVariables = 0;
MFRONT_SHAREDOBJ const char * const * Norton_ExternalStateVariables = nullptr;

MFRONT_SHAREDOBJ unsigned short Norton_nParameters = 5;
MFRONT_SHAREDOBJ const char * Norton_Parameters[5] = {"minimal_time_step_scaling_factor",
"maximal_time_step_scaling_factor","theta","epsilon","iterMax"};
MFRONT_SHAREDOBJ int Norton_ParametersTypes [] = {0,0,0,0,2};

MFRONT_SHAREDOBJ double Norton_minimal_time_step_scaling_factor_ParameterDefaultValue = 0.1;

MFRONT_SHAREDOBJ double Norton_maximal_time_step_scaling_factor_ParameterDefaultValue = 1.79769e+308;

MFRONT_SHAREDOBJ double Norton_theta_ParameterDefaultValue = 0.5;

MFRONT_SHAREDOBJ double Norton_epsilon_ParameterDefaultValue = 1e-08;

MFRONT_SHAREDOBJ unsigned short Norton_iterMax_ParameterDefaultValue  = 100;

MFRONT_SHAREDOBJ unsigned short Norton_requiresStiffnessTensor = 0;
MFRONT_SHAREDOBJ unsigned short Norton_requiresThermalExpansionCoefficientTensor = 0;

MFRONT_SHAREDOBJ void
Norton_setOutOfBoundsPolicy(const int p){
if(p==0){
Norton_getOutOfBoundsPolicy() = tfel::material::None;
} else if(p==1){
Norton_getOutOfBoundsPolicy() = tfel::material::Warning;
} else if(p==2){
Norton_getOutOfBoundsPolicy() = tfel::material::Strict;
} else {
std::cerr << "Norton_setOutOfBoundsPolicy: invalid argument\n";
}
}

MFRONT_SHAREDOBJ int
Norton_setParameter(const char *const key,const double value){
using tfel::material::NortonParametersInitializer;
auto& i = NortonParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int
Norton_setUnsignedShortParameter(const char *const key,const unsigned short value){
using tfel::material::NortonParametersInitializer;
auto& i = NortonParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int Norton_AxisymmetricalGeneralisedPlaneStrain(MFront_GB_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr const auto h = ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN;
using Behaviour = Norton<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, Norton_getOutOfBoundsPolicy());
return r;
} // end of Norton_AxisymmetricalGeneralisedPlaneStrain

MFRONT_SHAREDOBJ int Norton_Axisymmetrical(MFront_GB_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr const auto h = ModellingHypothesis::AXISYMMETRICAL;
using Behaviour = Norton<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, Norton_getOutOfBoundsPolicy());
return r;
} // end of Norton_Axisymmetrical

MFRONT_SHAREDOBJ int Norton_PlaneStrain(MFront_GB_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr const auto h = ModellingHypothesis::PLANESTRAIN;
using Behaviour = Norton<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, Norton_getOutOfBoundsPolicy());
return r;
} // end of Norton_PlaneStrain

MFRONT_SHAREDOBJ int Norton_GeneralisedPlaneStrain(MFront_GB_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr const auto h = ModellingHypothesis::GENERALISEDPLANESTRAIN;
using Behaviour = Norton<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, Norton_getOutOfBoundsPolicy());
return r;
} // end of Norton_GeneralisedPlaneStrain

MFRONT_SHAREDOBJ int Norton_Tridimensional(MFront_GB_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr const auto h = ModellingHypothesis::TRIDIMENSIONAL;
using Behaviour = Norton<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, Norton_getOutOfBoundsPolicy());
return r;
} // end of Norton_Tridimensional

#ifdef __cplusplus
}
#endif /* __cplusplus */

