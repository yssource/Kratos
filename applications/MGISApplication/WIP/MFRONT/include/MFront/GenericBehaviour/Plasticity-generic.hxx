/*!
* \file   Plasticity-generic.hxx
* \brief  This file declares the umat interface for the Plasticity behaviour law
* \author Helfer Thomas
* \date   23 / 11 / 06
*/

#ifndef LIB_GENERIC_PLASTICITY_HXX
#define LIB_GENERIC_PLASTICITY_HXX

#include"TFEL/Config/TFELConfig.hxx"
#include"MFront/GenericBehaviour/BehaviourData.h"

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

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */

MFRONT_SHAREDOBJ void
Plasticity_setOutOfBoundsPolicy(const int);

MFRONT_SHAREDOBJ int
Plasticity_setParameter(const char *const,const double);

MFRONT_SHAREDOBJ int
Plasticity_setUnsignedShortParameter(const char *const,const unsigned short);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int Plasticity_AxisymmetricalGeneralisedPlaneStrain(MFront_GB_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int Plasticity_Axisymmetrical(MFront_GB_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int Plasticity_PlaneStrain(MFront_GB_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int Plasticity_GeneralisedPlaneStrain(MFront_GB_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int Plasticity_Tridimensional(MFront_GB_BehaviourData* const);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* LIB_GENERIC_PLASTICITY_HXX */
