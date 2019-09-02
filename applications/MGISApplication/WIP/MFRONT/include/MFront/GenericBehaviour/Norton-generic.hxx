/*!
* \file   Norton-generic.hxx
* \brief  This file declares the umat interface for the Norton behaviour law
* \author Helfer Thomas
* \date   23 / 11 / 06
*/

#ifndef LIB_GENERIC_NORTON_HXX
#define LIB_GENERIC_NORTON_HXX

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
Norton_setOutOfBoundsPolicy(const int);

MFRONT_SHAREDOBJ int
Norton_setParameter(const char *const,const double);

MFRONT_SHAREDOBJ int
Norton_setUnsignedShortParameter(const char *const,const unsigned short);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int Norton_AxisymmetricalGeneralisedPlaneStrain(MFront_GB_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int Norton_Axisymmetrical(MFront_GB_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int Norton_PlaneStrain(MFront_GB_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int Norton_GeneralisedPlaneStrain(MFront_GB_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int Norton_Tridimensional(MFront_GB_BehaviourData* const);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* LIB_GENERIC_NORTON_HXX */
