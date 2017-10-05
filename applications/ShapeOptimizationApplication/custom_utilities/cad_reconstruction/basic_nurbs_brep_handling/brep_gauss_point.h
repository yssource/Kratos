// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef BREP_GAUSS_POINT
#define BREP_GAUSS_POINT

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

 */

class BREPGaussPoint
{
public:
	///@name Type Definitions
	///@{

	///@}

	/// Pointer definition of BREPGaussPoint
	//    KRATOS_CLASS_POINTER_DEFINITION[BREPGaussPoint];

	/// Default constructor.
	BREPGaussPoint(	unsigned int master_patch_id,
					unsigned int gauss_point_id,
					double integration_weight,
					array_1d<double,2> master_location,
					array_1d<double,2> master_tangent )
	: mMasterPatchId(master_patch_id),
	  mGaussPointId(gauss_point_id),
	  mIntegartionWeight(integration_weight),
	  mLocationOnMasterInParameterSpace(master_location),
	  mTangentOnMasterInParametersSpace(master_tangent)
	{
	}

	// Alternative constructor to carry neighbour relations
	BREPGaussPoint(	unsigned int master_patch_id,
					unsigned int slave_patch_id,
					unsigned int gauss_point_id,
					double integration_weight,
					array_1d<double,2> master_location,
					array_1d<double,2> master_tangent,
					array_1d<double,2> location_slave,
					array_1d<double,2> tangent_slave )
	: mMasterPatchId(master_patch_id),
	  mSlavePatchId(slave_patch_id),
	  mGaussPointId(gauss_point_id),
	  mIntegartionWeight(integration_weight),
	  mLocationOnMasterInParameterSpace(master_location),
	  mTangentOnMasterInParametersSpace(master_tangent),
	  mLocationOnSlaveInParameterSpace(location_slave),
	  mTangentOnSlaveInParametersSpace(tangent_slave)
	{
	}

	/// Destructor.
	virtual ~BREPGaussPoint()
	{
	}

	// --------------------------------------------------------------------------
	unsigned int GetId()
	{
		return mGaussPointId;
	}	

	// --------------------------------------------------------------------------
	unsigned int GetMasterPatchId()
	{
		return mMasterPatchId;
	}	

	// --------------------------------------------------------------------------
	unsigned int GetSlavePatchId()
	{
		return mSlavePatchId;
	}	

	// --------------------------------------------------------------------------
	double GetWeight()
	{
		return mIntegartionWeight;
	}			

	// --------------------------------------------------------------------------
	array_1d<double,2> GetLocationOnMasterInParameterSpace()
	{
		return mLocationOnMasterInParameterSpace;
	}		

	// --------------------------------------------------------------------------
	array_1d<double,2> GetLocationOnSlaveInParameterSpace()
	{
		return mLocationOnSlaveInParameterSpace;
	}	

	// --------------------------------------------------------------------------
	array_1d<double,2> GetTangentOnMasterInParameterSpace()
	{
		return mTangentOnMasterInParametersSpace;
	}	

	// --------------------------------------------------------------------------
	array_1d<double,2> GetTangentOnSlaveInParameterSpace()
	{
		return mTangentOnSlaveInParametersSpace;
	}				

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "BREPGaussPoint";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "BREPGaussPoint";
	}

	// ==============================================================================
	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}


private:
	// ==============================================================================
	// Initialized by class constructor
	// ==============================================================================
	unsigned int mMasterPatchId;
	unsigned int mSlavePatchId;
	unsigned int mGaussPointId;
	double mIntegartionWeight;
	array_1d<double,2> mLocationOnMasterInParameterSpace;
	array_1d<double,2> mTangentOnMasterInParametersSpace;
	array_1d<double,2> mLocationOnSlaveInParameterSpace;
	array_1d<double,2> mTangentOnSlaveInParametersSpace;

	// ==============================================================================
	// General working arrays
	// ==============================================================================
	/// Assignment operator.
	//      BREPGaussPoint& operator=[BREPGaussPoint const& rOther];

	/// Copy constructor.
	//      BREPGaussPoint[BREPGaussPoint const& rOther];

}; // Class BREPGaussPoint

} // namespace Kratos.

#endif // BREP_GAUSS_POINT
