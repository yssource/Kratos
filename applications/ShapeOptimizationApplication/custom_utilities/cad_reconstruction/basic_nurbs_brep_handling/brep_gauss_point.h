// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 [Released on march 05, 2007].

 Copyright [c] 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  [the
 "Software"], to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
//==============================================================================
//
//   Project Name:        KratosShape                            $
//   Created by:          $Author:    daniel.baumgaertner@tum.de $
//   Last modified by:    $Co-Author: daniel.baumgaertner@tum.de $
//   Date:                $Date:                   December 2016 $
//   Revision:            $Revision:                         0.0 $
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
					double weight,
					Vector location,
					Vector tangent )
	: m_master_patch_id(master_patch_id),
	  m_gauss_point_id(gauss_point_id),
	  m_weight(weight),
	  m_location(location),
	  m_tangent(tangent)
	{
	}

	// Alternative constructor to carry neighbour relations
	BREPGaussPoint(	unsigned int master_patch_id,
					unsigned int slave_patch_id,
					unsigned int gauss_point_id,
					double weight,
					Vector location,
					Vector tangent,
					Vector location_slave,
					Vector tangent_slave )
	: m_master_patch_id(master_patch_id),
	  m_slave_patch_id(slave_patch_id),
	  m_gauss_point_id(gauss_point_id),
	  m_weight(weight),
	  m_location(location),
	  m_tangent(tangent),
	  m_slave_location(location_slave),
	  m_slave_tangent(tangent_slave)
	{
	}

	/// Destructor.
	virtual ~BREPGaussPoint()
	{
	}

	// --------------------------------------------------------------------------
	unsigned int GetId()
	{
		return m_gauss_point_id;
	}	

	// --------------------------------------------------------------------------
	unsigned int GetPatchId()
	{
		return m_master_patch_id;
	}	

	// --------------------------------------------------------------------------
	unsigned int GetSlavePatchId()
	{
		return m_slave_patch_id;
	}	

	// --------------------------------------------------------------------------
	double GetWeight()
	{
		return m_weight;
	}			

	// --------------------------------------------------------------------------
	Vector GetLocation()
	{
		return m_location;
	}		

	// --------------------------------------------------------------------------
	Vector GetSlaveLocation()
	{
		return m_slave_location;
	}	

	// --------------------------------------------------------------------------
	Vector GetTangent()
	{
		return m_tangent;
	}	

	// --------------------------------------------------------------------------
	Vector GetSlaveTangent()
	{
		return m_slave_tangent;
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
	unsigned int m_master_patch_id;
	unsigned int m_slave_patch_id;
	unsigned int m_gauss_point_id;
	double m_weight;
	Vector m_location = ZeroVector(2);
	Vector m_tangent = ZeroVector(2);
	Vector m_slave_location = ZeroVector(2);
	Vector m_slave_tangent = ZeroVector(2);

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
