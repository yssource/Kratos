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

#ifndef BREP_ELEMENT_H
#define BREP_ELEMENT_H

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
#include "brep_gauss_point.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

	typedef std::vector<BREPGaussPoint> BREPGaussPointVector;

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

class BREPElement
{
public:
	///@name Type Definitions
	///@{

	///@}

	/// Pointer definition of BREPElement
	//    KRATOS_CLASS_POINTER_DEFINITION[BREPElement];

	/// Default constructor.
	BREPElement(unsigned int element_id, unsigned int edge_id, BREPGaussPointVector gauss_points, bool has_coupling_condition, bool has_dirichlet_condition) 
	: m_element_id(element_id),
	  m_edge_id(edge_id),
	  m_gauss_points(gauss_points),
	  m_has_coupling_condition(has_coupling_condition),
	  m_has_dirichlet_condition(has_dirichlet_condition)
	{
	}

	/// Destructor.
	virtual ~BREPElement()
	{
	}
	// --------------------------------------------------------------------------
	unsigned int GetElementId()
	{
		return m_element_id;
	}	

	// --------------------------------------------------------------------------
	unsigned int GetEdgeId()
	{
		return m_edge_id;
	}		

	// --------------------------------------------------------------------------
	BREPGaussPointVector GetGaussPoints()
	{
		return m_gauss_points;
	}	

	// --------------------------------------------------------------------------
	bool HasCouplingCondition()
	{
		return m_has_coupling_condition;
	}	

	// --------------------------------------------------------------------------
	bool HasDirichletCondition()
	{
		return m_has_dirichlet_condition;
	}			

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "BREPElement";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "BREPElement";
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
	unsigned int m_element_id;
	unsigned int m_edge_id;
	BREPGaussPointVector m_gauss_points;
	bool m_has_coupling_condition;
	bool m_has_dirichlet_condition;

	// ==============================================================================
	// General working arrays
	// ==============================================================================
	/// Assignment operator.
	//      BREPElement& operator=[BREPElement const& rOther];

	/// Copy constructor.
	//      BREPElement[BREPElement const& rOther];

}; // Class BREPElement

} // namespace Kratos.

#endif // BREP_ELEMENT_H
