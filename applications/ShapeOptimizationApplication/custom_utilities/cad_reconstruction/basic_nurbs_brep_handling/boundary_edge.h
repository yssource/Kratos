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

#ifndef BOUNDARY_EDGE_H
#define BOUNDARY_EDGE_H

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
#include "control_point.h"
#include "b_spline_utilities.h"

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

class BoundaryEdge
{
public:
	///@name Type Definitions
	///@{

	// For matrix / vector operations
	typedef std::vector<double> DoubleVector;
	typedef std::vector<ControlPoint> ControlPointVector;

	///@}

	/// Pointer definition of BoundaryEdge
	//    KRATOS_CLASS_POINTER_DEFINITION[BoundaryEdge];

	/// Default constructor.
	BoundaryEdge(DoubleVector knot_vector_u, unsigned int p, ControlPointVector control_points) 
	: m_knot_vector_u(knot_vector_u),
	  m_p(p),
	  m_control_points(control_points)
	{
		m_n_u = m_knot_vector_u.size() - m_p - 1;
	}

	/// Destructor.
	virtual ~BoundaryEdge()
	{
	}

	//  #####################################################################################
	// #######################################################################################
	//#
	//#                  +++++++++++++++++++++++++++++++++++++++
	//#                  +++  EvaluateCurvePoint +++
	//#                  +++++++++++++++++++++++++++++++++++++++
	//#
	///   \details    returns the cartesian coordinates (x,y,z) for a specific point 
	///               located on the NURBS curve C(u=fixed)
	///               Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
	///               Algorithm A4.1
	///
	/// ======================================================================================         
	///   \param[in]   rCurvePoint    point of interest in x,y,z coordinates
	///   \param[in]  _uoi    knot of interest (reals)
	///
	/// ======================================================================================
	///  \author     D. Baumgärtner 
	//
	//########################################################################################
	void EvaluateCurvePoint(Point<3>& rCurvePoint, double _uoi) 
	{
		const unsigned int R_Dim = 3;

		Vector tmpHomCtrlPt;
		Vector nBasis;
		Vector resulting_point;
		Vector homPoi;
		matrix<double> homCtrlPts;
		
		unsigned int span_u = find_Knot_Span(m_knot_vector_u,_uoi,m_p,m_n_u);
		eval_nonzero_basis_function(nBasis,m_knot_vector_u,_uoi,span_u,m_p);
		homCtrlPts.resize(R_Dim+1,m_p+1);
		for(unsigned int i=0;i<=m_p;i++)
		{
			int control_point_index = span_u-m_p+i;

			// tmpHomCtrlPt = m_control_points[control_point_index].GetWeight();
			// tmpHomCtrlPt = Ctrl_Pt[span-m_p+i]->Get_Ctrl_Pt_Coo_w();
			Vector tmpHomCtrlPt(R_Dim+1);

			double w = m_control_points[control_point_index].GetWeight();
			double x = m_control_points[control_point_index].GetX();
			double y = m_control_points[control_point_index].GetY();
			double z = m_control_points[control_point_index].GetZ();

			tmpHomCtrlPt[0] =  x*w;
			tmpHomCtrlPt[1] =  y*w;
			tmpHomCtrlPt[2] =  z*w;
			tmpHomCtrlPt[3] =  w;

			for(unsigned int j=0;j<=R_Dim;j++)
			{
				homCtrlPts(j,i) = tmpHomCtrlPt(j);
			}
		}
		homPoi = prod(homCtrlPts,nBasis);
		resulting_point = (1/homPoi(R_Dim))*homPoi;

		rCurvePoint.X() = resulting_point[0];
		rCurvePoint.Y() = resulting_point[1];
		rCurvePoint.Z() = resulting_point[2];
		
		if (std::abs(resulting_point(R_Dim)-1.00) > m_epsilon) 
		{
			KRATOS_THROW_ERROR( std::logic_error, "NURBS 1D: evalutation curve point failed!!!","");
		}
	}

	// --------------------------------------------------------------------------
	DoubleVector& GetKnotVectorU()
	{
		return m_knot_vector_u;
	}

	// --------------------------------------------------------------------------
	unsigned int GetPolynomialDegree()
	{
		return m_p;
	}
	

	// --------------------------------------------------------------------------
	ControlPointVector& GetControlPoints()
	{
		return m_control_points;
	}

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "BoundaryEdge";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "BoundaryEdge";
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
	DoubleVector m_knot_vector_u;
	unsigned int m_p;
	ControlPointVector m_control_points;
	unsigned int m_n_u; // number of control points in u-direction
	double m_epsilon = 1e-10; // Tolerance value

	// ==============================================================================
	// General working arrays
	// ==============================================================================
	/// Assignment operator.
	//      BoundaryEdge& operator=[BoundaryEdge const& rOther];

	/// Copy constructor.
	//      BoundaryEdge[BoundaryEdge const& rOther];

}; // Class BoundaryEdge

} // namespace Kratos.

#endif // BOUNDARY_EDGE_H
