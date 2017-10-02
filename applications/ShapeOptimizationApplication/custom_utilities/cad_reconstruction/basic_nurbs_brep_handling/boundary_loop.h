// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef BOUNDARY_LOOP_H
#define BOUNDARY_LOOP_H

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
#include "boundary_edge.h"

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

class BoundaryLoop
{
public:
	///@name Type Definitions
	///@{

	typedef std::vector<BoundaryEdge> BoundaryEdgeVector;

	///@}

	/// Pointer definition of BoundaryLoop
	//    KRATOS_CLASS_POINTER_DEFINITION[BoundaryLoop];

	/// Default constructor.
	BoundaryLoop(BoundaryEdgeVector boundary_edges, bool is_inner_loop) 
	: m_boundary_edges(boundary_edges),
	  m_is_inner_loop(is_inner_loop)
	{
		// Create polygon for inside / outside check
		CreatePolygon();
	}

	/// Destructor.
	virtual ~BoundaryLoop()
	{
	}

	// --------------------------------------------------------------------------
	void CreatePolygon()	
	{
		// Variables needed
		unsigned int np_per_edge = 500;
		unsigned int counter=0;

		// Loop over all boundary edges creating the closed boundary loop
		for (BoundaryEdgeVector::iterator edge_i = m_boundary_edges.begin(); edge_i != m_boundary_edges.end(); ++edge_i)
		{
			DoubleVector& knot_vec_u = edge_i->GetKnotVectorU();

			double u_min = knot_vec_u[0];
			double u_max = knot_vec_u[knot_vec_u.size()-1];
			double delta_u = (u_max-u_min) / np_per_edge;

			// Resize to incude points on current edge
			m_boundary_polygon.resize(m_boundary_polygon.size()+np_per_edge);

			// Add points of edge to polygon
			double u_i = u_min;
			for(unsigned int i=0;i<np_per_edge;i++)
			{
				u_i += delta_u;
				Point<3> curve_point;

				edge_i->EvaluateCurvePoint(curve_point,u_i);

				m_boundary_polygon[counter][0]=curve_point.X();
				m_boundary_polygon[counter][1]=curve_point.Y();

				counter++;
			}			
		}
	}

	// --------------------------------------------------------------------------
	BoundaryEdgeVector& GetBoundaryEdges()
	{
		return m_boundary_edges;
	}

	// --------------------------------------------------------------------------
	std::vector<array_1d<double, 2>>& GetBoundaryPolygon()
	{
		return m_boundary_polygon;
	}	

	// --------------------------------------------------------------------------
	bool IsInnerLoop()
	{
		return m_is_inner_loop;
	}
	

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "BoundaryLoop";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "BoundaryLoop";
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
	BoundaryEdgeVector m_boundary_edges;
	bool m_is_inner_loop;
	std::vector<array_1d<double, 2>> m_boundary_polygon;

	// ==============================================================================
	// General working arrays
	// ==============================================================================
	/// Assignment operator.
	//      BoundaryLoop& operator=[BoundaryLoop const& rOther];

	/// Copy constructor.
	//      BoundaryLoop[BoundaryLoop const& rOther];

}; // Class BoundaryLoop

} // namespace Kratos.

#endif // BOUNDARY_LOOP_H
