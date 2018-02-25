// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef PATCH_H
#define PATCH_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "control_point.h"
#include "nurbs_surface.h"
#include "boundary_loop.h"

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

class Patch
{
public:
	///@name Type Definitions
	///@{

	///@}

	// Pointer definition of Patch
	KRATOS_CLASS_POINTER_DEFINITION(Patch);

	/// Default constructor.
	Patch(unsigned int patch_id, NURBSSurface nurbs_surface, std::vector<BoundaryLoop> boundary_loops )
	: m_patch_id(patch_id),
	  m_nurbs_surface(nurbs_surface),
	  m_boundary_loops(boundary_loops)
	{
	}

	/// Destructor.
	virtual ~Patch()
	{
	}

	// --------------------------------------------------------------------------
	bool IsPointInside(array_1d<double, 2> &point_of_interest)
	{
		// Boost is used to check whether point of interest is inside given polygon or not

		// Type definitions to use boost functionalities
		typedef boost::geometry::model::d2::point_xy<double> point_type;
		typedef boost::geometry::model::polygon<point_type> polygon_type;

		// We assume point is inside, check all boundary loops if this is true. If this is not true for a single loop, then point is considered outside of this patch
		point_type poi(point_of_interest[0],point_of_interest[1]);
		bool is_inside = true;

		// Loop over all boundary loops of current patch
		for(auto & loop_i : m_boundary_loops)
		{
			// Initialize necessary variables
			polygon_type poly;
			std::vector<array_1d<double, 2>>& boundary_polygon = loop_i.GetBoundaryPolygon();

			// Prepare polygon for boost
			for(unsigned int i=0;i<boundary_polygon.size();i++)
				boost::geometry::append( boost::geometry::exterior_ring(poly),
										 boost::geometry::make<point_type>(boundary_polygon[i][0], boundary_polygon[i][1]));
			if(boundary_polygon.size()>0)
				boost::geometry::append( boost::geometry::exterior_ring(poly),
										 boost::geometry::make<point_type>(boundary_polygon[0][0], boundary_polygon[0][1]));

			// Check inside or outside
			is_inside = boost::geometry::within(poi, poly);

			// Boost does not consider the polygon direction, it always assumes inside as in the interior of the closed polygon.
			// If the CAD loop is an inner loop, however, the area which is considered inner is the unbounded side of the closed polygon.
			// So we toggle the results from the within search to be correct.
			if(loop_i.IsInnerLoop())
				is_inside = !is_inside;

			// If a point is considered outside in one loop, it is outside in general, hence we can break the loop
			if(!is_inside)
				break;
		}

		return is_inside;
	}

	// --------------------------------------------------------------------------
	void EvaluateSurfacePoint( array_1d<double,2>& parameter_values, Point<3>& rSurfacePoint )
	{
		m_nurbs_surface.EvaluateSurfacePoint( parameter_values[0], parameter_values[1], rSurfacePoint );
	}

	// --------------------------------------------------------------------------
	void EvaluateSurfaceDisplacement( array_1d<double,2>& parameter_values, array_1d<double,3>& rSurfaceDisplacement )
	{
		m_nurbs_surface.EvaluateSurfaceDisplacement( parameter_values[0], parameter_values[1], rSurfaceDisplacement );
	}

	// --------------------------------------------------------------------------
	std::vector<double> EvaluateNURBSFunctions( array_1d<int,2>& parameter_spans, array_1d<double,2>& parameter_values )
	{
		return m_nurbs_surface.EvaluateNURBSFunctions( parameter_spans[0], parameter_spans[1], parameter_values[0], parameter_values[1] );
	}

	// --------------------------------------------------------------------------
	void ComputeLocalCSY( array_1d<int,2>& parameter_spans,
									 array_1d<double,2>& parameter_values,
									 array_1d<double,2>& _par_g1,
								     Vector& _t1,
								     Vector& _t2,
								     Vector& _t3 )
	{
		m_nurbs_surface.ComputeLocalCSY( parameter_spans[0],
													parameter_spans[1],
													parameter_values[0],
													parameter_values[1],
													_par_g1,
													_t1, _t2, _t3 );
	}

	// --------------------------------------------------------------------------
	void ComputeVariationOfLocalCSY( array_1d<int,2>& parameter_spans,
									 array_1d<double,2>& parameter_values,
									 array_1d<double,2>& _par_g1,
								     Vector& _t1,
								     Vector& _t2,
								     Vector& _t3,
								     std::vector<Vector>& _t1r,
								     std::vector<Vector>& _t2r,
									 std::vector<Vector>& _t3r )
	{
		m_nurbs_surface.ComputeVariationOfLocalCSY( parameter_spans[0],
													parameter_spans[1],
													parameter_values[0],
													parameter_values[1],
													_par_g1,
													_t1, _t2, _t3,
													_t1r, _t2r, _t3r );
	}

	// --------------------------------------------------------------------------
	void ComputeSecondVariationOfLocalCSY( array_1d<int,2>& parameter_spans,
									       array_1d<double,2>& parameter_values,
									       array_1d<double,2>& _par_g1,
									  	   Vector& _t1,
										   Vector& _t2,
										   Vector& _t3,
										   Vector& _t1_der,
										   Vector& _t2_der,
										   Vector& _t3_der,
										   std::vector<Vector >& _t1_r,
										   std::vector<Vector >& _t2_r,
										   std::vector<Vector >& _t3_r,
										   std::vector<Vector >& _t1_der_r,
										   std::vector<Vector >& _t2_der_r,
										   std::vector<Vector >& _t3_der_r,
										   std::vector<std::vector<Vector > >& _t1_rs,
										   std::vector<std::vector<Vector > >& _t2_rs,
										   std::vector<std::vector<Vector > >& _t3_rs,
										   std::vector<std::vector<Vector > >& _t1_der_rs,
										   std::vector<std::vector<Vector > >& _t2_der_rs,
										   std::vector<std::vector<Vector > >& _t3_der_rs )
	{
		m_nurbs_surface.ComputeSecondVariationOfLocalCSY( parameter_spans[0],
													      parameter_spans[1],
													      parameter_values[0],
													      parameter_values[1],
													      _par_g1,
												          _t1, _t2, _t3,
												          _t1_der, _t2_der, _t3_der,
												          _t1_r, _t2_r, _t3_r,
												          _t1_der_r, _t2_der_r, _t3_der_r,
												          _t1_rs, _t2_rs, _t3_rs,
												          _t1_der_rs, _t2_der_rs, _t3_der_rs );
	}

	// --------------------------------------------------------------------------
	array_1d<int,2> ComputeSurfaceKnotSpans( array_1d<double,2>& parameter_values )
	{
		return m_nurbs_surface.ComputeKnotSpans( parameter_values[0], parameter_values[1] );
	}

	// --------------------------------------------------------------------------
	Matrix ComputeBaseVectors( array_1d<int,2>& parameter_spans, array_1d<double,2>& parameter_values )
	{
		return m_nurbs_surface.ComputeBaseVectors( parameter_spans[0], parameter_spans[1], parameter_values[0], parameter_values[1] );
	}

	// --------------------------------------------------------------------------
	void ComputeGrevilleAbscissae( std::vector<double>& rGrevilleAbscissaeInUDirection, std::vector<double>& rGrevilleAbscissaeInVDirection )
	{
		m_nurbs_surface.ComputeGrevilleAbscissae( rGrevilleAbscissaeInUDirection, rGrevilleAbscissaeInVDirection );
	}

	// --------------------------------------------------------------------------
	void RefineGrevilleAbscissae( std::vector<double>& rGrevilleAbscissaeInUDirection,
								  std::vector<double>& rGrevilleAbscissaeInVDirection,
								  std::vector<double>& rRefinedGrevilleAbscissaeInUDirection,
								  std::vector<double>& rRefinedGrevilleAbscissaeInVDirection )
	{
		m_nurbs_surface.RefineGrevilleAbscissae( rGrevilleAbscissaeInUDirection,
												 rGrevilleAbscissaeInVDirection,
												 rRefinedGrevilleAbscissaeInUDirection,
												 rRefinedGrevilleAbscissaeInVDirection );
	}

	// --------------------------------------------------------------------------
	void EvaluateGradientsForClosestPointSearch( Vector& distance, Matrix& rHessian, Vector& rGradient , array_1d<double,2>& parameter_values )
	{
		return m_nurbs_surface.EvaluateGradientsForClosestPointSearch( distance, rHessian, rGradient, parameter_values[0], parameter_values[1] );
	}

	// --------------------------------------------------------------------------
	void FlagAffectedControlPointsForReconstruction( array_1d<int,2>& parameter_spans, array_1d<double,2>& parameter_values )
	{
		return m_nurbs_surface.FlagAffectedControlPointsForReconstruction( parameter_spans[0], parameter_spans[1], parameter_values[0], parameter_values[1] );
	}

	// --------------------------------------------------------------------------
	unsigned int GetId()
	{
		return m_patch_id;
	}

	// --------------------------------------------------------------------------
	std::vector<ControlPoint>& GetSurfaceControlPoints()
	{
		return m_nurbs_surface.GetControlPoints();
	}

	// --------------------------------------------------------------------------
	std::vector<double>& GetSurfaceKnotVectorU()
	{
		return m_nurbs_surface.GetKnotVectorU();
	}

	// --------------------------------------------------------------------------
	std::vector<double>& GetSurfaceKnotVectorV()
	{
		return m_nurbs_surface.GetKnotVectorV();
	}

	// --------------------------------------------------------------------------
	std::vector<ControlPoint*> GetPointersToAffectedControlPoints( array_1d<int,2>& parameter_spans, array_1d<double,2>& parameter_values )
	{
		return m_nurbs_surface.GetPointersToAffectedControlPoints( parameter_spans[0], parameter_spans[1], parameter_values[0], parameter_values[1] );
	}

	// --------------------------------------------------------------------------
	std::vector<int> GetEquationIdsOfAffectedControlPoints( array_1d<int,2>& parameter_spans, array_1d<double,2>& parameter_values )
	{
		return m_nurbs_surface.GetEquationIdsOfAffectedControlPoints( parameter_spans[0], parameter_spans[1], parameter_values[0], parameter_values[1] );
	}

	// --------------------------------------------------------------------------
	std::vector<BoundaryLoop>& GetBoundaryLoops()
	{
		return m_boundary_loops;
	}

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "Patch";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "Patch";
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
	unsigned int m_patch_id;
	NURBSSurface m_nurbs_surface;
	std::vector<BoundaryLoop> m_boundary_loops;
	std::vector<array_1d<double, 2>> m_closed_boundary_polygon;

	// ==============================================================================
	// General working arrays
	// ==============================================================================
    /// Assignment operator.
    //      Patch& operator=(Patch const& rOther);

    /// Copy constructor.
    //      Patch(Patch const& rOther);

}; // Class Patch

} // namespace Kratos.

#endif // PATCH_H
