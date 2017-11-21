// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef CAD_PROJECTION_BASE_H
#define CAD_PROJECTION_BASE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"

// ==============================================================================

namespace Kratos
{
class CADProjectionBase
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef std::vector<Patch> PatchVector; 

    /// Pointer definition of CADProjectionBase
    KRATOS_CLASS_POINTER_DEFINITION(CADProjectionBase);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CADProjectionBase()
    {
    }

    /// Destructor.
    virtual ~CADProjectionBase()
    {
    }

    // --------------------------------------------------------------------------
    virtual void Initialize() = 0;

    // --------------------------------------------------------------------------
    virtual void DetermineNearestCADPoint( NodeType& PointOfInterest, array_1d<double,2>& parameter_values_of_nearest_point, int& patch_index_of_nearest_point ) = 0;

    /// Turn back information as a string.
    virtual std::string Info() const
    {
    return "CADProjectionBase";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
    rOStream << "CADProjectionBase";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

    // ==============================================================================
protected:

    void OptimizeGuessWithNewtonRaphson( NodeType& PointOfInterest,
                                         NodeType& nearest_point,
                                         array_1d<double,2>& parameter_values_of_nearest_point,
                                         Patch& patch_of_nearest_point,
                                         Parameters& projection_parameters )
    {        
        // Initialize what's needed in the Newton-Raphson iteration				
        Vector Distance = ZeroVector(3); 
        Matrix hessian = ZeroMatrix(2,2);
        Vector gradient = ZeroVector(2);
        double determinant_of_hessian = 0;
        Matrix inverse_of_hessian = ZeroMatrix(2,2);
        Point<3> current_nearest_point;
        current_nearest_point[0] = nearest_point.X();
        current_nearest_point[1] = nearest_point.Y();
        current_nearest_point[2] = nearest_point.Z();

        // Variables neeed by the Netwon Raphson algorithm
        double norm_delta_u = 100000000;
        
        int max_iterations = projection_parameters["max_projection_iterations"].GetInt();
        double tolerance = projection_parameters["projection_tolerance"].GetDouble(); 

        // Newton-Raphson algorithm if iterations are specified
        for(int k=0; k<max_iterations; k++)
        {
            // The distance between point on CAD surface point on the FE-mesh
            Distance(0) = current_nearest_point[0] - PointOfInterest.X();
            Distance(1) = current_nearest_point[1] - PointOfInterest.Y();
            Distance(2) = current_nearest_point[2] - PointOfInterest.Z();
            
            // The distance is used to compute hessian and gradient
            patch_of_nearest_point.EvaluateGradientsForClosestPointSearch( Distance, hessian, gradient , parameter_values_of_nearest_point );

            // u_k and v_k are updated
            MathUtils<double>::InvertMatrix( hessian, inverse_of_hessian, determinant_of_hessian );
            Vector delta_u = prod(inverse_of_hessian,gradient);
            parameter_values_of_nearest_point[0] -= delta_u(0);
            parameter_values_of_nearest_point[1] -= delta_u(1);

            // Point on CAD surface is udpated
            patch_of_nearest_point.EvaluateSurfacePoint( parameter_values_of_nearest_point, current_nearest_point );
            
            // Check convergence
            norm_delta_u = norm_2(delta_u);
            if(norm_delta_u<tolerance)
                break;
            else if(k+1==max_iterations)
            {
                std::cout << "WARNING!!! Newton-Raphson in projection did not converge in the following number of iterations: " << k+1 << std::endl;
                KRATOS_WATCH(current_nearest_point)
                KRATOS_WATCH(PointOfInterest)
                std::cout << "Do you want to continue (Y/N)?\n";
                char ans = 'N';
                std::cin >> ans;
                if(ans == 'Y' || ans == 'y')
                    continue;
                else
                    exit(EXIT_SUCCESS);            
            }
        }
    }  

    // ==============================================================================
private:

    /// Assignment operator.
    //      CADProjectionBase& operator=(CADProjectionBase const& rOther);

    /// Copy constructor.
    //      CADProjectionBase(CADProjectionBase const& rOther);

}; // Class CADProjectionBase
} // namespace Kratos.

#endif // CAD_PROJECTION_BASE_H
