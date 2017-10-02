// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef CONTROL_POINT_H_
#define CONTROL_POINT_H_

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ==============================================================================

namespace Kratos
{
class ControlPoint
{
public:

    /// Pointer definition of ControlPoint
    KRATOS_CLASS_POINTER_DEFINITION(ControlPoint);

    /// Default constructor.
    ControlPoint( double X, double Y, double Z, double W, unsigned int global_id)
    {
    	m_coordinates.push_back( X );
    	m_coordinates.push_back( Y );
    	m_coordinates.push_back( Z );
        m_displacements.push_back( 0.0 );
    	m_displacements.push_back( 0.0 );
    	m_displacements.push_back( 0.0 );
		m_w = W;
        m_global_id = global_id;
    }

    /// Destructor.
    virtual ~ControlPoint()
    {
    }
    
    // --------------------------------------------------------------------------
    double GetX()
    {
        return (m_coordinates[0] + m_displacements[0]);
    }

    // --------------------------------------------------------------------------
    double GetY()
    {
     	return (m_coordinates[1] + m_displacements[1]);
    }

    // --------------------------------------------------------------------------
    double GetZ()
    {
    	return (m_coordinates[2] + m_displacements[2]);
    }

    // --------------------------------------------------------------------------
    double GetX0()
    {
    	return m_coordinates[0];
    }

    // --------------------------------------------------------------------------
    double GetY0()
    {
     	return m_coordinates[1];
    }

    // --------------------------------------------------------------------------
    double GetZ0()
    {
    	return m_coordinates[2];
    }    

    // --------------------------------------------------------------------------
    double GetdX()
    {
    	return m_displacements[0];
    }

    // --------------------------------------------------------------------------
    double GetdY()
    {
     	return m_displacements[1];
    }

    // --------------------------------------------------------------------------
    double GetdZ()
    {
    	return m_displacements[2];
    }  

    // --------------------------------------------------------------------------
    void SetdX(double dX)
    {
        m_displacements[0] = dX;
    }  

    // --------------------------------------------------------------------------
    void SetdY(double dY)
    {
        m_displacements[1] = dY;
    }  

    // --------------------------------------------------------------------------
    void SetdZ(double dZ)
    {
        m_displacements[2] = dZ;
    }          

    // --------------------------------------------------------------------------
    double GetWeight()
    {
    	return m_w;
    }

    // --------------------------------------------------------------------------
    int GetGlobalId()
    {
    	return m_global_id;
    }

    // --------------------------------------------------------------------------
    void SetEquationId(int id)
    {
        m_equation_id = id;
    }    

    // --------------------------------------------------------------------------
    int GetEquationId()
    {
        if(m_equation_id<0)
            KRATOS_THROW_ERROR(std::logic_error, "No mapping matrix ID specified for current control point", "");

    	return m_equation_id;
    }

    // --------------------------------------------------------------------------
    void SetRelevantForReconstruction()
    {
    	m_is_relevant_for_reconstruction = true;
    }   

    // --------------------------------------------------------------------------
    bool IsRelevantForReconstruction()
    {
    	return m_is_relevant_for_reconstruction;
    }    

    // ==============================================================================


//    evaluat

    // ==============================================================================


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ControlPoint";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ControlPoint";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


private:
    // ==============================================================================
    // General working arrays
    // ==============================================================================
    std::vector<double> m_coordinates;
    std::vector<double> m_displacements;
    double m_w;
    unsigned int m_global_id;
    int m_equation_id = -1;
    bool m_is_relevant_for_reconstruction = false;


}; // Class ControlPoint

}  // namespace Kratos.

#endif // CONTROL_POINT_H_
