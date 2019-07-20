//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_MANUFACTURED_SOLUTION_H_INCLUDED
#define KRATOS_MANUFACTURED_SOLUTION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
///@addtogroup ManufacturedFluidSolutionApplication
///@{

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

/// Base class for manufactured solutions.
/** The base class for the manufactured solutions.
*/
class ManufacturedSolution
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ManufacturedSolution
    KRATOS_CLASS_POINTER_DEFINITION(ManufacturedSolution);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ManufacturedSolution(Properties::Pointer pProperties, Parameters::Pointer pParameters)
        : mpProperties(pProperties)
        , mpParameters(pParameters)
    {}

    /// Destructor.
    virtual ~ManufacturedSolution(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // Fields

    virtual array_1d<double, 3> BodyForce(array_1d<double, 3>& rCoords, double& rTime);

    virtual array_1d<double, 3> Velocity(array_1d<double, 3>& rCoords, double& rTime);

    virtual double Pressure(array_1d<double, 3>& rCoords, double& rTime);

    // Operators

    virtual array_1d<double, 3> TimeDerivative(array_1d<double, 3>& rCoords, double& rTime);

    virtual array_1d<double, 3> ConvectiveTerm(array_1d<double, 3>& rCoords, double& rTime);

    virtual array_1d<double, 3> ViscousTerm(array_1d<double, 3>& rCoords, double& rTime);

    virtual array_1d<double, 3> PressureGradient(array_1d<double, 3>& rCoords, double& rTime);

    virtual BoundedMatrix<double, 3, 3> VelocityGradient(array_1d<double, 3>& rCoords, double& rTime);

    virtual array_1d<double, 3> VelocityLaplacian(array_1d<double, 3>& rCoords, double& rTime);

    // CFD quantities

    virtual double Reynolds();

    virtual double Strouhal();

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    virtual Parameters& GetParameters();

    virtual Properties& GetProperties();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "ManufacturedSolution" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "ManufacturedSolution";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    Properties::Pointer mpProperties;
    Parameters::Pointer mpParameters;

    double mDensity;
    double mInvDensity;
    double mKinematicViscosity;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    virtual bool IsInsideDomain(array_1d<double, 3>& rCoords) {return false;}

    // Velocity components
    virtual double U1(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double U2(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double U3(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}

    // Velocity time derivative
    virtual double DU1_DT(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DU2_DT(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DU3_DT(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}

    // Velocity first derivative
    virtual double DU1_DX1(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DU1_DX2(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DU1_DX3(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DU2_DX1(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DU2_DX2(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DU2_DX3(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DU3_DX1(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DU3_DX2(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DU3_DX3(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}

    // Velocity second derivatives
    virtual double DDU1_DX11(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DDU1_DX22(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DDU1_DX33(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DDU2_DX11(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DDU2_DX22(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DDU2_DX33(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DDU3_DX11(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DDU3_DX22(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DDU3_DX33(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}

    // Pressure and derivatives
    virtual double P(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DP_DX1(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DP_DX2(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}
    virtual double DP_DX3(array_1d<double, 3>& rCoords, double& rTime) {return 0.0;}

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ManufacturedSolution& operator=(ManufacturedSolution const& rOther);

    /// Copy constructor.
    ManufacturedSolution(ManufacturedSolution const& rOther);

    ///@}

}; // Class ManufacturedSolution

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                ManufacturedSolution& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const ManufacturedSolution& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MANUFACTURED_SOLUTION_H_INCLUDED  defined
