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

#ifndef KRATOS_CODINA_VORTEX_H_INCLUDED
#define KRATOS_CODINA_VORTEX_H_INCLUDED


// System includes


// External includes


// Project includes
#include "manufactured_solution.h"


namespace Kratos
{
///@addtogroup ManufacturedFluidSolutionsApplication
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

/// Short class definition.
/** Detail class definition.
*/
class CodinaVortex : public ManufacturedSolution
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CodinaVortex
    KRATOS_CLASS_POINTER_DEFINITION(CodinaVortex);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    CodinaVortex(Properties::Pointer pProperties, Parameters::Pointer pParameters);

    /// Destructor.
    virtual ~CodinaVortex() override = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    double Reynolds() override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "CodinaVortex" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "CodinaVortex";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}

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

    double mA;
    double mVelocity;
    double mLength;
    double mOmega;
    double mDamp;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    virtual bool IsInsideDomain(array_1d<double, 3>& rCoords) override;

    // Velocity components
    virtual double U1(array_1d<double, 3>& rCoords, double& rTime) override;
    virtual double U2(array_1d<double, 3>& rCoords, double& rTime) override;

    // Velocity time derivative
    virtual double DU1_DT(array_1d<double, 3>& rCoords, double& rTime) override;
    virtual double DU2_DT(array_1d<double, 3>& rCoords, double& rTime) override;

    // Velocity first derivative
    virtual double DU1_DX1(array_1d<double, 3>& rCoords, double& rTime) override;
    virtual double DU1_DX2(array_1d<double, 3>& rCoords, double& rTime) override;
    virtual double DU2_DX1(array_1d<double, 3>& rCoords, double& rTime) override;
    virtual double DU2_DX2(array_1d<double, 3>& rCoords, double& rTime) override;

    // Velocity second derivatives
    virtual double DDU1_DX11(array_1d<double, 3>& rCoords, double& rTime) override;
    virtual double DDU1_DX22(array_1d<double, 3>& rCoords, double& rTime) override;
    virtual double DDU2_DX11(array_1d<double, 3>& rCoords, double& rTime) override;
    virtual double DDU2_DX22(array_1d<double, 3>& rCoords, double& rTime) override;

    // Auxiliary functions
    double G(double& rT);
    double DG(double& rT);
    double F(double& rX);
    double DF(double& rX);
    double DDF(double& rX);
    double DDDF(double& rX);

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
    CodinaVortex& operator=(CodinaVortex const& rOther);

    /// Copy constructor.
    CodinaVortex(CodinaVortex const& rOther);


    ///@}

}; // Class CodinaVortex

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                CodinaVortex& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const CodinaVortex& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CODINA_VORTEX_H_INCLUDED  defined
