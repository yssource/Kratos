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
//                   Guillermo Casas
//

#ifndef KRATOS_ECA_FLOW_H_INCLUDED
#define KRATOS_ECA_FLOW_H_INCLUDED


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
class EcaFlow : public ManufacturedSolution
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EcaFlow
    KRATOS_CLASS_POINTER_DEFINITION(EcaFlow);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EcaFlow(Properties::Pointer pProperties, Parameters::Pointer pParameters);

    /// Destructor.
    virtual ~EcaFlow() override = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    double Reynolds() override;

    double Strouhal() override;

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
        buffer << "EcaFlow" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "EcaFlow";}

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
    double mB;
    double mVelocity;
    double mLength;
    double mSigma;
    bool mIsPeriodic;
    double mOmega;

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

    // Pressure and derivatives
    virtual double P(array_1d<double, 3>& rCoords, double& rTime) override;
    virtual double DP_DX1(array_1d<double, 3>& rCoords, double& rTime) override;
    virtual double DP_DX2(array_1d<double, 3>& rCoords, double& rTime) override;

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
    EcaFlow& operator=(EcaFlow const& rOther);

    /// Copy constructor.
    EcaFlow(EcaFlow const& rOther);


    ///@}

}; // Class EcaFlow

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                EcaFlow& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const EcaFlow& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ECA_FLOW_H_INCLUDED  defined
