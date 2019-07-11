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

#ifndef KRATOS_MANUFACTURED_SOLUTION_UTILITY_H_INCLUDED
#define KRATOS_MANUFACTURED_SOLUTION_UTILITY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "custom_manufactured/manufactured_solution.h"


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
class KRATOS_API(MANUFACTURED_FLUID_SOLUTIONS_APPLICATION) ManufacturedSolutionUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ManufacturedSolutionUtility
    KRATOS_CLASS_POINTER_DEFINITION(ManufacturedSolutionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ManufacturedSolutionUtility(ModelPart& rModelPart, ManufacturedSolution& rManufactured);

    /// Destructor.
    virtual ~ManufacturedSolutionUtility() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void SetBodyForce();

    void SetVelocity();

    void SetPressure();

    void ComputeExactVelocity();

    void ComputeExactPressure();

    void ComputeExactMaterialAcceleration();

    void ComputeVelocityRelativeError();

    void ComputePressureRelativeError();

    void ComputeMaterialAccelerationError();

    template<class TVarType>
    double ComputeMean(TVarType& rVariable)
    {
        double err = 0;
        #pragma omp parallel for reduction(+:err)
        for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it_node = mrModelPart.NodesBegin() + i;
            err += Norm(it_node->GetValue(rVariable));
        }
        err /= mrModelPart.NumberOfNodes();
        return err;
    }

    template<class TVarType>
    double ComputeRootMeanSquare(TVarType& rVariable)
    {
        double err = 0;
        #pragma omp parallel for reduction(+:err)
        for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it_node = mrModelPart.NodesBegin() + i;
            err += std::pow(Norm(it_node->GetValue(rVariable)), 2);
        }
        err /= mrModelPart.NumberOfNodes();
        return std::sqrt(err);
    }

    template<class TVarType>
    void ComputeError(
        TVarType& rExactVar,
        TVarType& rComputationVar,
        TVarType& rDestinationVar)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it_node = mrModelPart.NodesBegin() + i;
            auto exact = it_node->GetValue(rExactVar);
            auto fem = it_node->FastGetSolutionStepValue(rComputationVar);
            it_node->SetValue(rDestinationVar , exact - fem);
        }
    }

    template<class TVarType>
    void ComputeRelativeError(
        TVarType& rExactVar,
        TVarType& rComputationVar,
        TVarType& rDestinationVar)
    {
        double mean = ComputeMean<TVarType>(rExactVar) + mEpsilon;
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it_node = mrModelPart.NodesBegin() + i;
            auto exact = it_node->GetValue(rExactVar);
            auto fem = it_node->FastGetSolutionStepValue(rComputationVar);
            it_node->SetValue(rDestinationVar , (exact - fem)/mean);
        }
    }

    void RecoverMaterialAcceleration();

    template<class TVarType>
    void BDF1(
        TVarType& rPrimaryVariable,
        TVarType& rDerivativeVariable)
    {
        double dt_inv = 1 / mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it_node = mrModelPart.NodesBegin() + i;
            auto value_n = it_node->FastGetSolutionStepValue(rPrimaryVariable);
            auto value_nn = it_node->FastGetSolutionStepValue(rPrimaryVariable, 1);
            auto derivative = dt_inv * (value_n - value_nn);
            it_node->FastGetSolutionStepValue(rDerivativeVariable) = derivative;
        }
    }

    void ComputeNodalCFL();

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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "ManufacturedSolutionUtility";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "ManufacturedSolutionUtility";}

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

    ModelPart& mrModelPart;
    ManufacturedSolution& mrManufactured;
    const double mEpsilon = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    template<class TVectorType>
    double Norm(TVectorType& rValue) {return norm_2(rValue);}

    double Norm(double& rValue) {return std::abs(rValue);}

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
    ManufacturedSolutionUtility& operator=(ManufacturedSolutionUtility const& rOther);

    /// Copy constructor.
    ManufacturedSolutionUtility(ManufacturedSolutionUtility const& rOther);


    ///@}

}; // Class ManufacturedSolutionUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

inline std::ostream& operator << (std::ostream& rOStream, const ManufacturedSolutionUtility& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MANUFACTURED_SOLUTION_UTILITY_H_INCLUDED  defined
