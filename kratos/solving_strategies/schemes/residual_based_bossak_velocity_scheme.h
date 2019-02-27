//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RESIDUAL_BASED_BOSSAK_VELOCITY_SCHEME_H_INCLUDED)
#define KRATOS_RESIDUAL_BASED_BOSSAK_VELOCITY_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "utilities/element_derivatives_extension.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// A scheme for dynamic adjoint equations, using Bossak time integration.
/**
 * It can be used for either first- or second-order time derivatives. Elements
 * and conditions must provide a specialization of AdjointExtensions via their
 * data value container, which allows the scheme to operate independently of
 * the variable arrangements in the element or condition.
 */
template <class TSparseSpace, class TDenseSpace>
class ResidualBasedBossakVelocityScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBossakVelocityScheme);

    typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TSystemMatrixType SystemMatrixType;

    typedef typename BaseType::TSystemVectorType SystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.

    ResidualBasedBossakVelocityScheme(const double NewAlphaBossak)
    {
        mBossak.Alpha = NewAlphaBossak;
        mIsSteady = false;

        if (mIsSteady)
            KRATOS_INFO("ResidualBasedBossakVelocityScheme")
                << "Using steady bossak velocity scheme\n";
        else
            KRATOS_INFO("ResidualBasedBossakVelocityScheme")
                << "Using transient bossak velocity scheme\n";
    }

    ResidualBasedBossakVelocityScheme(Parameters Settings)
    {
        Parameters default_parameters(R"({
            "scheme_type": "bossak",
            "alpha_bossak": -0.3
        })");
        Settings.ValidateAndAssignDefaults(default_parameters);
        mBossak.Alpha = Settings["alpha_bossak"].GetDouble();

        mIsSteady = false;
        if (Settings["scheme_type"].GetString() == "steady")
            mIsSteady = true;

        if (mIsSteady)
            KRATOS_INFO("ResidualBasedBossakVelocityScheme")
                << "Using steady bossak velocity scheme\n";
        else
            KRATOS_INFO("ResidualBasedBossakVelocityScheme")
                << "Using transient bossak velocity scheme\n";
    }

    /// Destructor.
    ~ResidualBasedBossakVelocityScheme() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        BaseType::Initialize(rModelPart);

        // Allocate auxiliary memory.
        int num_threads = OpenMPUtils::GetNumThreads();

        mMassMatrix.resize(num_threads);
        mDampingMatrix.resize(num_threads);
        mForceVector.resize(num_threads);
        mAuxVector.resize(num_threads);

        mValuesVector.resize(num_threads);
        mFirstDerivativeValuesVector.resize(num_threads);
        mSecondDerivativeValuesVector.resize(num_threads);
        mSecondDerivativeValuesVectorOld.resize(num_threads);

        InitializeNodeNeighbourCount(rModelPart.Nodes());

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep(ModelPart& rModelPart,
                                SystemMatrixType& rA,
                                SystemVectorType& rDx,
                                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        const auto& r_current_process_info = rModelPart.GetProcessInfo();
        mBossak = CalculateBossakConstants(mBossak.Alpha, GetTimeStep(r_current_process_info));

        this->CalculateNodeNeighbourCount(rModelPart);

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep(ModelPart& rModelPart,
                              SystemMatrixType& rA,
                              SystemVectorType& rDx,
                              SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        KRATOS_CATCH("");
    }

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                SystemMatrixType& rA,
                SystemVectorType& rDx,
                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        // Update degrees of freedom: adjoint variables associated to the
        // residual of the physical problem.
        this->mpDofUpdater->UpdateDofs(rDofSet, rDx);

        if (mIsSteady)
            return;
        // update the second derivatives
        auto& elements_array = rModelPart.Elements();

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(elements_array.size()); ++i)
        {
            auto it_element = elements_array.begin() + i;
            auto& element_derivatives_dofs_extension =
                *it_element->GetValue(ELEMENT_DERIVATIVES_DOFS_EXTENSION);

            DofsVectorType r_second_derivatives;
            element_derivatives_dofs_extension.GetSecondDerivativesDofList(
                r_second_derivatives, rModelPart.GetProcessInfo());

            Vector r_current_velocity;
            it_element->GetFirstDerivativesVector(r_current_velocity);
            Vector r_old_velocity;
            it_element->GetFirstDerivativesVector(r_old_velocity, 1);

            for (unsigned int j = 0; j < r_second_derivatives.size(); ++j)
            {
                // Skip the dof, if it is represented by shared null ptr
                if (r_second_derivatives[j] == nullptr)
                    continue;

                double& current_value = r_second_derivatives[j]->GetSolutionStepValue();
                const double old_value =
                    r_second_derivatives[j]->GetSolutionStepValue(1);

#pragma omp critical
                current_value =
                    (r_current_velocity[j] - r_old_velocity[j]) * mBossak.C2 -
                    mBossak.C3 * old_value;
            }
        }

        KRATOS_CATCH("");
    }

    // predicts the solution at the current step as
    // v = vold
    void Predict(ModelPart& rModelPart,
                 DofsArrayType& rDofSet,
                 TSystemMatrixType& A,
                 TSystemVectorType& Dv,
                 TSystemVectorType& b) override
    {
        auto& elements_array = rModelPart.Elements();

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(elements_array.size()); ++i)
        {
            auto it_element = elements_array.begin() + i;
            auto& element_derivatives_dofs_extension =
                *it_element->GetValue(ELEMENT_DERIVATIVES_DOFS_EXTENSION);

            DofsVectorType r_first_derivatives;
            element_derivatives_dofs_extension.GetFirstDerivativesDofList(
                r_first_derivatives, rModelPart.GetProcessInfo());
            
            DofsVectorType r_second_derivatives;
            element_derivatives_dofs_extension.GetSecondDerivativesDofList(
                r_second_derivatives, rModelPart.GetProcessInfo());
            
            const std::size_t dof_size = r_first_derivatives.size();

            for (std::size_t i = 0; i < dof_size; ++i)
            {
                double& r_current_velocity = r_first_derivatives[i]->GetSolutionStepValue();
                const double r_old_velocity = r_first_derivatives[i]->GetSolutionStepValue(1);

                if (r_first_derivatives[i]->IsFree())
                {
#pragma omp atomic
                    r_current_velocity = r_old_velocity;           
                }

                if (mIsSteady)
                    continue;

                // Skip the dof, if it is represented by shared null ptr
                if (r_second_derivatives[i] == nullptr)
                    continue;

                const double delta_velocity = r_current_velocity - r_old_velocity;
                double& current_aceleration = r_second_derivatives[i]->GetSolutionStepValue();
#pragma omp critical
                current_aceleration =
                    (r_current_velocity[j] - r_old_velocity[j]) * mBossak.C2 -
                    mBossak.C3 * old_value;                
            }
        }
    }

    void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        auto& r_current_element = *pCurrentElement;
        const auto k = OpenMPUtils::ThisThread();

        r_current_element.InitializeNonLinearIteration(rCurrentProcessInfo);

        r_current_element.CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution,rCurrentProcessInfo);
        r_current_element.CalculateLocalVelocityContribution(mDampingMatrix[k], rRHS_Contribution, rCurrentProcessInfo);

        if (!mIsSteady)
        {
            r_current_element.CalculateMassMatrix(mMassMatrix[k], rCurrentProcessInfo);
           
            AddDynamicsToLHS(rLHS_Contribution, mDampingMatrix[k], mMassMatrix[k], rCurrentProcessInfo);
            AddDynamicsToRHS(r_current_element, rRHS_Contribution, mDampingMatrix[k], mMassMatrix[k], rCurrentProcessInfo);
        }

        this->CalculateResidualLocalContributions(
            r_current_element, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        r_current_element.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Calculate_LHS_Contribution(Element::Pointer pCurrentElement,
                                    LocalSystemMatrixType& rLHS_Contribution,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        LocalSystemVectorType RHS_Contribution;
        CalculateSystemContributions(pCurrentElement, rLHS_Contribution, RHS_Contribution,
                                     rEquationId, rCurrentProcessInfo);
        KRATOS_CATCH("");
    }

    void Condition_CalculateSystemContributions(Condition::Pointer pCurrentCondition,
                                                LocalSystemMatrixType& rLHS_Contribution,
                                                LocalSystemVectorType& rRHS_Contribution,
                                                Condition::EquationIdVectorType& rEquationId,
                                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        // NOT TESTED !!!
        pCurrentCondition->CalculateLocalSystem(
            rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
        KRATOS_CATCH("");
    }

    void Condition_Calculate_LHS_Contribution(Condition::Pointer pCurrentCondition,
                                              LocalSystemMatrixType& rLHS_Contribution,
                                              Condition::EquationIdVectorType& rEquationId,
                                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        LocalSystemVectorType RHS_Contribution;
        Condition_CalculateSystemContributions(pCurrentCondition,
                                               rLHS_Contribution, RHS_Contribution,
                                               rEquationId, rCurrentProcessInfo);
        KRATOS_CATCH("");
    }

    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

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
    std::string Info() const override
    {
        return "ResidualBasedBossakVelocityScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }
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

    //****************************************************************************

    /**
    Kdyn = am*M + D + a1*K
     */
    void AddDynamicsToLHS(LocalSystemMatrixType& LHS_Contribution,
                          LocalSystemMatrixType& D,
                          LocalSystemMatrixType& M,
                          ProcessInfo& CurrentProcessInfo)
    {
        // multipling time scheme factor
        LHS_Contribution *= mBossak.C1;

        // adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += mBossak.C0 * M;
        }

        // adding  damping contribution
        if (D.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += D;
        }
    }

    //****************************************************************************

    /// Add Bossak contributions from the inertial term to the RHS vector.
    /** This essentially performs bdyn = b - M*acc for the current element.
     *  Note that viscous/pressure contributions to the RHS are expected to be added by the element itself.
     *  @param[in] rCurrentElement The fluid element we are assembling.
     *  @param[in/out] rRHS_Contribution The right hand side term where the contribution will be added.
     *  @param[in] rD The elemental velocity/pressure LHS matrix.
     *  @param[in] rM The elemental acceleration LHS matrix.
     *  @param[in] rCurrentProcessInfo ProcessInfo instance for the containing ModelPart.
     */
    void AddDynamicsToRHS(Element::Pointer rCurrentElement,
                          LocalSystemVectorType& rRHS_Contribution,
                          LocalSystemMatrixType& rD,
                          LocalSystemMatrixType& rM,
                          ProcessInfo& rCurrentProcessInfo)
    {
        // adding inertia contribution
        if (rM.size1() != 0)
        {
            int k = OpenMPUtils::ThisThread();
            rCurrentElement->GetSecondDerivativesVector(
                mSecondDerivativeValuesVector[k], 0);
            (mSecondDerivativeValuesVector[k]) *= (1.00 - mBossak.Alpha);
            rCurrentElement->GetSecondDerivativesVector(
                mSecondDerivativeValuesVectorOld[k], 1);
            noalias(mSecondDerivativeValuesVector[k]) +=
                mBossak.Alpha * mSecondDerivativeValuesVectorOld[k];
            noalias(rRHS_Contribution) -= prod(rM, mSecondDerivativeValuesVector[k]);
        }
    }

    /// Add Bossak contributions from the inertial term to the RHS vector.
    /** This essentially performs bdyn = b - M*acc for the current condition.
     *  Note that viscous/pressure contributions to the RHS are expected to be added by the element condition.
     *  @param[in] rCurrentCondition The fluid condition we are assembling.
     *  @param[in/out] rRHS_Contribution The right hand side term where the contribution will be added.
     *  @param[in] rD The elemental velocity/pressure LHS matrix.
     *  @param[in] rM The elemental acceleration LHS matrix.
     *  @param[in] rCurrentProcessInfo ProcessInfo instance for the containing ModelPart.
     */
    void AddDynamicsToRHS(Condition::Pointer rCurrentCondition,
                          LocalSystemVectorType& rRHS_Contribution,
                          LocalSystemMatrixType& D,
                          LocalSystemMatrixType& rM,
                          ProcessInfo& rCurrentProcessInfo)
    {
        // adding inertia contribution
        if (rM.size1() != 0)
        {
            int k = OpenMPUtils::ThisThread();
            rCurrentCondition->GetSecondDerivativesVector(
                mSecondDerivativeValuesVector[k], 0);
            (mSecondDerivativeValuesVector[k]) *= (1.00 - mBossak.Alpha);
            rCurrentCondition->GetSecondDerivativesVector(
                mSecondDerivativeValuesVectorOld[k], 1);
            noalias(mSecondDerivativeValuesVector[k]) +=
                mBossak.Alpha * mSecondDerivativeValuesVectorOld[k];

            noalias(rRHS_Contribution) -= prod(rM, mSecondDerivativeValuesVector[k]);
        }
    }

    void UpdateAcceleration(array_1d<double, 3>& CurrentAcceleration,
                            const array_1d<double, 3>& DeltaVel,
                            const array_1d<double, 3>& OldAcceleration)
    {
        noalias(CurrentAcceleration) = ma0 * DeltaVel + ma2 * OldAcceleration;
    }
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
    struct BossakConstants
    {
        double Alpha;
        double Gamma;
        double Beta;
        double C0;
        double C1;
        double C2;
        double C3;
        double C4;
    };

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    BossakConstants mBossak;
    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater =
        TSparseSpace::CreateDofUpdater();

    std::vector<LocalSystemVectorType> mFirstDerivativeValuesVector;
    std::vector<LocalSystemVectorType> mSecondDerivativeValuesVectorOld;
    std::vector<LocalSystemVectorType> mSecondDerivativeValuesVector;
    std::vector<LocalSystemVectorType> mValuesVector;

    std::vector<LocalSystemMatrixType> mMassMatrix;
    std::vector<LocalSystemMatrixType> mDampingMatrix;
    std::vector<LocalSystemVectorType> mForceVector;
    std::vector<LocalSystemVectorType> mAuxVector;

    bool mIsSteady;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateResidualLocalContributions(Element& rCurrentElement,
                                             LocalSystemMatrixType& rLHS_Contribution,
                                             LocalSystemVectorType& rRHS_Contribution,
                                             ProcessInfo& rCurrentProcessInfo)
    {
        int k = OpenMPUtils::ThisThread();
        auto& r_residual_adjoint = mValuesVector[k];
        rCurrentElement.GetValuesVector(r_residual_adjoint);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, r_residual_adjoint);
    }

    void InitializeNodeNeighbourCount(ModelPart::NodesContainerType& rNodes)
    {
        // This loop should not be omp parallel
        // The operation is not threadsafe if the value is uninitialized
        for (auto& r_node : rNodes)
            r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS, 0.0);
    }

    void CalculateNodeNeighbourCount(ModelPart& rModelPart)
    {
        // Calculate number of neighbour elements for each node.
        const int num_nodes = rModelPart.NumberOfNodes();
#pragma omp parallel for
        for (int i = 0; i < num_nodes; ++i)
        {
            Node<3>& r_node = *(rModelPart.Nodes().begin() + i);
            r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS, 0.0);
        }

        const int num_elements = rModelPart.NumberOfElements();
#pragma omp parallel for
        for (int i = 0; i < num_elements; ++i)
        {
            Element& r_element = *(rModelPart.Elements().begin() + i);
            Geometry<Node<3>>& r_geometry = r_element.GetGeometry();
            for (unsigned j = 0; j < r_geometry.PointsNumber(); ++j)
            {
                double& r_num_neighbour =
                    r_geometry[j].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
#pragma omp atomic
                r_num_neighbour += 1.0;
            }
        }

        rModelPart.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_ELEMENTS);
    }

    void CheckAndResizeThreadStorage(unsigned SystemSize)
    {
        const int k = OpenMPUtils::ThisThread();

        if (mMassMatrix[k].size1() != SystemSize || mMassMatrix[k].size2() != SystemSize)
        {
            mMassMatrix[k].resize(SystemSize, SystemSize, false);
        }

        if (mDampingMatrix[k].size1() != SystemSize || mDampingMatrix[k].size2() != SystemSize)
        {
            mDampingMatrix[k].resize(SystemSize, SystemSize, false);
        }

        if (mForceVector[k].size() != SystemSize)
        {
            mForceVector[k].resize(SystemSize, false);
        }

        if (mAuxVector[k].size() != SystemSize)
        {
            mAuxVector[k].resize(SystemSize, false);
        }

        if (mFirstDerivativeValuesVector[k].size() != SystemSize)
        {
            mFirstDerivativeValuesVector[k].resize(SystemSize, false);
        }

        if (mSecondDerivativeValuesVector[k].size() != SystemSize)
        {
            mSecondDerivativeValuesVector[k].resize(SystemSize, false);
        }

        if (mSecondDerivativeValuesVectorOld[k].size() != SystemSize)
        {
            mSecondDerivativeValuesVectorOld[k].resize(SystemSize, false);
        }
    }

    static BossakConstants CalculateBossakConstants(double Alpha, double DeltaTime)
    {
        BossakConstants bc;
        bc.Alpha = Alpha;
        bc.Gamma = 0.5 - bc.Alpha;
        bc.Beta = 0.25 * pow((1.00 - bc.Alpha), 2);
        bc.C0 = (1.0 - bc.Alpha) / (bc.Gamma * DeltaTime);
        bc.C1 = DeltaTime / (bc.Beta * bc.Gamma);
        bc.C2 = 1.0 / (bc.Gamma * DeltaTime);
        bc.C3 = (1.0 - bc.Gamma) / bc.Gamma;
        return bc;
    }

    static double GetTimeStep(const ProcessInfo& rCurrentProcessInfo)
    {
        const ProcessInfo& r_last_process_info =
            rCurrentProcessInfo.GetPreviousSolutionStepInfo(1);

        double time_step =
            rCurrentProcessInfo.GetValue(TIME) - r_last_process_info.GetValue(TIME);
        KRATOS_ERROR_IF(time_step <= 0.0)
            << "Forwards in time solution is not increasing time from last "
               "step."
            << std::endl;
        return time_step;
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class ResidualBasedBossakVelocityScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BOSSAK_VELOCITY_SCHEME_H_INCLUDED defined */
