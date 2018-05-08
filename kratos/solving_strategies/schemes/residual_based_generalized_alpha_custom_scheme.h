//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Main authors:  Andreas Winterstein
//
//

#if !defined(KRATOS_RESIDUAL_BASED_GENERALIZED_ALPHA_CUSTOM_SCHEME )
#define  KRATOS_RESIDUAL_BASED_GENERALIZED_ALPHA_CUSTOM_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/residual_based_implicit_time_scheme.h"
#include "includes/variables.h"
#include "includes/checks.h"

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

/**
 * @class ResidualBasedGeneralizedAlphaCustomScheme
 * @ingroup KratosCore
 * @brief generalized alpha integration scheme (for dynamic problems)
 * @details This is a dynamic implicit scheme based of the generalized alpha algorithm.
 * @according to the paper of J. Chung and G. M. Hulbert, A Time Integration Algorithms for Structural Dynamics
 * @With Improved Numerical Dissipation: The Generalized-alpha Method
 * @author Andreas Winterstein
 */
template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedGeneralizedAlphaCustomScheme
    : public ResidualBasedImplicitTimeScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedGeneralizedAlphaCustomScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                                  BaseType;

    typedef ResidualBasedImplicitTimeScheme<TSparseSpace,TDenseSpace> ImplicitBaseType;

    typedef typename ImplicitBaseType::TDataType                             TDataType;

    typedef typename ImplicitBaseType::DofsArrayType                     DofsArrayType;

    typedef typename Element::DofsVectorType                            DofsVectorType;

    typedef typename ImplicitBaseType::TSystemMatrixType             TSystemMatrixType;

    typedef typename ImplicitBaseType::TSystemVectorType             TSystemVectorType;

    typedef typename ImplicitBaseType::LocalSystemVectorType     LocalSystemVectorType;

    typedef typename ImplicitBaseType::LocalSystemMatrixType     LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                               NodesArrayType;

    typedef ModelPart::ElementsContainerType                         ElementsArrayType;

    typedef ModelPart::ConditionsContainerType                     ConditionsArrayType;

    typedef typename BaseType::Pointer                                 BaseTypePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @detail The generalized alpha method
     * @rAlpham and rAlphaf The generalized alpha parameters. Default values are 0, which is the Newmark method
     */
    explicit ResidualBasedGeneralizedAlphaCustomScheme(double const rAlphaM = 0.0, double const rAlphaF = 0.0)
        :ImplicitBaseType()
    {

        mGenAlpha.alpha_m = rAlphaM;
        mGenAlpha.alpha_f = rAlphaF;

        // Calcualte generalized alpha coefficients
        mGenAlpha.beta  = 0.25 * (1 - rAlphaM + rAlphaF) * (1 - rAlphaM + rAlphaF);
        mGenAlpha.gamma = 0.5 - rAlphaM + rAlphaF;

        // Allocate auxiliary memory
        const std::size_t num_threads = OpenMPUtils::GetNumThreads();

        mVector.previous_velocity.resize(num_threads);
        mVector.previous_acceleration.resize(num_threads);
        mVector.previous_displacement.resize(num_threads);

        KRATOS_DETAIL("MECHANICAL SCHEME: The Generalized Alpha Time Integration Scheme ")
        << "[alpha_m = " << mGenAlpha.alpha_m
        << " alpha_f = " << mGenAlpha.alpha_f
        << " beta = " << mGenAlpha.beta
        << " gamma = " << mGenAlpha.gamma << "]" <<std::endl;
    }

    /**
     * @brief Copy Constructor.
     */
    ResidualBasedGeneralizedAlphaCustomScheme(ResidualBasedGeneralizedAlphaCustomScheme& rOther)
        :ImplicitBaseType(rOther)
        ,mGenAlpha(rOther.mGenAlpha)
        ,mVector(rOther.mVector)
    {
    }

    /**
     * @brief Clone method
     */
    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new ResidualBasedGeneralizedAlphaCustomScheme(*this) );
    }

    /** Destructor.
     */
    ~ResidualBasedGeneralizedAlphaCustomScheme
    () override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /**
     * @brief Performing the update of the solution
     * @details Incremental update within newton iteration. It updates the state variables at the end of the time step u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param rModelPart The model of the problem to solve
     * @param rDofSet Set of all primary variables
     * @param A LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */

    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;

        // Update of displacement (by DOF)
        const int num_dof = static_cast<int>(rDofSet.size());

        #pragma omp parallel for
        for(int i = 0;  i < num_dof; ++i) {
            auto it_dof = rDofSet.begin() + i;

            if (it_dof->IsFree())
                it_dof->GetSolutionStepValue() += TSparseSpace::GetValue(Dx,it_dof->EquationId());
        }

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>(rModelPart.NumberOfNodes());

        #pragma omp parallel for
        for(int i = 0;  i < num_nodes; ++i) {
            auto it_node = rModelPart.Nodes().begin() + i;

            array_1d<double, 3 > delta_displacement;

            noalias(delta_displacement) = it_node->FastGetSolutionStepValue(DISPLACEMENT) - it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);

            array_1d<double, 3>& current_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& previous_velocity = it_node->FastGetSolutionStepValue(VELOCITY, 1);

            array_1d<double, 3>& current_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION);
            const array_1d<double, 3>& previous_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION, 1);

            UpdateVelocity(current_velocity, delta_displacement, previous_velocity, previous_acceleration);
            UpdateAcceleration(current_acceleration, delta_displacement, previous_velocity, previous_acceleration);
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Performing the prediction of the solution
     * @details It predicts the solution for the current step x = xold + vold * Dt
     * @param rModelPart The model of the problem to solve
     * @param rDofSet set of all primary variables
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */

    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;

        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>(rModelPart.NumberOfNodes());

        array_1d<double, 3 > delta_displacement;

        #pragma omp parallel for private(delta_displacement)
        for(int i = 0;  i < num_nodes; ++i) {
            auto it_node = rModelPart.Nodes().begin() + i;

            //Predicting: NewDisplacement = previous_displacement + previous_velocity * delta_time;
            //ATTENTION::: the prediction is performed only on free nodes

            const array_1d<double, 3 > & previous_acceleration = (it_node)->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double, 3 > & previous_velocity     = (it_node)->FastGetSolutionStepValue(VELOCITY,     1);
            const array_1d<double, 3 > & previous_displacement = (it_node)->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double, 3 > & current_acceleration        = (it_node)->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3 > & current_velocity            = (it_node)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & current_displacement        = (it_node)->FastGetSolutionStepValue(DISPLACEMENT);

            if (it_node -> IsFixed(ACCELERATION_X)) {
                current_displacement[0] = previous_displacement[0] + delta_time * previous_velocity[0] + std::pow(delta_time, 2) * ( 0.5 * (1.0 -  2.0 * mGenAlpha.beta) * previous_acceleration[0] + mGenAlpha.beta * current_acceleration[0]);
            } else if (it_node -> IsFixed(VELOCITY_X)) {
                current_displacement[0] = previous_displacement[0] + 0.5 * delta_time * (previous_velocity[0] + current_velocity[0]) + 0.5 * std::pow(delta_time, 2) * previous_acceleration[0];
            } else if (it_node -> IsFixed(DISPLACEMENT_X) == false) {
                current_displacement[0] = previous_displacement[0] + delta_time * previous_velocity[0] + 0.5 * std::pow(delta_time, 2) * previous_acceleration[0];
            }

            if (it_node -> IsFixed(ACCELERATION_Y)) {
                current_displacement[1] = previous_displacement[1] + delta_time * previous_velocity[1] + std::pow(delta_time, 2) * ( 0.5 * (1.0 -  2.0 * mGenAlpha.beta) * previous_acceleration[1] + mGenAlpha.beta * current_acceleration[1]);
            } else if (it_node -> IsFixed(VELOCITY_Y)) {
                current_displacement[1] = previous_displacement[1] + 0.5 * delta_time * (previous_velocity[1] + current_velocity[1]) + 0.5 * std::pow(delta_time, 2) * previous_acceleration[1] ;
            } else if (it_node -> IsFixed(DISPLACEMENT_Y) == false) {
                current_displacement[1] = previous_displacement[1] + delta_time * previous_velocity[1] + 0.5 * std::pow(delta_time, 2) * previous_acceleration[1];
            }

            // For 3D cases
            if (it_node -> HasDofFor(DISPLACEMENT_Z)) {
                if (it_node -> IsFixed(ACCELERATION_Z)) {
                    current_displacement[2] = previous_displacement[2] + delta_time * previous_velocity[2] + std::pow(delta_time, 2) * ( 0.5 * (1.0 -  2.0 * mGenAlpha.beta) * previous_acceleration[2] + mGenAlpha.beta * current_acceleration[2]);
                } else if (it_node -> IsFixed(VELOCITY_Z)) {
                    current_displacement[2] = previous_displacement[2] + 0.5 * delta_time * (previous_velocity[2] + current_velocity[2]) + 0.5 * std::pow(delta_time, 2) * previous_acceleration[2] ;
                } else if (it_node -> IsFixed(DISPLACEMENT_Z) == false) {
                    current_displacement[2] = previous_displacement[2] + delta_time * previous_velocity[2] + 0.5 * std::pow(delta_time, 2) * previous_acceleration[2];
                }
            }


            // Updating time derivatives ::: Please note that displacements and its time derivatives can not be consistently fixed separately
            noalias(delta_displacement) = current_displacement - previous_displacement;

            UpdateVelocity(current_velocity, delta_displacement, previous_velocity, previous_acceleration);

            UpdateAcceleration(current_acceleration, delta_displacement, previous_velocity, previous_acceleration);
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief It initializes time step solution. Only for reasons if the time step solution is restarted
     * @param rModelPart The model of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;

        ProcessInfo& current_process_info= rModelPart.GetProcessInfo();

        ImplicitBaseType::InitializeSolutionStep(rModelPart, A, Dx, b);

        const double delta_time = current_process_info[DELTA_TIME];

        // Initializing generalized alpha constants
        mGenAlpha.c0 = ( 1.0 / (mGenAlpha.beta * std::pow(delta_time, 2)) );
        mGenAlpha.c1 = ( 1.0 / (mGenAlpha.beta * delta_time) );
        mGenAlpha.c2 = ( 1.0 / (mGenAlpha.beta * 2.0) );
        mGenAlpha.c3 = delta_time;
        KRATOS_CATCH( "" );
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided.
     * @details Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model of the problem to solve
     * @return Zero means  all ok
     */

    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        const int err = ImplicitBaseType::Check(rModelPart);
        if(err != 0) return err;

        // Check for variables keys
        // Verify that the variables are correctly initialized
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
        KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)

        // Check that variables are correctly allocated
        for(auto& rnode : rModelPart.Nodes()) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
        }

        // Check for minimum value of the buffer index
        // Verify buffer size
        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2) << "Insufficient buffer size. Buffer size should be greater than 2. Current size is" << rModelPart.GetBufferSize() << std::endl;

        // Check for admissible value of the AlphaBossak
        KRATOS_ERROR_IF(mGenAlpha.alpha_m > 0.5 || mGenAlpha.alpha_f > 0.5)
        << "Maximum value for alpha_m and alpha_f is less than 0.5.  Current values are: akpha_m = "
        << mGenAlpha.alpha_m << " and alpha_f = " << mGenAlpha.alpha_f << std::endl;

        return 0;
        KRATOS_CATCH( "" );
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

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    struct GeneralizedAlphaMethod
    {
        double alpha_f;
        double alpha_m;

        double beta;
        double gamma;

        double c0, c1, c2, c3;
    };

    /**
     * @brief Vector containing the velocity and acceleration used on integration
     */
    struct GeneralVectors
    {
        std::vector< Vector > previous_displacement;
        std::vector< Vector > previous_velocity;
        std::vector< Vector > previous_acceleration;
    };

    GeneralizedAlphaMethod mGenAlpha; /// The structure containing the Generalized alpha components
    GeneralVectors mVector;        /// The structure containing the velocities and accelerations

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Updating first time Derivative
     * @param CurrentVelocity The current velocity
     * @param DeltaDisplacement The increment of displacement
     * @param PreviousVelocity The previous velocity
     * @param PreviousAcceleration The previous acceleration
     */

    inline void UpdateVelocity(
        array_1d<double, 3>& CurrentVelocity,
        const array_1d<double, 3>& DeltaDisplacement,
        const array_1d<double, 3>& PreviousVelocity,
        const array_1d<double, 3>& PreviousAcceleration
        )
    {

        noalias(CurrentVelocity) = PreviousVelocity + mGenAlpha.c3 * ( (1 - mGenAlpha.gamma)* PreviousAcceleration
        + mGenAlpha.gamma * mGenAlpha.c0 * DeltaDisplacement - mGenAlpha.gamma * mGenAlpha.c1 * PreviousVelocity +
        (mGenAlpha.gamma - mGenAlpha.gamma * mGenAlpha.c2) * PreviousAcceleration );
    }

    /**
     * @brief Updating second time Derivative
     * @param CurrentAcceleration The current velocity
     * @param DeltaDisplacement The increment of displacement
     * @param PreviousVelocity The previous velocity
     * @param PreviousAcceleration The previous acceleration
     */

    inline void UpdateAcceleration(
        array_1d<double, 3>& CurrentAcceleration,
        const array_1d<double, 3>& DeltaDisplacement,
        const array_1d<double, 3>& PreviousVelocity,
        const array_1d<double, 3>& PreviousAcceleration
        )
    {
        noalias(CurrentAcceleration) = mGenAlpha.c0 * DeltaDisplacement
        - mGenAlpha.c1 * PreviousVelocity + (1 - mGenAlpha.c2) * PreviousAcceleration;
    }

    /**
     * @brief It adds the dynamic LHS contribution of the elements
     * M*(c0 - alpha_m*c0) + D*(gamma * c1 - alpha_f * gamma * c1) + K * (1 - alpha_f)
     * @param LHS_Contribution The dynamic contribution for the LHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */

    void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        LocalSystemMatrixType& K,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // Adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
            noalias(LHS_Contribution) += M * (mGenAlpha.c0
            - mGenAlpha.alpha_m * mGenAlpha.c0);


        // Adding damping contribution
        if (D.size1() != 0) // if D matrix declared
            noalias(LHS_Contribution) += D * ( mGenAlpha.gamma * mGenAlpha.c1
            - mGenAlpha.alpha_f * mGenAlpha.gamma * mGenAlpha.c1 );

        //Adding stiffness contribution
        if (K.size1() != 0) // if D matrix declared
            noalias(LHS_Contribution) += K * (1 - mGenAlpha.alpha_f);

    }

    /**
     * @brief It adds the dynamic RHS contribution of the elements b - M*a - D*v
     * @param pElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */

    void AddDynamicsToRHS(
        Element::Pointer pElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        LocalSystemMatrixType& K,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (M.size1() != 0) {

            pElement->GetSecondDerivativesVector(mVector.previous_acceleration[this_thread], 1);
            mVector.previous_acceleration[this_thread] *=
            (mGenAlpha.c2 - mGenAlpha.alpha_m * mGenAlpha.c2 - 1);

            pElement->GetFirstDerivativesVector(mVector.previous_velocity[this_thread], 1);
            mVector.previous_velocity[this_thread] *= (mGenAlpha.c1 - mGenAlpha.alpha_m * mGenAlpha.c1);

            pElement->GetValuesVector(mVector.previous_displacement[this_thread], 1);
            mVector.previous_displacement[this_thread] *= (mGenAlpha.c0 - mGenAlpha.alpha_m * mGenAlpha.c0);

            noalias(RHS_Contribution) += prod(M, mVector.previous_acceleration[this_thread]);
            noalias(RHS_Contribution) += prod(M, mVector.previous_velocity[this_thread]);
            noalias(RHS_Contribution) += prod(M, mVector.previous_displacement[this_thread]);
        }

        // Adding damping contribution
        if (D.size1() != 0) {

            pElement->GetSecondDerivativesVector(mVector.previous_acceleration[this_thread], 1);
            mVector.previous_acceleration[this_thread] *=
            (mGenAlpha.c3 * mGenAlpha.alpha_f
            + mGenAlpha.c3 * mGenAlpha.gamma * mGenAlpha.c2
            - mGenAlpha.c3
            - mGenAlpha.c3 * mGenAlpha.alpha_f * mGenAlpha.gamma * mGenAlpha.c2);

            pElement->GetFirstDerivativesVector(mVector.previous_velocity[this_thread], 1);
            mVector.previous_velocity[this_thread] *= (mGenAlpha.gamma / mGenAlpha.beta
            - mGenAlpha.alpha_f * mGenAlpha.gamma / mGenAlpha.beta
            - 1);

            pElement->GetValuesVector(mVector.previous_displacement[this_thread], 1);
            mVector.previous_displacement[this_thread] *= (mGenAlpha.gamma * mGenAlpha.c1
            - mGenAlpha.alpha_f * mGenAlpha.gamma * mGenAlpha.c1);

            noalias(RHS_Contribution) += prod(D, mVector.previous_acceleration[this_thread]);
            noalias(RHS_Contribution) += prod(D, mVector.previous_velocity[this_thread]);
            noalias(RHS_Contribution) += prod(D, mVector.previous_displacement[this_thread]);
        }
        // Adding stiffness contribution
        if (K.size1() != 0) {
            pElement->GetValuesVector(mVector.previous_displacement[this_thread], 1);
            mVector.previous_displacement[this_thread] *= mGenAlpha.alpha_f;

            noalias(RHS_Contribution) -= prod(K, mVector.previous_displacement[this_thread]);

        }

    }

    bool GetMassMatrixNeeded() override
    {
        return true;
    }

    bool GetDampingMatrixNeeded() override
    {
        return true;
    }

    bool GetStiffnessMatrixNeeded() override
    {
        return true;
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
    ///@{

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

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}
}; /* Class ResidualBasedGeneralizedAlphaCustomScheme */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_GENERALIZED_ALPHA_CUSTOM_SCHEME defined */
