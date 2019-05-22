//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT)
#define KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT

// System includes

// External includes

// Project includes
#include "includes/element.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utils.h"
#include "rans_modelling_application_variables.h"
#include "stabilized_convection_diffusion_reaction_utilities.h"
#include "utilities/geometrical_sensitivity_utility.h"

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

/// Base class for all Elements.

/**
 * This is the base class for all elements used in KRATOS
 * Elements inherited from this class have to reimplement
 * all public functions that are needed to perform their designated
 * tasks. Due to a dummy implementation of every function though,
 * not all of them have to be implemented if they are not needed for
 * the actual problem
 */
template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionAdjointData>
class StabilizedConvectionDiffusionReactionAdjointElement : public Element
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition of Element

    constexpr static unsigned int TVelPrBlockSize = TDim + 1;

    constexpr static unsigned int TVelPrLocalSize = TNumNodes * TVelPrBlockSize;

    KRATOS_CLASS_POINTER_DEFINITION(StabilizedConvectionDiffusionReactionAdjointElement);

    /// base type: an GeometricalObject that automatically has a unique number
    typedef Element BaseType;

    /// definition of node type (default is: Node<3>)
    typedef Node<3> NodeType;

    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    typedef Properties PropertiesType;

    /// definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;

    /// definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    /// Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef GeometryData GeometryDataType;

    typedef BoundedMatrix<double, TDim, TDim> BoundedMatrixDD;

    typedef BoundedMatrix<double, TNumNodes, TDim> BoundedMatrixND;

    typedef BoundedVector<double, TNumNodes> BoundedVectorN;
    ///@}

    ///@name Life Cycle
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * contructors, copy constructors and destructor: MANDATORY
     */

    /**
     * Constructor.
     */
    StabilizedConvectionDiffusionReactionAdjointElement(IndexType NewId = 0)
        : Element(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    StabilizedConvectionDiffusionReactionAdjointElement(IndexType NewId,
                                                        const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    StabilizedConvectionDiffusionReactionAdjointElement(IndexType NewId,
                                                        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    StabilizedConvectionDiffusionReactionAdjointElement(IndexType NewId,
                                                        GeometryType::Pointer pGeometry,
                                                        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    StabilizedConvectionDiffusionReactionAdjointElement(
        StabilizedConvectionDiffusionReactionAdjointElement const& rOther)
        : Element(rOther)
    {
    }

    ~StabilizedConvectionDiffusionReactionAdjointElement() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * assignment operator: MANDATORY
     */

    /// Assignment operator.

    StabilizedConvectionDiffusionReactionAdjointElement& operator=(
        StabilizedConvectionDiffusionReactionAdjointElement const& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator=(rOther);
        // mpProperties = rOther.mpProperties;
        return *this;
    }

    ///@}
    ///@name Informations
    ///@{

    /** Dimensional space of the element geometry
    @return SizeType, working space dimension of this geometry.
    */

    ///@}
    ///@name Operations
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * Create and Clone methods: MANDATORY
     */

    /**
     * @brief It creates a new element pointer
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        KRATOS_ERROR
            << "Attempting to Create base "
               "StabilizedConvectionDiffusionReactionAdjointElement instances."
            << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * @brief It creates a new element pointer
     * @param NewId the ID of the new element
     * @param pGeom the geometry to be employed
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        KRATOS_ERROR
            << "Attempting to Create base "
               "StabilizedConvectionDiffusionReactionAdjointElement instances."
            << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<StabilizedConvectionDiffusionReactionAdjointElement>(
            NewId, GetGeometry().Create(ThisNodes), pGetProperties());
        KRATOS_CATCH("");
    }

    /**
     * ELEMENTS inherited from this class have to implement next
     * EquationIdVector and GetDofList methods: MANDATORY
     */

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult the elemental equation ID vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "EquationIdVector method. Please implement it in the "
                        "derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        KRATOS_ERROR
            << "Attempting to call base "
               "StabilizedConvectionDiffusionReactionAdjointElement GetDofList "
               "method. Please implement it in the derrived class."
            << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide
     * methods they can be managed internally with a private method to do the
     * same calculations only once: MANDATORY
     */

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rRightHandSideVector the elemental right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        // KRATOS_TRY

        // KRATOS_THROW_ERROR(std::runtime_error,
        //                    "this function is not implemented.", "")

        // KRATOS_CATCH("")
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

        rLeftHandSideMatrix.clear();
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        // KRATOS_TRY

        // KRATOS_THROW_ERROR(std::runtime_error,
        //                    "this function is not implemented.", "")

        // KRATOS_CATCH("")
    }

    /**
     * ELEMENTS inherited from this class must implement this methods
     * if they need to add dynamic element contributions
     * note: first derivatives means the velocities if the displacements are the dof of the analysis
     * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
     * CalculateFirstDerivativesContributions,
     * CalculateFirstDerivativesLHS, CalculateFirstDerivativesRHS methods are : OPTIONAL
     */

    /**
     * this is called during the assembling process in order
     * to calculate the first derivatives contributions for the LHS and RHS
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rRightHandSideVector the elemental right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                                                VectorType& rRightHandSideVector,
                                                ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2())
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix for the first derivatives constributions
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        const Variable<double>& r_derivative_variable = this->GetPrimalVariable();
        CalculateElementTotalResidualScalarDerivatives(
            rLeftHandSideMatrix, r_derivative_variable, rCurrentProcessInfo);
        AddPrimalDampingMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector for the first derivatives constributions
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        if (rRightHandSideVector.size() != 0)
            rRightHandSideVector.resize(0, false);
    }

    /**
     * ELEMENTS inherited from this class must implement this methods
     * if they need to add dynamic element contributions
     * note: second derivatives means the accelerations if the displacements are the dof of the analysis
     * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
     * CalculateSecondDerivativesContributions,
     * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
     */

    /**
     * this is called during the assembling process in order
     * to calculate the second derivative contributions for the LHS and RHS
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rRightHandSideVector the elemental right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                                                 VectorType& rRightHandSideVector,
                                                 ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != 0)
            rLeftHandSideMatrix.resize(0, 0, false);
        if (rRightHandSideVector.size() != 0)
            rRightHandSideVector.resize(0, false);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix for the second derivatives constributions
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

        rLeftHandSideMatrix.clear();

        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];

        AddPrimalMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
        noalias(rLeftHandSideMatrix) = rLeftHandSideMatrix * ((bossak_alpha - 1.0));
        AddPrimalSteadyTermScalarRateDerivatives(rLeftHandSideMatrix, rCurrentProcessInfo);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector for the second derivatives constributions
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        if (rRightHandSideVector.size() != 0)
            rRightHandSideVector.resize(0, false);
    }

    /**
     * ELEMENTS inherited from this class must implement this methods
     * if they need to add dynamic element contributions
     * CalculateMassMatrix and CalculateDampingMatrix methods are: OPTIONAL
     */

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix the elemental mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        if (rMassMatrix.size1() != 0)
            rMassMatrix.resize(0, 0, false);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental damping matrix
     * @param rDampingMatrix the elemental damping matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        if (rDampingMatrix.size1() != 0)
            rDampingMatrix.resize(0, 0, false);
    }

    void Calculate(const Variable<Matrix>& rVariable,
                   Matrix& Output,
                   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rVariable == RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE)
        {
            CalculateElementTotalResidualVelocityDerivatives(Output, rCurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR << "Unsupported variable "
                         << rVariable.Name() << " requested at StabilizedConvectionDiffusionReactionAdjoint::Calculate.";
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates the sensitivity matrix.
     *
     * \f[
     *    \partial_{\mathbf{s}}\mathbf{f}(\mathbf{w}^n)^T
     *  - \partial_{\mathbf{s}}(\mathbf{M}^n \dot{\mathbf{w}}^{n-\alpha})^T
     * \f]
     */
    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rSensitivityVariable == SHAPE_SENSITIVITY)
        {
            this->CalculateElementTotalResidualShapeSensitivity(rOutput, rCurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                         << " not supported." << std::endl;
        }

        KRATOS_CATCH("")
    }

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */

    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(this->Id() < 1) << "StabilizedConvectionDiffusionReacti"
                                           "onAdjointElement found with Id 0 "
                                           "or negative"
                                        << std::endl;

        KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0)
            << "On StabilizedConvectionDiffusionReactionAdjointElement -> "
            << this->Id() << "; Area cannot be less than or equal to 0" << std::endl;

        const Variable<double>& r_primal_variable = this->GetPrimalVariable();
        const Variable<double>& r_primal_relaxed_rate_variable =
            this->GetPrimalVariableRelaxedRate();

        KRATOS_CHECK_VARIABLE_KEY(r_primal_variable);
        KRATOS_CHECK_VARIABLE_KEY(r_primal_relaxed_rate_variable);
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(BOSSAK_ALPHA);
        KRATOS_CHECK_VARIABLE_KEY(NEWMARK_GAMMA);
        KRATOS_CHECK_VARIABLE_KEY(DELTA_TIME);

        for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
        {
            NodeType& r_node = this->GetGeometry()[iNode];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_primal_variable, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_primal_relaxed_rate_variable, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        }

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "StabilizedConvectionDiffusionReactionAdjointElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StabilizedConvectionDiffusionReactionAdjointElement #" << Id();
    }

    /// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
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

    virtual const Variable<double>& GetPrimalVariable() const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base GetPrimalVariable method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        return RANS_ADJOINT_SCALAR_1;

        KRATOS_CATCH("");
    }

    virtual const Variable<double>& GetPrimalVariableRelaxedRate() const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base GetPrimalVariableRelaxedRate method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        return RANS_ADJOINT_SCALAR_1;

        KRATOS_CATCH("");
    }

    virtual double CalculateRelaxedScalarRate(const TConvectionDiffusionReactionAdjointData& rCurrentData,
                                              const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "CalculateRelaxedScalarRate method. "
                        "Please implement it in the derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    virtual void CalculateConvectionDiffusionReactionAdjointData(
        TConvectionDiffusionReactionAdjointData& rData,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;
        KRATOS_ERROR
            << "Attempting to call base "
               "StabilizedConvectionDiffusionReactionAdjointElement "
               "CalculateConvectionDiffusionReactionAdjointData method. "
               "Please implement it in the derrived class."
            << std::endl;
        KRATOS_CATCH("");
    }

    virtual double CalculateEffectiveKinematicViscosity(
        const TConvectionDiffusionReactionAdjointData& rCurrentData,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateEffectiveKinematicViscosity "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    virtual double CalculateReactionTerm(const TConvectionDiffusionReactionAdjointData& rCurrentData,
                                         const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateReactionTerm "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    virtual double CalculateSourceTerm(const TConvectionDiffusionReactionAdjointData& rCurrentData,
                                       const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateSourceTerm "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    virtual void CalculateEffectiveKinematicViscosityScalarDerivatives(
        Vector& rOutput,
        const Variable<double>& rDerivativeVariable,
        const TConvectionDiffusionReactionAdjointData& rCurrentData,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base "
                        "CalculateEffectiveKinematicViscosityScalarDerivatives "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    virtual void CalculateReactionTermDerivatives(Vector& rOutput,
                                                  const Variable<double>& rDerivativeVariable,
                                                  const TConvectionDiffusionReactionAdjointData& rCurrentData,
                                                  const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateReactionTermDerivatives "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    virtual void CalculateSourceTermDerivatives(Vector& rOutput,
                                                const Variable<double>& rDerivativeVariable,
                                                const TConvectionDiffusionReactionAdjointData& rCurrentData,
                                                const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateSourceTermDerivatives "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    virtual void CalculateEffectiveKinematicViscosityVelocityDerivatives(
        Matrix& rOutput,
        const TConvectionDiffusionReactionAdjointData& rCurrentData,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR
            << "Calling base "
               "CalculateEffectiveKinematicViscosityVelocityDerivatives "
               "method in "
               "StabilizedConvectionDiffusionReactionAdjointElement "
               "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    virtual void CalculateReactionTermVelocityDerivatives(
        Matrix& rOutput,
        const TConvectionDiffusionReactionAdjointData& rCurrentData,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateReactionTermVelocityDerivatives "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    virtual void CalculateSourceTermVelocityDerivatives(
        Matrix& rOutput,
        const TConvectionDiffusionReactionAdjointData& rCurrentData,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateSourceTermVelocityDerivatives "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    virtual double CalculateEffectiveKinematicViscosityShapeSensitivity(
        const TConvectionDiffusionReactionAdjointData& rCurrentData,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base "
                        "CalculateEffectiveKinematicViscosityShapeSensitivity "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    virtual double CalculateReactionTermShapeSensitivity(
        const TConvectionDiffusionReactionAdjointData& rCurrentData,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateReactionTermShapeSensitivity "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    virtual double CalculateSourceTermShapeSensitivity(
        const TConvectionDiffusionReactionAdjointData& rCurrentData,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateSourceTermShapeSensitivity "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    void CalculateElementTotalResidualScalarDerivatives(Matrix& rResidualDerivatives,
                                                        const Variable<double>& rDerivativeVariable,
                                                        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rResidualDerivatives.size1() != TNumNodes || rResidualDerivatives.size2() != TNumNodes)
            rResidualDerivatives.resize(TNumNodes, TNumNodes, false);

        rResidualDerivatives.clear();

        AddPrimalSteadyTermScalarDerivatives(
            rResidualDerivatives, rDerivativeVariable, rCurrentProcessInfo);
        AddMassTermScalarDerivatives(rResidualDerivatives, rDerivativeVariable,
                                     rCurrentProcessInfo);
    }

    /// Determine integration point weights and shape funcition derivatives from the element's geometry.
    virtual void CalculateGeometryData(Vector& rGaussWeights,
                                       Matrix& rNContainer,
                                       ShapeFunctionDerivativesArrayType& rDN_DX) const
    {
        const GeometryType& r_geometry = this->GetGeometry();

        RansCalculationUtilities().CalculateGeometryData(
            r_geometry, this->GetIntegrationMethod(), rGaussWeights, rNContainer, rDN_DX);
    }

    ShapeFunctionDerivativesArrayType GetGeometryParameterDerivatives() const
    {
        const GeometryType& r_geometry = this->GetGeometry();
        return RansCalculationUtilities().CalculateGeometryParameterDerivatives(
            r_geometry, this->GetIntegrationMethod());
    }

    double EvaluateInPoint(const Variable<double>& rVariable,
                           const Vector& rShapeFunction,
                           const int Step = 0) const
    {
        return RansCalculationUtilities().EvaluateInPoint(
            this->GetGeometry(), rVariable, rShapeFunction, Step);
    }

    array_1d<double, 3> EvaluateInPoint(const Variable<array_1d<double, 3>>& rVariable,
                                        const Vector& rShapeFunction,
                                        const int Step = 0) const
    {
        return RansCalculationUtilities().EvaluateInPoint(
            this->GetGeometry(), rVariable, rShapeFunction, Step);
    }

    void GetConvectionOperator(BoundedVector<double, TNumNodes>& rOutput,
                               const array_1d<double, 3>& rVector,
                               const Matrix& rShapeDerivatives) const
    {
        rOutput.clear();
        for (unsigned int i = 0; i < TNumNodes; ++i)
            for (unsigned int j = 0; j < TDim; j++)
            {
                rOutput[i] += rVector[j] * rShapeDerivatives(i, j);
            }
    }

    double GetDivergenceOperator(const Variable<array_1d<double, 3>>& rVariable,
                                 const Matrix& rShapeDerivatives,
                                 const int Step = 0) const
    {
        double value = 0.0;
        const GeometryType& r_geometry = this->GetGeometry();

        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            const array_1d<double, 3>& r_value =
                r_geometry[i].FastGetSolutionStepValue(rVariable, Step);
            for (unsigned int j = 0; j < TDim; ++j)
            {
                value += r_value[j] * rShapeDerivatives(i, j);
            }
        }

        return value;
    }

    void CalculateGradient(BoundedMatrix<double, TDim, TDim>& rOutput,
                           const Variable<array_1d<double, 3>>& rVariable,
                           const Matrix& rShapeDerivatives,
                           const int Step = 0) const
    {
        const GeometryType& r_geometry = this->GetGeometry();

        RansCalculationUtilities().CalculateGradient<TDim>(
            rOutput, r_geometry, rVariable, rShapeDerivatives, Step);
    }

    void CalculateGradient(array_1d<double, 3>& rOutput,
                           const Variable<double>& rVariable,
                           const Matrix& rShapeDerivatives,
                           const int Step = 0) const
    {
        const GeometryType& r_geometry = this->GetGeometry();
        RansCalculationUtilities().CalculateGradient(
            rOutput, r_geometry, rVariable, rShapeDerivatives, Step);
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

    void CalculateStabilizationTauScalarDerivatives(Vector& rOutput,
                                                    const double tau,
                                                    const double effective_kinematic_viscosity,
                                                    const double reaction,
                                                    const double element_length,
                                                    const Matrix& rContravariantMetricTensor,
                                                    const Vector& rEffectiveKinematicViscosityScalarDerivatives,
                                                    const Vector& rReactionScalarDerivatives)
    {
        noalias(rOutput) =
            (rEffectiveKinematicViscosityScalarDerivatives *
                 (144 * effective_kinematic_viscosity / std::pow(element_length, 4)) +
             rReactionScalarDerivatives * (reaction)) *
            (-1.0 * std::pow(tau, 3));
    }

    void CalculateStabilizationTauVelocityDerivatives(
        Matrix& rOutput,
        const double tau,
        const double effective_kinematic_viscosity,
        const double reaction,
        const double element_length,
        const array_1d<double, 3>& rVelocity,
        const Matrix& rContravariantMetricTensor,
        const Matrix& rEffectiveKinematicViscosityVelocityDerivatives,
        const Matrix& rReactionVelocityDerivatives,
        const Matrix& rElementLengthDerivatives,
        const Vector& rGaussShapeFunctions)
    {
        Vector contravariant_metric_velocity(TDim);
        const Vector& velocity = RansCalculationUtilities().GetVector<TDim>(rVelocity);

        noalias(contravariant_metric_velocity) =
            prod(rContravariantMetricTensor, velocity) +
            prod(trans(rContravariantMetricTensor), velocity);

        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
            for (std::size_t i_dim = 0; i_dim < TDim; ++i_dim)
                rOutput(i_node, i_dim) = 0.5 * rGaussShapeFunctions[i_node] *
                                         contravariant_metric_velocity[i_dim];

        noalias(rOutput) +=
            rEffectiveKinematicViscosityVelocityDerivatives *
            (144.0 * effective_kinematic_viscosity / std::pow(element_length, 4));
        noalias(rOutput) -= rElementLengthDerivatives *
                            (288.0 * std::pow(effective_kinematic_viscosity, 2) /
                             std::pow(element_length, 5));
        noalias(rOutput) += rReactionVelocityDerivatives * (reaction);
        noalias(rOutput) = rOutput * (-1.0 * std::pow(tau, 3));
    }

    double CalculateStabilizationTauShapeSensitivity(const double tau,
                                                     const double velocity_magnitude,
                                                     const double element_length,
                                                     const double element_length_deriv,
                                                     const double effective_kinematic_viscosity,
                                                     const double effective_kinematic_viscosity_deriv,
                                                     const double reaction,
                                                     const double reaction_deriv)
    {
        double shape_sensitivity = 0.0;

        shape_sensitivity += 4.0 * std::pow(velocity_magnitude, 2) *
                             element_length_deriv / std::pow(element_length, 3);
        shape_sensitivity -= 144.0 * effective_kinematic_viscosity *
                             effective_kinematic_viscosity_deriv /
                             std ::pow(element_length, 4);
        shape_sensitivity += 288.0 * std::pow(effective_kinematic_viscosity, 2) *
                             element_length_deriv / std::pow(element_length, 5);
        shape_sensitivity -= reaction * reaction_deriv;

        shape_sensitivity *= std::pow(tau, 3);

        return shape_sensitivity;
    }

    void CalculateVelocityMagnitudeVelocityDerivative(Matrix& rOutput,
                                                      const double velocity_magnitude,
                                                      const array_1d<double, 3>& rVelocity,
                                                      const Vector& rGaussShapeFunctions)
    {
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
                rOutput(i_node, i_dim) =
                    rVelocity[i_dim] * rGaussShapeFunctions[i_node] / velocity_magnitude;
    }

    void CalculateElementLengthH2VelocityDerivative(Matrix& rOutput,
                                                    const double velocity_magnitude,
                                                    const array_1d<double, 3>& rVelocity,
                                                    const Matrix& rVelocityMagnitudeVelocityDerivatives,
                                                    const Matrix& rContravariantMetricTensor,
                                                    const Vector& rGaussShapeFunctions)
    {
        const Vector& velocity = RansCalculationUtilities().GetVector<TDim>(rVelocity);

        const double sqrt_u_e_u =
            std::sqrt(inner_prod(velocity, prod(rContravariantMetricTensor, velocity)));

        Vector contravariant_metric_velocity(TDim);
        noalias(contravariant_metric_velocity) =
            prod(rContravariantMetricTensor, velocity) +
            prod(trans(rContravariantMetricTensor), velocity);

        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
            for (std::size_t i_dim = 0; i_dim < TDim; ++i_dim)
                rOutput(i_node, i_dim) = rGaussShapeFunctions[i_node] *
                                         contravariant_metric_velocity[i_dim];

        noalias(rOutput) =
            rOutput * (-1.0 * velocity_magnitude / std::pow(sqrt_u_e_u, 3));
        noalias(rOutput) += (rVelocityMagnitudeVelocityDerivatives) * (2.0 / sqrt_u_e_u);
    }

    double CalculateElementLengthH2ShapeSensitivity(const double velocity_magnitude,
                                                    const array_1d<double, 3>& rVelocity,
                                                    const Matrix& rContravariantMetricTensor,
                                                    const Matrix& rContravariantMetricTensorShapeSensitivity)
    {
        const Vector& velocity = RansCalculationUtilities().GetVector<TDim>(rVelocity);

        const double u_e_u = std::pow(
            inner_prod(velocity, prod(rContravariantMetricTensor, velocity)), 1.5);

        return -velocity_magnitude *
               (inner_prod(velocity, prod(rContravariantMetricTensorShapeSensitivity, velocity))) /
               u_e_u;
    }

    void CalculateChiScalarDerivatives(Vector& rOutput,
                                       const double chi,
                                       const double element_length,
                                       const double bossak_alpha,
                                       const double bossak_gamma,
                                       const double delta_time,
                                       const double reaction,
                                       const Vector& rReactionScalarDerivatives)
    {
        const double reaction_tilde =
            reaction + (1 - bossak_alpha) / (bossak_gamma * delta_time);

        CalculateAbsoluteScalarValueScalarDerivatives(
            rOutput, reaction_tilde, rReactionScalarDerivatives);
        noalias(rOutput) = rOutput * (-0.5 * std::pow(chi, 2) * element_length);
    }

    void CalculateChiVelocityDerivatives(Matrix& rOutput,
                                         const double chi,
                                         const double element_length,
                                         const double bossak_alpha,
                                         const double bossak_gamma,
                                         const double delta_time,
                                         const double reaction,
                                         const Matrix& rReactionDerivatives,
                                         const Matrix& rVelocityMagnitudeDerivatives,
                                         const Matrix& rElementLengthDerivatives)
    {
        const double reaction_tilde =
            reaction + (1 - bossak_alpha) / (bossak_gamma * delta_time);
        const double abs_reaction_tilde = std::abs(reaction_tilde);

        CalculateAbsoluteScalarValueVectorDerivatives(rOutput, reaction_tilde,
                                                      rReactionDerivatives);

        noalias(rOutput) = (rOutput * element_length + rElementLengthDerivatives * abs_reaction_tilde +
                            rVelocityMagnitudeDerivatives * 2.0) *
                           (-0.5 * std::pow(chi, 2));
    }

    double CalculateChiShapeSensitivity(const double chi,
                                        const double reaction,
                                        const double reaction_deriv,
                                        const double element_length,
                                        const double element_length_deriv,
                                        const double bossak_alpha,
                                        const double bossak_gamma,
                                        const double delta_time)
    {
        const double reaction_tilde =
            reaction + (1 - bossak_alpha) / (bossak_gamma * delta_time);
        const double abs_reaction_tilde = std::abs(reaction_tilde);

        return -0.5 * std::pow(chi, 2) *
               (abs_reaction_tilde * element_length_deriv +
                reaction_tilde * element_length * reaction_deriv / abs_reaction_tilde);
    }

    void CalculateAbsoluteScalarGradientScalarDerivative(Vector& rOutput,
                                                         const array_1d<double, 3> rScalarGradient,
                                                         const Matrix& rShapeFunctionDerivatives)
    {
        const double scalar_gradient_norm = norm_2(rScalarGradient);

        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
        {
            const Vector& shape_function_gradient = row(rShapeFunctionDerivatives, i_node);
            rOutput[i_node] =
                CalculateScalarProduct(shape_function_gradient, rScalarGradient) / scalar_gradient_norm;
        }
    }

    double CalculateAbsoluteScalarGradientShapeSensitivity(const array_1d<double, 3>& rScalarGradient,
                                                           const Matrix& rShapeFunctionDerivShapeSensitivity,
                                                           const Vector& rNodalScalarValues)
    {
        const double scalar_gradient_norm = norm_2(rScalarGradient);

        Vector scalar_gradient_shape_sensitivity(TDim);
        noalias(scalar_gradient_shape_sensitivity) =
            prod(trans(rShapeFunctionDerivShapeSensitivity), rNodalScalarValues);

        const Vector& scalar_gradient =
            RansCalculationUtilities().GetVector<TDim>(rScalarGradient);
        return inner_prod(scalar_gradient_shape_sensitivity, scalar_gradient) / scalar_gradient_norm;
    }

    void CalculateResidualScalarDerivative(Vector& rOutput,
                                           const double scalar_value,
                                           const double reaction,
                                           const array_1d<double, 3>& rVelocity,
                                           const Vector& rReactionScalarDerivatives,
                                           const Vector& rSourceScalarDerivatives,
                                           const Vector& rShapeFunctions,
                                           const Matrix& rShapeFunctionDerivatives,
                                           const Variable<double>& rDerivativeVariable)
    {
        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
        {
            const Vector& shape_function_gradient = row(rShapeFunctionDerivatives, i_node);
            double value = 0.0;

            value += scalar_value * rReactionScalarDerivatives[i_node];
            value -= rSourceScalarDerivatives[i_node];

            if (this->GetPrimalVariable() == rDerivativeVariable)
            {
                value += reaction * rShapeFunctions[i_node];
                value += CalculateScalarProduct(shape_function_gradient, rVelocity);
            }

            rOutput[i_node] = value;
        }
    }

    void CalculateResidualVelocityDerivative(Matrix& rOutput,
                                             const double primal_variable_value,
                                             const array_1d<double, 3>& rPrimalVariableGradient,
                                             const Matrix& rReactionDerivatives,
                                             const Matrix& rSourceDerivatives,
                                             const Vector& rGaussShapeFunctions)
    {
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
                rOutput(i_node, i_dim) =
                    rGaussShapeFunctions[i_node] * rPrimalVariableGradient[i_dim];

        noalias(rOutput) = rOutput + rReactionDerivatives * primal_variable_value - rSourceDerivatives;
    }

    double CalculateResidualShapeSensitivity(const double residual,
                                             const array_1d<double, 3>& rVelocity,
                                             const Matrix& rShapeFunctionDerivShapeSensitivity,
                                             const double scalar_value,
                                             const Vector& rNodalScalarValues,
                                             const double reaction_deriv,
                                             const double source_deriv)
    {
        const double abs_residual = std::abs(residual);

        const Vector& r_velocity = RansCalculationUtilities().GetVector<TDim>(rVelocity);
        Vector primal_variable_gradient_shape_sensitivity(TDim);
        noalias(primal_variable_gradient_shape_sensitivity) =
            prod(trans(rShapeFunctionDerivShapeSensitivity), rNodalScalarValues);

        return residual *
               (inner_prod(r_velocity, primal_variable_gradient_shape_sensitivity) +
                reaction_deriv * scalar_value - source_deriv) /
               abs_residual;
    }

    void CalculatePositivityPreservationCoefficientScalarDerivatives(
        Vector& rOutput,
        const double chi,
        const double residual,
        const double scalar_gradient_norm,
        const double velocity_norm_square,
        const Vector& rChiScalarDerivatives,
        const Vector& rAbsoluteResidualScalarDerivatives,
        const Vector& rAbsoluteScalarGradientScalarDerivative,
        const Variable<double>& rDerivativeVariable)
    {
        const double abs_residual = std::abs(residual);

        noalias(rOutput) = rAbsoluteResidualScalarDerivatives *
                           (chi / (velocity_norm_square * scalar_gradient_norm));
        noalias(rOutput) +=
            rChiScalarDerivatives *
            (abs_residual / (velocity_norm_square * scalar_gradient_norm));

        if (this->GetPrimalVariable() == rDerivativeVariable)
            noalias(rOutput) -=
                rAbsoluteScalarGradientScalarDerivative *
                (chi * abs_residual / (std::pow(scalar_gradient_norm, 2) * velocity_norm_square));
    }

    void CalculatePositivityPreservationCoefficientVelocityDerivatives(
        Matrix& rOutput,
        const double absolute_residual,
        const double primal_variable_gradient_norm,
        const double velocity_magnitude,
        const double chi,
        const Matrix& rChiDerivatives,
        const Matrix& rAbsoluteResidualDerivatives,
        const Matrix& rVelocityMagnitudeDerivatives)
    {
        const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);

        noalias(rOutput) =
            (rVelocityMagnitudeDerivatives * (-2.0 * chi / velocity_magnitude) + rChiDerivatives) *
                (absolute_residual / (velocity_magnitude_square * primal_variable_gradient_norm)) +
            rAbsoluteResidualDerivatives *
                (chi / (primal_variable_gradient_norm * velocity_magnitude_square));
    }

    double CalculatePositivityPreservationCoefficientShapeSensitivity(
        const double chi,
        const double chi_deriv,
        const double abs_residual,
        const double abs_residual_deriv,
        const double velocity_magnitude_square,
        const double scalar_gradient_norm,
        const double scalar_gradient_norm_deriv)
    {
        return chi_deriv * abs_residual / (scalar_gradient_norm * velocity_magnitude_square) +
               chi * abs_residual_deriv / (scalar_gradient_norm * velocity_magnitude_square) -
               chi * abs_residual * scalar_gradient_norm_deriv /
                   (std::pow(scalar_gradient_norm, 2) * velocity_magnitude_square);
    }

    void CalculatePsiOneScalarDerivatives(Vector& rOutput,
                                          const double velocity_norm,
                                          const double reaction_tilde,
                                          const double tau,
                                          const Vector& rTauScalarDerivatives,
                                          const Vector& rAbsoluteReactionTildeScalarDerivatives)
    {
        const double absolute_reaction_tilde = std::abs(reaction_tilde);

        noalias(rOutput) = rTauScalarDerivatives * (velocity_norm * absolute_reaction_tilde);
        noalias(rOutput) += rAbsoluteReactionTildeScalarDerivatives * (tau * velocity_norm);
    }

    void CalculatePsiOneVelocityDerivatives(Matrix& rOutput,
                                            const double velocity_norm,
                                            const double reaction_tilde,
                                            const double tau,
                                            const Matrix& rTauDerivatives,
                                            const Matrix& rAbsoluteReactionTildeDerivatives,
                                            const Matrix& rVelocityMagnitudeDerivatives)
    {
        noalias(rOutput) = rVelocityMagnitudeDerivatives +
                           rTauDerivatives * (velocity_norm * reaction_tilde) +
                           rVelocityMagnitudeDerivatives * (tau * reaction_tilde) +
                           rAbsoluteReactionTildeDerivatives * (tau * velocity_norm);
    }

    double CalculatePsiOneShapeSensitivity(const double tau,
                                           const double tau_deriv,
                                           const double velocity_magnitude,
                                           const double reaction,
                                           const double reaction_deriv,
                                           const double bossak_alpha,
                                           const double bossak_gamma,
                                           const double delta_time)
    {
        const double reaction_dynamics =
            reaction + (1 - bossak_alpha) / (bossak_gamma * delta_time);
        const double abs_reaction_dynamics = std::abs(reaction_dynamics);

        return tau_deriv * velocity_magnitude * abs_reaction_dynamics +
               tau * velocity_magnitude * reaction_dynamics * reaction_deriv / abs_reaction_dynamics;
    }

    void CalculatePsiTwoScalarDerivatives(Vector& rOutput,
                                          const double element_length,
                                          const double tau,
                                          const double reaction_tilde,
                                          const Vector& rTauScalarDerivatives,
                                          const Vector& rReactionTildeDerivatives,
                                          const Vector& rAbsoluteReactionTildeScalarDerivatives)
    {
        const double absolute_reaction_tilde = std::abs(reaction_tilde);

        noalias(rOutput) = rReactionTildeDerivatives;
        noalias(rOutput) +=
            rTauScalarDerivatives * (reaction_tilde * absolute_reaction_tilde);
        noalias(rOutput) += rReactionTildeDerivatives * (tau * absolute_reaction_tilde);
        noalias(rOutput) += rAbsoluteReactionTildeScalarDerivatives * (tau * reaction_tilde);
        noalias(rOutput) = rOutput * (std::pow(element_length, 2) / 6.0);
    }

    void CalculatePsiTwoVelocityDerivatives(Matrix& rOutput,
                                            const double reaction_tilde,
                                            const double tau,
                                            const double element_length,
                                            const Matrix& rTauDerivatives,
                                            const Matrix& rReactionTildeDerivatives,
                                            const Matrix& rAbsoluteReactionTildeDerivatives,
                                            const Matrix& rElementLengthDerivatives)
    {
        const double abs_reaction_tilde = std::abs(reaction_tilde);

        noalias(rOutput) =
            (rReactionTildeDerivatives + rTauDerivatives * (reaction_tilde * abs_reaction_tilde) +
             rReactionTildeDerivatives * (tau * abs_reaction_tilde) +
             rAbsoluteReactionTildeDerivatives * (tau * reaction_tilde)) *
                std::pow(element_length, 2) / 6.0 +
            rElementLengthDerivatives *
                (element_length *
                 (reaction_tilde + tau * reaction_tilde * abs_reaction_tilde) / 3.0);
    }

    double CalculatePsiTwoShapeSensitivity(const double psi_two,
                                           const double element_length,
                                           const double element_length_deriv,
                                           const double reaction,
                                           const double reaction_deriv,
                                           const double tau,
                                           const double tau_deriv,
                                           const double bossak_alpha,
                                           const double bossak_gamma,
                                           const double delta_time)
    {
        double shape_sensitivity = 0.0;

        const double reaction_dynamics =
            reaction + (1 - bossak_alpha) / (bossak_gamma * delta_time);
        const double abs_reaction_dynamics = std::abs(reaction_dynamics);

        shape_sensitivity += reaction_deriv;
        shape_sensitivity += tau_deriv * reaction_dynamics * abs_reaction_dynamics;
        shape_sensitivity += tau * reaction_deriv * abs_reaction_dynamics;
        shape_sensitivity += tau * reaction_dynamics * reaction_dynamics *
                             reaction_deriv / abs_reaction_dynamics;

        shape_sensitivity *= std::pow(element_length, 2) / 6.0;

        shape_sensitivity += 2.0 * psi_two * element_length_deriv / element_length;

        return shape_sensitivity;
    }

    void CalculateStreamLineDiffusionCoeffScalarDerivatives(
        Vector& rOutput,
        const double element_length,
        const double tau,
        const double velocity_norm,
        const double reaction_tilde,
        const double psi_one,
        const double psi_two,
        const Vector& rPsiOneScalarDerivatives,
        const Vector& rPsiTwoScalarDerivatives,
        const Vector& rTauScalarDerivatives,
        const Vector& rReactionTildeScalarDerivatives,
        const Vector& rEffectiveViscosityScalarDerivatives)
    {
        noalias(rOutput) = rPsiOneScalarDerivatives;
        noalias(rOutput) -= rTauScalarDerivatives * (velocity_norm * reaction_tilde);
        noalias(rOutput) -= rReactionTildeScalarDerivatives * (tau * velocity_norm);

        const double coeff = psi_one - tau * velocity_norm * reaction_tilde;
        noalias(rOutput) = rOutput * (0.5 * element_length * (coeff) / std::abs(coeff));

        noalias(rOutput) += rPsiTwoScalarDerivatives;
        noalias(rOutput) -= rEffectiveViscosityScalarDerivatives;
        noalias(rOutput) -= rTauScalarDerivatives * std::pow(velocity_norm, 2);
    }

    void CalculateStreamLineDiffusionCoeffVelocityDerivatives(
        Matrix& rOutput,
        const double element_length,
        const double tau,
        const double velocity_norm,
        const double reaction_tilde,
        const double psi_one,
        const double psi_two,
        const Matrix& rVelocityMagnitudeDerivatives,
        const Matrix& rPsiOneDerivatives,
        const Matrix& rPsiTwoDerivatives,
        const Matrix& rTauDerivatives,
        const Matrix& rReactionTildeDerivatives,
        const Matrix& rEffectiveViscosityDerivatives,
        const Matrix& rElementLengthDerivatives)
    {
        noalias(rOutput) =
            (rPsiOneDerivatives - rTauDerivatives * (velocity_norm * reaction_tilde) -
             rVelocityMagnitudeDerivatives * (tau * reaction_tilde) -
             rReactionTildeDerivatives * (tau * velocity_norm));

        const double coeff = psi_one - tau * velocity_norm * reaction_tilde;
        noalias(rOutput) = rOutput * (0.5 * element_length * (coeff) / std::abs(coeff));

        noalias(rOutput) += rElementLengthDerivatives * (0.5 * std::abs(coeff));

        noalias(rOutput) += rPsiTwoDerivatives;
        noalias(rOutput) -= rEffectiveViscosityDerivatives;
        noalias(rOutput) -= rTauDerivatives * std::pow(velocity_norm, 2);
        noalias(rOutput) -= rVelocityMagnitudeDerivatives * (2.0 * tau * velocity_norm);
    }

    double CalculateStreamLineDiffusionCoeffShapeSensitivity(
        const double psi_one,
        const double psi_one_deriv,
        const double tau,
        const double tau_deriv,
        const double velocity_magnitude,
        const double reaction,
        const double reaction_deriv,
        const double element_length,
        const double element_length_deriv,
        const double effective_kinematic_viscosity_deriv,
        const double psi_two_deriv,
        const double bossak_alpha,
        const double bossak_gamma,
        const double delta_time)
    {
        const double reaction_dynamics =
            reaction + (1 - bossak_alpha) / (bossak_gamma * delta_time);
        const double coeff = psi_one - tau * velocity_magnitude * reaction_dynamics;
        const double abs_coeff = std::abs(coeff);
        double shape_sensitivity = 0.0;

        shape_sensitivity += psi_one_deriv - tau_deriv * velocity_magnitude * reaction_dynamics -
                             tau * velocity_magnitude * reaction_deriv;
        shape_sensitivity *= 0.5 * coeff * element_length / abs_coeff;
        shape_sensitivity += 0.5 * abs_coeff * element_length_deriv;
        shape_sensitivity -= effective_kinematic_viscosity_deriv +
                             tau_deriv * std::pow(velocity_magnitude, 2);
        shape_sensitivity += psi_two_deriv;

        return shape_sensitivity;
    }

    void CalculateCrossWindDiffusionCoeffScalarDerivatives(
        Vector& rOutput,
        const double psi_one,
        const double element_length,
        const Vector& rPsiOneScalarDerivatives,
        const Vector& rPsiTwoScalarDerivatives,
        const Vector& rEffectiveKinematicViscosityScalarDerivatives)
    {
        noalias(rOutput) = rPsiOneScalarDerivatives *
                           (0.5 * psi_one * element_length / std::abs(psi_one));
        noalias(rOutput) -= rEffectiveKinematicViscosityScalarDerivatives;
        noalias(rOutput) += rPsiTwoScalarDerivatives;
    }

    void CalculateCrossWindDiffusionCoeffVelocityDerivatives(
        Matrix& rOutput,
        const double psi_one,
        const double element_length,
        const Matrix& rPsiOneDerivatives,
        const Matrix& rPsiTwoDerivatives,
        const Matrix& rEffectiveKinematicViscosityDerivatives,
        const Matrix& rElementLengthDerivatives)
    {
        const double abs_psi_one = std::abs(psi_one);

        noalias(rOutput) =
            rPsiOneDerivatives * (0.5 * psi_one * element_length / abs_psi_one) +
            rElementLengthDerivatives * (0.5 * abs_psi_one) -
            rEffectiveKinematicViscosityDerivatives + rPsiTwoDerivatives;
    }

    double CalculateCrossWindDiffusionCoeffShapeSensitivity(const double psi_one,
                                                            const double psi_one_deriv,
                                                            const double element_length,
                                                            const double element_length_deriv,
                                                            const double effective_kinematic_viscosity_deriv,
                                                            const double psi_two_deriv)
    {
        const double abs_psi_one = std::abs(psi_one);
        return 0.5 * psi_one * element_length * psi_one_deriv / abs_psi_one +
               0.5 * abs_psi_one * element_length_deriv -
               effective_kinematic_viscosity_deriv + psi_two_deriv;
    }

    void CalculateAbsoluteScalarValueScalarDerivatives(Vector& rOutput,
                                                       const double scalar_value,
                                                       const Vector& rScalarValueDerivatives)
    {
        noalias(rOutput) =
            rScalarValueDerivatives * (scalar_value / std::abs(scalar_value));
    }

    void CalculateAbsoluteScalarValueVectorDerivatives(Matrix& rOutput,
                                                       const double scalar_value,
                                                       const Matrix& rScalarValueDerivatives)
    {
        noalias(rOutput) =
            rScalarValueDerivatives * (scalar_value / std::abs(scalar_value));
    }

    double CalculateScalarProduct(const Vector& rVector1, const array_1d<double, 3>& rVector2)
    {
        double result = 0.0;
        for (std::size_t i_dim = 0; i_dim < rVector1.size(); ++i_dim)
            result += rVector1[i_dim] * rVector2[i_dim];
        return result;
    }

    void AddPrimalSteadyTermScalarDerivatives(MatrixType& rLeftHandSideMatrix,
                                              const Variable<double>& rDerivativeVariable,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma = rCurrentProcessInfo[NEWMARK_GAMMA];

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        BoundedVector<double, TNumNodes> velocity_convective_terms, scalar_convective_terms;

        Matrix contravariant_metric_tensor(TDim, TDim);

        Vector effective_kinematic_viscosity_derivatives(TNumNodes),
            reaction_derivatives(TNumNodes), source_derivatives(TNumNodes),
            tau_derivatives(TNumNodes), s_derivatives(TNumNodes),
            chi_derivatives(TNumNodes), scalar_gradient_norm_derivative(TNumNodes),
            residual_derivatives(TNumNodes), absolute_residual_derivatives(TNumNodes),
            positivity_preserving_coeff_derivatives(TNumNodes),
            absolute_reaction_tilde_derivatives(TNumNodes),
            psi_one_derivatives(TNumNodes), psi_two_derivatives(TNumNodes),
            streamline_diffusion_coeff_derivatives(TNumNodes),
            crosswind_diffusion_coeff_derivatives(TNumNodes);

        array_1d<double, 3> scalar_gradient;

        TConvectionDiffusionReactionAdjointData current_data;

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            this->CalculateConvectionDiffusionReactionAdjointData(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double scalar_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);

            this->CalculateGradient(scalar_gradient, primal_variable, r_shape_derivatives);
            this->GetConvectionOperator(scalar_convective_terms,
                                        scalar_gradient, r_shape_derivatives);

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);
            this->CalculateEffectiveKinematicViscosityScalarDerivatives(
                effective_kinematic_viscosity_derivatives, rDerivativeVariable,
                current_data, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);
            this->CalculateReactionTermDerivatives(
                reaction_derivatives, rDerivativeVariable, current_data, rCurrentProcessInfo);

            const double source =
                this->CalculateSourceTerm(current_data, rCurrentProcessInfo);
            this->CalculateSourceTermDerivatives(source_derivatives, rDerivativeVariable,
                                                 current_data, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor, reaction,
                effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time);
            this->CalculateStabilizationTauScalarDerivatives(
                tau_derivatives, tau, effective_kinematic_viscosity, reaction,
                element_length, contravariant_metric_tensor,
                effective_kinematic_viscosity_derivatives, reaction_derivatives);

            const double s = std::abs(reaction);
            CalculateAbsoluteScalarValueScalarDerivatives(
                s_derivatives, reaction, reaction_derivatives);

            const double velocity_dot_scalar_gradient =
                inner_prod(velocity, scalar_gradient);

            const double velocity_magnitude = norm_2(velocity);
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
            const double scalar_gradient_norm = norm_2(scalar_gradient);

            double chi, k1, k2, residual, positivity_preserving_coeff;
            if (scalar_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
            {
                residual = this->CalculateRelaxedScalarRate(current_data, rCurrentProcessInfo);
                residual += velocity_dot_scalar_gradient;
                residual += reaction * scalar_value;
                residual -= source;
                residual = std::abs(residual);

                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau, effective_kinematic_viscosity,
                    reaction, bossak_alpha, bossak_gamma, delta_time, element_length);

                positivity_preserving_coeff =
                    residual * chi / (velocity_magnitude_square * scalar_gradient_norm);
            }

            this->CalculateChiScalarDerivatives(
                chi_derivatives, chi, element_length, bossak_alpha,
                bossak_gamma, delta_time, reaction, reaction_derivatives);

            this->CalculateAbsoluteScalarGradientScalarDerivative(
                scalar_gradient_norm_derivative, scalar_gradient, r_shape_derivatives);

            this->CalculateResidualScalarDerivative(
                residual_derivatives, scalar_value, reaction, velocity, reaction_derivatives,
                source_derivatives, gauss_shape_functions, r_shape_derivatives, rDerivativeVariable);

            this->CalculateAbsoluteScalarValueScalarDerivatives(
                absolute_residual_derivatives, residual, residual_derivatives);

            this->CalculatePositivityPreservationCoefficientScalarDerivatives(
                positivity_preserving_coeff_derivatives, chi, residual,
                scalar_gradient_norm, velocity_magnitude_square, chi_derivatives,
                absolute_residual_derivatives, scalar_gradient_norm_derivative, rDerivativeVariable);

            const double reaction_tilde =
                reaction + (1 - bossak_alpha) / (bossak_gamma * delta_time);
            this->CalculateAbsoluteScalarValueScalarDerivatives(
                absolute_reaction_tilde_derivatives, reaction_tilde, reaction_derivatives);

            const double psi_one =
                StabilizedConvectionDiffusionReactionUtilities::CalculatePsiOne(
                    velocity_magnitude, tau, reaction_tilde);
            this->CalculatePsiOneScalarDerivatives(
                psi_one_derivatives, velocity_magnitude, reaction_tilde, tau,
                tau_derivatives, absolute_reaction_tilde_derivatives);

            const double psi_two =
                StabilizedConvectionDiffusionReactionUtilities::CalculatePsiTwo(
                    reaction_tilde, tau, element_length);
            this->CalculatePsiTwoScalarDerivatives(
                psi_two_derivatives, element_length, tau, reaction_tilde, tau_derivatives,
                reaction_derivatives, absolute_reaction_tilde_derivatives);

            this->CalculateStreamLineDiffusionCoeffScalarDerivatives(
                streamline_diffusion_coeff_derivatives, element_length, tau,
                velocity_magnitude, reaction_tilde, psi_one, psi_two,
                psi_one_derivatives, psi_two_derivatives, tau_derivatives,
                reaction_derivatives, effective_kinematic_viscosity_derivatives);

            this->CalculateCrossWindDiffusionCoeffScalarDerivatives(
                crosswind_diffusion_coeff_derivatives, psi_one, element_length, psi_one_derivatives,
                psi_two_derivatives, effective_kinematic_viscosity_derivatives);

            // calculating primal damping matrix scalar derivatives
            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    double dNa_dNc = 0.0;
                    for (unsigned int i = 0; i < TDim; i++)
                        dNa_dNc += r_shape_derivatives(a, i) * r_shape_derivatives(c, i);

                    double value = 0.0;

                    // adding derivative of the diffusion term
                    value += scalar_convective_terms[a] *
                             effective_kinematic_viscosity_derivatives[c];

                    // adding reaction term derivatives
                    value += gauss_shape_functions[a] * reaction_derivatives[c] * scalar_value;

                    // adding SUPG stabilization derivatives
                    value += tau_derivatives[c] *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             velocity_dot_scalar_gradient;
                    value += tau * s_derivatives[c] * gauss_shape_functions[a] *
                             velocity_dot_scalar_gradient;

                    value += tau_derivatives[c] *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             reaction * scalar_value;
                    value += tau * (s_derivatives[c] * gauss_shape_functions[a]) *
                             reaction * scalar_value;
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             reaction_derivatives[c] * scalar_value;

                    // Adding cross wind dissipation derivatives
                    value += positivity_preserving_coeff_derivatives[c] * k2 *
                             scalar_convective_terms[a] * velocity_magnitude_square;
                    value += positivity_preserving_coeff *
                             crosswind_diffusion_coeff_derivatives[c] *
                             scalar_convective_terms[a] * velocity_magnitude_square;

                    value -= positivity_preserving_coeff_derivatives[c] * k2 *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;
                    value -= positivity_preserving_coeff *
                             crosswind_diffusion_coeff_derivatives[c] *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;

                    // Adding stream line dissipation derivatives
                    value += positivity_preserving_coeff_derivatives[c] * k1 *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;
                    value += positivity_preserving_coeff *
                             streamline_diffusion_coeff_derivatives[c] *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;

                    // putting it in the transposed matrix
                    rLeftHandSideMatrix(c, a) += -1.0 * gauss_weights[g] * value;
                }
            }

            // calculating right hand side scalar derivatives
            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    double value = 0.0;

                    value += gauss_shape_functions[a] * source_derivatives[c];
                    value += tau_derivatives[c] *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             source;
                    value += tau * (s_derivatives[c] * gauss_shape_functions[a]) * source;
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             source_derivatives[c];

                    // putting it in the transposed matrix
                    rLeftHandSideMatrix(c, a) += gauss_weights[g] * value;
                }
            }
        }
    }

    void AddPrimalSteadyTermScalarRateDerivatives(MatrixType& rLeftHandSideMatrix,
                                                  ProcessInfo& rCurrentProcessInfo)
    {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma = rCurrentProcessInfo[NEWMARK_GAMMA];

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        BoundedVector<double, TNumNodes> velocity_convective_terms, scalar_convective_terms;

        Matrix contravariant_metric_tensor(TDim, TDim);

        Vector positivity_preserving_coeff_derivatives(TNumNodes);

        array_1d<double, 3> scalar_gradient;

        TConvectionDiffusionReactionAdjointData current_data;

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            this->CalculateConvectionDiffusionReactionAdjointData(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double scalar_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);

            this->CalculateGradient(scalar_gradient, primal_variable, r_shape_derivatives);
            this->GetConvectionOperator(scalar_convective_terms,
                                        scalar_gradient, r_shape_derivatives);

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);

            const double source =
                this->CalculateSourceTerm(current_data, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor, reaction,
                effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time);

            const double velocity_dot_scalar_gradient =
                inner_prod(velocity, scalar_gradient);

            const double velocity_magnitude = norm_2(velocity);
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
            const double scalar_gradient_norm = norm_2(scalar_gradient);

            double chi, k1, k2, residual;
            if (scalar_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
            {
                residual = this->CalculateRelaxedScalarRate(current_data, rCurrentProcessInfo);
                residual += velocity_dot_scalar_gradient;
                residual += reaction * scalar_value;
                residual -= source;

                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau, effective_kinematic_viscosity,
                    reaction, bossak_alpha, bossak_gamma, delta_time, element_length);
            }

            noalias(positivity_preserving_coeff_derivatives) =
                gauss_shape_functions *
                ((1 - bossak_alpha) * (residual / std::abs(residual)) * chi /
                 (scalar_gradient_norm * velocity_magnitude_square));

            // calculating primal damping matrix scalar derivatives
            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    double value = 0.0;

                    // Adding cross wind dissipation derivatives
                    value += positivity_preserving_coeff_derivatives[c] * k2 *
                             scalar_convective_terms[a] * velocity_magnitude_square;

                    value -= positivity_preserving_coeff_derivatives[c] * k2 *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;

                    // Adding stream line dissipation derivatives
                    value += positivity_preserving_coeff_derivatives[c] * k1 *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;

                    // putting it in the transposed matrix
                    rLeftHandSideMatrix(c, a) += -1.0 * gauss_weights[g] * value;
                }
            }
        }
    }

    void AddPrimalDampingMatrix(MatrixType& rPrimalDampingMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();
        const unsigned int num_gauss_points = gauss_weights.size();

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma = rCurrentProcessInfo[NEWMARK_GAMMA];

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        BoundedVector<double, TNumNodes> velocity_convective_terms;

        Matrix contravariant_metric_tensor(TDim, TDim);

        array_1d<double, 3> variable_gradient;

        TConvectionDiffusionReactionAdjointData r_current_data;

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];

            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);
            const double velocity_magnitude = norm_2(velocity);

            this->CalculateConvectionDiffusionReactionAdjointData(
                r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            this->CalculateGradient(variable_gradient, primal_variable, r_shape_derivatives);
            const double variable_gradient_norm = norm_2(variable_gradient);

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(r_current_data, rCurrentProcessInfo);
            const double relaxed_variable_acceleration =
                this->CalculateRelaxedScalarRate(r_current_data, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(r_current_data, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor, reaction,
                effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time);

            // Calculate residual for cross wind dissipation coefficient
            double cross_wind_diffusion{0.0}, stream_line_diffusion{0.0};
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);

            const double velocity_dot_variable_gradient =
                inner_prod(velocity, variable_gradient);
            const double variable_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);

            if (variable_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
            {
                const double source =
                    this->CalculateSourceTerm(r_current_data, rCurrentProcessInfo);

                double residual = relaxed_variable_acceleration;
                residual += velocity_dot_variable_gradient;
                residual += reaction * variable_value;
                residual -= source;
                residual = std::abs(residual);
                residual /= variable_gradient_norm;

                double chi, k1, k2;
                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau, effective_kinematic_viscosity,
                    reaction, bossak_alpha, bossak_gamma, delta_time, element_length);

                stream_line_diffusion = residual * chi * k1 / velocity_magnitude_square;
                cross_wind_diffusion = residual * chi * k2 / velocity_magnitude_square;
            }

            const double s = std::abs(reaction);

            for (unsigned int a = 0; a < TNumNodes; a++)
            {
                for (unsigned int c = 0; c < TNumNodes; c++)
                {
                    double dNa_dNc = 0.0;
                    for (unsigned int i = 0; i < TDim; i++)
                        dNa_dNc += r_shape_derivatives(a, i) * r_shape_derivatives(c, i);

                    double value = 0.0;

                    value += gauss_shape_functions[a] * velocity_convective_terms[c];
                    value += gauss_shape_functions[a] * reaction *
                             gauss_shape_functions[c]; // * positive_values_list[c];
                    value += effective_kinematic_viscosity * dNa_dNc;

                    // Adding SUPG stabilization terms
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             velocity_convective_terms[c];
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             reaction * gauss_shape_functions[c]; // * positive_values_list[c];

                    // Adding cross wind dissipation
                    value += cross_wind_diffusion * dNa_dNc * velocity_magnitude_square;
                    value -= cross_wind_diffusion * velocity_convective_terms[a] *
                             velocity_convective_terms[c];

                    // Adding stream line dissipation
                    value += stream_line_diffusion * velocity_convective_terms[a] *
                             velocity_convective_terms[c];

                    rPrimalDampingMatrix(c, a) += -1.0 * gauss_weights[g] * value;
                }
            }
        }
    }

    void AddLumpedMassMatrix(MatrixType& rMassMatrix, const double Mass)
    {
        for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
            rMassMatrix(iNode, iNode) += Mass;
    }

    void AddPrimalMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();
        const unsigned int num_gauss_points = gauss_weights.size();

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma = rCurrentProcessInfo[NEWMARK_GAMMA];

        TConvectionDiffusionReactionAdjointData current_data;

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
            Matrix contravariant_metric_tensor(r_parameter_derivatives_g.size1(),
                                               r_parameter_derivatives_g.size2());
            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            const double mass = gauss_weights[g] / TNumNodes;
            this->AddLumpedMassMatrix(rMassMatrix, mass);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            BoundedVector<double, TNumNodes> velocity_convective_terms;
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            this->CalculateConvectionDiffusionReactionAdjointData(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor, reaction,
                effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time);

            const double s = std::abs(reaction);

            // Add mass stabilization terms
            for (unsigned int i = 0; i < TNumNodes; ++i)
                for (unsigned int j = 0; j < TNumNodes; ++j)
                    rMassMatrix(j, i) +=
                        gauss_weights[g] * tau *
                        (velocity_convective_terms[i] + s * gauss_shape_functions[i]) *
                        gauss_shape_functions[j];
        }

        KRATOS_CATCH("");
    }

    void AddMassTermScalarDerivatives(MatrixType& rLeftHandSideMatrix,
                                      const Variable<double>& rDerivativeVariable,
                                      const ProcessInfo& rCurrentProcessInfo)
    {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma = rCurrentProcessInfo[NEWMARK_GAMMA];

        Matrix contravariant_metric_tensor(TDim, TDim);

        BoundedVector<double, TNumNodes> velocity_convective_terms;

        TConvectionDiffusionReactionAdjointData current_data;

        Vector effective_kinematic_viscosity_derivatives(TNumNodes),
            reaction_derivatives(TNumNodes), tau_derivatives(TNumNodes),
            s_derivatives(TNumNodes);

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];

            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            this->CalculateConvectionDiffusionReactionAdjointData(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);
            this->CalculateEffectiveKinematicViscosityScalarDerivatives(
                effective_kinematic_viscosity_derivatives, rDerivativeVariable,
                current_data, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);
            this->CalculateReactionTermDerivatives(
                reaction_derivatives, rDerivativeVariable, current_data, rCurrentProcessInfo);

            const double relaxed_scalar_rate =
                this->CalculateRelaxedScalarRate(current_data, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor, reaction,
                effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time);
            this->CalculateStabilizationTauScalarDerivatives(
                tau_derivatives, tau, effective_kinematic_viscosity, reaction,
                element_length, contravariant_metric_tensor,
                effective_kinematic_viscosity_derivatives, reaction_derivatives);

            const double s = std::abs(reaction);
            this->CalculateAbsoluteScalarValueScalarDerivatives(
                s_derivatives, reaction, reaction_derivatives);

            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    double value = 0.0;

                    value += tau_derivatives[c] *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             relaxed_scalar_rate;
                    value += tau * (s_derivatives[c] * gauss_shape_functions[a]) * relaxed_scalar_rate;

                    rLeftHandSideMatrix(c, a) += -1.0 * gauss_weights[g] * value;
                }
            }
        }
    }

    void CalculateElementTotalResidualShapeSensitivity(Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
    {
        constexpr unsigned int TLocalSize = TNumNodes * TDim;
        if (rOutput.size1() != TLocalSize || rOutput.size2() != TNumNodes)
            rOutput.resize(TLocalSize, TNumNodes, false);

        rOutput.clear();

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        // const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        Geometry<NodeType>& r_geometry = this->GetGeometry();

        ShapeParameter deriv;

        Geometry<Point>::JacobiansType J;
        r_geometry.Jacobian(J, this->GetIntegrationMethod());

        Geometry<Point>::ShapeFunctionsGradientsType DN_De;
        DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

        RansCalculationUtilities rans_calculation_utilities;
        RansVariableUtils rans_variable_utils;

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma = rCurrentProcessInfo[NEWMARK_GAMMA];

        const GeometryType::ShapeFunctionsGradientsType& r_dn_de =
            this->GetGeometry().ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

        Matrix contravariant_metric_tensor(TDim, TDim),
            parameter_derivatives_shape_derivs(TDim, TDim),
            contravariant_metric_tensor_deriv(TDim, TDim);

        TConvectionDiffusionReactionAdjointData current_data;

        BoundedVector<double, TNumNodes> velocity_convective_terms,
            primal_variable_gradient_convective_terms,
            convective_primal_variable_gradient_terms_deriv,
            convective_deriv_primal_variable_gradient_terms,
            velocity_convective_terms_deriv,
            primal_variable_gradient_convective_terms_deriv,
            primal_variable_gradient_deriv_convective_terms;

        Vector primal_variable_relaxed_rate_nodal_values(TNumNodes),
            primal_variable_nodal_values(TNumNodes);

        GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_deriv;

        array_1d<double, 3> primal_variable_gradient, primal_variable_gradient_deriv;

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& gauss_r_dn_de = r_dn_de[g];
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);
            const double gauss_weight = gauss_weights[g];

            const Matrix& rJ = J[g];
            const Matrix& rDN_De = DN_De[g];
            const double inv_detJ = 1.0 / MathUtils<double>::DetMat(rJ);
            GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];

            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            this->CalculateConvectionDiffusionReactionAdjointData(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            const double velocity_magnitude = norm_2(velocity);
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            const double primal_variable_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);
            this->CalculateGradient(primal_variable_gradient, primal_variable,
                                    r_shape_derivatives);
            const double primal_variable_gradient_norm = norm_2(primal_variable_gradient);
            const double velocity_dot_primal_variable_gradient =
                inner_prod(velocity, primal_variable_gradient);
            this->GetConvectionOperator(primal_variable_gradient_convective_terms,
                                        primal_variable_gradient, r_shape_derivatives);
            const double primal_variable_relaxed_rate =
                this->CalculateRelaxedScalarRate(current_data, rCurrentProcessInfo);
            rans_variable_utils.GetNodalArray(
                primal_variable_relaxed_rate_nodal_values, *this,
                this->GetPrimalVariableRelaxedRate());

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);
            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);
            const double source =
                this->CalculateSourceTerm(current_data, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor, reaction,
                effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time);

            double chi{0.0}, stream_line_diffusion_coeff{0.0},
                cross_wind_diffusion_coeff{0.0}, residual{0.0};

            if (primal_variable_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
            {
                residual = primal_variable_relaxed_rate;
                residual += velocity_dot_primal_variable_gradient;
                residual += reaction * primal_variable_value;
                residual -= source;
                residual = std::abs(residual);

                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, stream_line_diffusion_coeff, cross_wind_diffusion_coeff,
                    velocity_magnitude, tau, effective_kinematic_viscosity, reaction,
                    bossak_alpha, bossak_gamma, delta_time, element_length);
            }

            rans_variable_utils.GetNodalArray(primal_variable_nodal_values,
                                              *this, primal_variable);

            const double positivity_preserving_coeff =
                chi * residual / (primal_variable_gradient_norm * velocity_magnitude_square);

            const double psi_one =
                StabilizedConvectionDiffusionReactionUtilities::CalculatePsiOne(
                    velocity_magnitude, tau,
                    reaction + (1 - bossak_alpha) / (bossak_gamma * delta_time));

            const double psi_two = StabilizedConvectionDiffusionReactionUtilities::CalculatePsiTwo(
                reaction + (1 - bossak_alpha) / (bossak_gamma * delta_time), tau, element_length);

            const double s = std::abs(reaction);

            for (unsigned int c = 0; c < TNumNodes; ++c)
            {
                const unsigned int block_size = c * TDim;
                for (unsigned int k = 0; k < TDim; ++k)
                {
                    deriv.NodeIndex = c;
                    deriv.Direction = k;

                    double detJ_deriv;
                    geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, DN_DX_deriv);
                    const double gauss_weight_deriv =
                        detJ_deriv * inv_detJ * gauss_weights[g];

                    rans_calculation_utilities.CalculateGeometryParameterDerivativesShapeSensitivity(
                        parameter_derivatives_shape_derivs, deriv,
                        gauss_r_dn_de, r_parameter_derivatives_g);

                    noalias(contravariant_metric_tensor_deriv) =
                        prod(trans(parameter_derivatives_shape_derivs), r_parameter_derivatives_g) +
                        prod(trans(r_parameter_derivatives_g), parameter_derivatives_shape_derivs);

                    this->CalculateGradient(primal_variable_gradient_deriv,
                                            primal_variable, DN_DX_deriv);
                    const double velocity_dot_primal_variable_gradient_deriv =
                        inner_prod(velocity, primal_variable_gradient_deriv);

                    const double reaction_deriv = this->CalculateReactionTermShapeSensitivity(
                        current_data, deriv, detJ_deriv, DN_DX_deriv, rCurrentProcessInfo);
                    const double effective_kinematic_viscosity_deriv =
                        this->CalculateEffectiveKinematicViscosityShapeSensitivity(
                            current_data, deriv, detJ_deriv, DN_DX_deriv, rCurrentProcessInfo);

                    this->GetConvectionOperator(
                        convective_primal_variable_gradient_terms_deriv,
                        primal_variable_gradient_deriv, r_shape_derivatives);
                    this->GetConvectionOperator(convective_deriv_primal_variable_gradient_terms,
                                                primal_variable_gradient, DN_DX_deriv);

                    const double element_length_deriv =
                        this->CalculateElementLengthH2ShapeSensitivity(
                            velocity_magnitude, velocity, contravariant_metric_tensor,
                            contravariant_metric_tensor_deriv);

                    const double tau_deriv = this->CalculateStabilizationTauShapeSensitivity(
                        tau, velocity_magnitude, element_length,
                        element_length_deriv, effective_kinematic_viscosity,
                        effective_kinematic_viscosity_deriv, reaction, reaction_deriv);

                    const double chi_deriv = this->CalculateChiShapeSensitivity(
                        chi, reaction, reaction_deriv, element_length,
                        element_length_deriv, bossak_alpha, bossak_gamma, delta_time);

                    const double primal_variable_gradient_norm_deriv =
                        this->CalculateAbsoluteScalarGradientShapeSensitivity(
                            primal_variable_gradient, DN_DX_deriv, primal_variable_nodal_values);

                    const double source_deriv = this->CalculateSourceTermShapeSensitivity(
                        current_data, deriv, detJ_deriv, DN_DX_deriv, rCurrentProcessInfo);

                    const double residual_deriv = this->CalculateResidualShapeSensitivity(
                        residual, velocity, DN_DX_deriv, primal_variable_value,
                        primal_variable_nodal_values, reaction_deriv, source_deriv);

                    const double positivity_preserving_coeff_deriv =
                        CalculatePositivityPreservationCoefficientShapeSensitivity(
                            chi, chi_deriv, residual, residual_deriv,
                            velocity_magnitude_square, primal_variable_gradient_norm,
                            primal_variable_gradient_norm_deriv);

                    const double psi_one_deriv = CalculatePsiOneShapeSensitivity(
                        tau, tau_deriv, velocity_magnitude, reaction,
                        reaction_deriv, bossak_alpha, bossak_gamma, delta_time);

                    const double psi_two_deriv = CalculatePsiTwoShapeSensitivity(
                        psi_two, element_length, element_length_deriv, reaction, reaction_deriv,
                        tau, tau_deriv, bossak_alpha, bossak_gamma, delta_time);

                    const double stream_line_diffusion_coeff_deriv =
                        CalculateStreamLineDiffusionCoeffShapeSensitivity(
                            psi_one, psi_one_deriv, tau, tau_deriv, velocity_magnitude,
                            reaction, reaction_deriv, element_length,
                            element_length_deriv, effective_kinematic_viscosity_deriv,
                            psi_two_deriv, bossak_alpha, bossak_gamma, delta_time);

                    const double cross_wind_diffusion_coeff_deriv =
                        CalculateCrossWindDiffusionCoeffShapeSensitivity(
                            psi_one, psi_one_deriv, element_length, element_length_deriv,
                            effective_kinematic_viscosity_deriv, psi_two_deriv);

                    const double s_deriv = reaction * reaction_deriv / s;

                    this->GetConvectionOperator(velocity_convective_terms_deriv,
                                                velocity, DN_DX_deriv);
                    this->GetConvectionOperator(primal_variable_gradient_convective_terms_deriv,
                                                primal_variable_gradient, DN_DX_deriv);
                    this->GetConvectionOperator(
                        primal_variable_gradient_deriv_convective_terms,
                        primal_variable_gradient_deriv, r_shape_derivatives);

                    for (unsigned int a = 0; a < TNumNodes; ++a)
                    {
                        const double tau_operator = velocity_convective_terms[a] +
                                                    s * gauss_shape_functions[a];
                        const double tau_operator_deriv =
                            velocity_convective_terms_deriv[a] +
                            s_deriv * gauss_shape_functions[a];

                        double value = 0.0;

                        // adding convective term shape sensitivity
                        value += gauss_shape_functions[a] *
                                 velocity_dot_primal_variable_gradient_deriv * gauss_weight;
                        value += gauss_shape_functions[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight_deriv;

                        // adding reaction term shape sensitivity
                        value += gauss_shape_functions[a] * primal_variable_value *
                                 reaction_deriv * gauss_weight;
                        value += gauss_shape_functions[a] * primal_variable_value *
                                 reaction * gauss_weight_deriv;

                        // adding diffusion term shape sensitivity
                        value += effective_kinematic_viscosity_deriv *
                                 primal_variable_gradient_convective_terms[a] * gauss_weight;
                        value += effective_kinematic_viscosity *
                                 convective_deriv_primal_variable_gradient_terms[a] *
                                 gauss_weight;
                        value += effective_kinematic_viscosity *
                                 convective_primal_variable_gradient_terms_deriv[a] *
                                 gauss_weight;
                        value += effective_kinematic_viscosity *
                                 primal_variable_gradient_convective_terms[a] *
                                 gauss_weight_deriv;

                        // adding SUPG term derivatives
                        value += tau_deriv * tau_operator *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value += tau * tau_operator_deriv *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value += tau * tau_operator *
                                 velocity_dot_primal_variable_gradient_deriv * gauss_weight;
                        value += tau * tau_operator *
                                 velocity_dot_primal_variable_gradient * gauss_weight_deriv;

                        value += tau_deriv * tau_operator * reaction *
                                 primal_variable_value * gauss_weight;
                        value += tau * tau_operator_deriv * reaction *
                                 primal_variable_value * gauss_weight;
                        value += tau * tau_operator * reaction_deriv *
                                 primal_variable_value * gauss_weight;
                        value += tau * tau_operator * reaction *
                                 primal_variable_value * gauss_weight_deriv;

                        // adding cross wind dissipation term derivatives
                        value += velocity_magnitude_square * positivity_preserving_coeff_deriv *
                                 cross_wind_diffusion_coeff *
                                 primal_variable_gradient_convective_terms[a] * gauss_weight;
                        value += velocity_magnitude_square * positivity_preserving_coeff *
                                 cross_wind_diffusion_coeff_deriv *
                                 primal_variable_gradient_convective_terms[a] * gauss_weight;
                        value += velocity_magnitude_square * positivity_preserving_coeff *
                                 cross_wind_diffusion_coeff *
                                 primal_variable_gradient_convective_terms_deriv[a] *
                                 gauss_weight;
                        value += velocity_magnitude_square * positivity_preserving_coeff *
                                 cross_wind_diffusion_coeff *
                                 primal_variable_gradient_deriv_convective_terms[a] *
                                 gauss_weight;
                        value += velocity_magnitude_square * positivity_preserving_coeff *
                                 cross_wind_diffusion_coeff *
                                 primal_variable_gradient_convective_terms[a] *
                                 gauss_weight_deriv;

                        value -= positivity_preserving_coeff_deriv * cross_wind_diffusion_coeff *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value -= positivity_preserving_coeff * cross_wind_diffusion_coeff_deriv *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value -= positivity_preserving_coeff * cross_wind_diffusion_coeff *
                                 velocity_convective_terms_deriv[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value -= positivity_preserving_coeff * cross_wind_diffusion_coeff *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient_deriv * gauss_weight;
                        value -= positivity_preserving_coeff * cross_wind_diffusion_coeff *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight_deriv;

                        // adding stream line diffusion term derivatives
                        value += positivity_preserving_coeff_deriv * stream_line_diffusion_coeff *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value += positivity_preserving_coeff * stream_line_diffusion_coeff_deriv *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value += positivity_preserving_coeff * stream_line_diffusion_coeff *
                                 velocity_convective_terms_deriv[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value += positivity_preserving_coeff * stream_line_diffusion_coeff *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient_deriv * gauss_weight;
                        value += positivity_preserving_coeff * stream_line_diffusion_coeff *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight_deriv;

                        // adding right hand side terms shape gradient terms

                        value -= gauss_shape_functions[a] * source_deriv * gauss_weight;
                        value -= gauss_shape_functions[a] * source * gauss_weight_deriv;

                        value -= tau_deriv * tau_operator * source * gauss_weight;
                        value -= tau * tau_operator_deriv * source * gauss_weight;
                        value -= tau * tau_operator * source_deriv * gauss_weight;
                        value -= tau * tau_operator * source * gauss_weight_deriv;

                        // adding mass shape gradient terms
                        value += gauss_weight_deriv *
                                 primal_variable_relaxed_rate_nodal_values[a] / TNumNodes;
                        value += tau_deriv * tau_operator *
                                 primal_variable_relaxed_rate * gauss_weight;
                        value += tau * tau_operator_deriv *
                                 primal_variable_relaxed_rate * gauss_weight;
                        value += tau * tau_operator * primal_variable_relaxed_rate * gauss_weight_deriv;

                        rOutput(block_size + k, a) -= value;
                    }
                }
            }
        }
    }

    void CalculateElementTotalResidualVelocityDerivatives(Matrix& rResidualDerivatives,
                                                          const ProcessInfo& rCurrentProcessInfo)
    {
        if (rResidualDerivatives.size1() != TVelPrLocalSize ||
            rResidualDerivatives.size2() != TNumNodes)
            rResidualDerivatives.resize(TVelPrLocalSize, TNumNodes, false);

        rResidualDerivatives.clear();

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma = rCurrentProcessInfo[NEWMARK_GAMMA];

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        BoundedVector<double, TNumNodes> primal_variable_gradient_convective_terms,
            velocity_convective_terms;

        Matrix effective_kinematic_viscosity_derivatives(TNumNodes, TDim),
            reaction_derivatives(TNumNodes, TDim),
            velocity_magnitude_derivatives(TNumNodes, TDim),
            element_length_derivatives(TNumNodes, TDim),
            tau_derivatives(TNumNodes, TDim), source_derivatives(TNumNodes, TDim),
            chi_derivatives(TNumNodes, TDim), residual_derivatives(TNumNodes, TDim),
            absolute_residual_derivatives(TNumNodes, TDim),
            positivity_preservation_coeff_derivatives(TNumNodes, TDim),
            absolute_reaction_tilde_derivatives(TNumNodes, TDim),
            psi_one_derivatives(TNumNodes, TDim), psi_two_derivatives(TNumNodes, TDim),
            stream_line_diffusion_coeff_derivatives(TNumNodes, TDim),
            cross_wind_diffusion_coeff_derivatives(TNumNodes, TDim),
            s_derivatives(TNumNodes, TDim), contravariant_metric_tensor(TDim, TDim);

        array_1d<double, 3> primal_variable_gradient;

        TConvectionDiffusionReactionAdjointData current_data;

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            this->CalculateConvectionDiffusionReactionAdjointData(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double primal_variable_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);
            this->CalculateGradient(primal_variable_gradient, primal_variable,
                                    r_shape_derivatives);
            this->GetConvectionOperator(primal_variable_gradient_convective_terms,
                                        primal_variable_gradient, r_shape_derivatives);

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);
            this->CalculateEffectiveKinematicViscosityVelocityDerivatives(
                effective_kinematic_viscosity_derivatives, current_data, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);
            this->CalculateReactionTermVelocityDerivatives(
                reaction_derivatives, current_data, rCurrentProcessInfo);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            const double velocity_magnitude = norm_2(velocity);
            this->CalculateVelocityMagnitudeVelocityDerivative(
                velocity_magnitude_derivatives, velocity_magnitude, velocity,
                gauss_shape_functions);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor, reaction,
                effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time);

            this->CalculateElementLengthH2VelocityDerivative(
                element_length_derivatives, velocity_magnitude, velocity,
                velocity_magnitude_derivatives, contravariant_metric_tensor,
                gauss_shape_functions);

            this->CalculateStabilizationTauVelocityDerivatives(
                tau_derivatives, tau, effective_kinematic_viscosity, reaction,
                element_length, velocity, contravariant_metric_tensor,
                effective_kinematic_viscosity_derivatives, reaction_derivatives,
                element_length_derivatives, gauss_shape_functions);

            const double source =
                this->CalculateSourceTerm(current_data, rCurrentProcessInfo);
            this->CalculateSourceTermVelocityDerivatives(
                source_derivatives, current_data, rCurrentProcessInfo);

            const double velocity_dot_primal_variable_gradient =
                inner_prod(velocity, primal_variable_gradient);

            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
            const double primal_variable_gradient_norm = norm_2(primal_variable_gradient);

            double chi, k1, k2;

            const double primal_variable_relaxed_rate =
                this->CalculateRelaxedScalarRate(current_data, rCurrentProcessInfo);

            double residual = primal_variable_relaxed_rate;
            residual += velocity_dot_primal_variable_gradient;
            residual += reaction * primal_variable_value;
            residual -= source;

            StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                chi, k1, k2, velocity_magnitude, tau, effective_kinematic_viscosity,
                reaction, bossak_alpha, bossak_gamma, delta_time, element_length);

            const double positivity_preservation_coeff =
                std::abs(residual) * chi /
                (velocity_magnitude_square * primal_variable_gradient_norm);

            this->CalculateChiVelocityDerivatives(
                chi_derivatives, chi, element_length, bossak_alpha,
                bossak_gamma, delta_time, reaction, reaction_derivatives,
                velocity_magnitude_derivatives, element_length_derivatives);

            this->CalculateResidualVelocityDerivative(
                residual_derivatives, primal_variable_value, primal_variable_gradient,
                reaction_derivatives, source_derivatives, gauss_shape_functions);

            this->CalculateAbsoluteScalarValueVectorDerivatives(
                absolute_residual_derivatives, residual, residual_derivatives);

            this->CalculatePositivityPreservationCoefficientVelocityDerivatives(
                positivity_preservation_coeff_derivatives, std::abs(residual),
                primal_variable_gradient_norm, velocity_magnitude, chi, chi_derivatives,
                absolute_residual_derivatives, velocity_magnitude_derivatives);

            const double reaction_tilde =
                reaction + (1 - bossak_alpha) / (bossak_gamma * delta_time);
            this->CalculateAbsoluteScalarValueVectorDerivatives(
                absolute_reaction_tilde_derivatives, reaction_tilde, reaction_derivatives);

            const double psi_one =
                StabilizedConvectionDiffusionReactionUtilities::CalculatePsiOne(
                    velocity_magnitude, tau, reaction_tilde);
            this->CalculatePsiOneVelocityDerivatives(
                psi_one_derivatives, velocity_magnitude, reaction_tilde, tau, tau_derivatives,
                absolute_reaction_tilde_derivatives, velocity_magnitude_derivatives);

            const double psi_two =
                StabilizedConvectionDiffusionReactionUtilities::CalculatePsiTwo(
                    reaction_tilde, tau, element_length);
            this->CalculatePsiTwoVelocityDerivatives(
                psi_two_derivatives, reaction_tilde, tau, element_length,
                tau_derivatives, reaction_derivatives,
                absolute_reaction_tilde_derivatives, element_length_derivatives);

            this->CalculateStreamLineDiffusionCoeffVelocityDerivatives(
                stream_line_diffusion_coeff_derivatives, element_length, tau,
                velocity_magnitude, reaction_tilde, psi_one, psi_two,
                velocity_magnitude_derivatives, psi_one_derivatives,
                psi_two_derivatives, tau_derivatives, reaction_derivatives,
                effective_kinematic_viscosity_derivatives, element_length_derivatives);

            this->CalculateCrossWindDiffusionCoeffVelocityDerivatives(
                cross_wind_diffusion_coeff_derivatives, psi_one, element_length,
                psi_one_derivatives, psi_two_derivatives,
                effective_kinematic_viscosity_derivatives, element_length_derivatives);

            const double s = std::abs(reaction);
            this->CalculateAbsoluteScalarValueVectorDerivatives(
                s_derivatives, reaction, reaction_derivatives);

            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    unsigned int block_size = c * TVelPrBlockSize;
                    for (unsigned int k = 0; k < TDim; ++k)
                    {
                        // adding damping matrix derivatives
                        double value = 0.0;

                        // convection term derivatives
                        value += gauss_shape_functions[a] * gauss_shape_functions[c] *
                                 primal_variable_gradient[k];

                        // adding diffusion term derivatives
                        value += effective_kinematic_viscosity_derivatives(c, k) *
                                 primal_variable_gradient_convective_terms[a];

                        // adding reaction term derivatives
                        value += reaction_derivatives(c, k) *
                                 gauss_shape_functions[a] * primal_variable_value;

                        // adding SUPG term derivatives
                        value += tau_derivatives(c, k) *
                                 (velocity_convective_terms[a] +
                                  s * gauss_shape_functions[a]) *
                                 velocity_dot_primal_variable_gradient;
                        value += tau *
                                 (gauss_shape_functions[c] * r_shape_derivatives(a, k) +
                                  s_derivatives(c, k) * gauss_shape_functions[a]) *
                                 velocity_dot_primal_variable_gradient;
                        value +=
                            tau *
                            (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                            gauss_shape_functions[c] * primal_variable_gradient[k];

                        value += tau_derivatives(c, k) *
                                 (velocity_convective_terms[a] +
                                  s * gauss_shape_functions[a]) *
                                 reaction * primal_variable_value;
                        value += tau *
                                 (gauss_shape_functions[c] * r_shape_derivatives(a, k) +
                                  s_derivatives(c, k) * gauss_shape_functions[a]) *
                                 reaction * primal_variable_value;
                        value += tau *
                                 (velocity_convective_terms[a] +
                                  s * gauss_shape_functions[a]) *
                                 reaction_derivatives(c, k) * primal_variable_value;

                        // adding cross wind dissipation derivatives
                        value += positivity_preservation_coeff_derivatives(c, k) *
                                 k2 * velocity_magnitude_square *
                                 primal_variable_gradient_convective_terms[a];
                        value += positivity_preservation_coeff *
                                 cross_wind_diffusion_coeff_derivatives(c, k) *
                                 velocity_magnitude_square *
                                 primal_variable_gradient_convective_terms[a];
                        value += positivity_preservation_coeff * k2 * 2 * velocity_magnitude *
                                 velocity_magnitude_derivatives(c, k) *
                                 primal_variable_gradient_convective_terms[a];
                        value -= positivity_preservation_coeff_derivatives(c, k) *
                                 k2 * velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient;
                        value -= positivity_preservation_coeff *
                                 cross_wind_diffusion_coeff_derivatives(c, k) *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient;
                        value -= positivity_preservation_coeff * k2 *
                                 (gauss_shape_functions[c] * r_shape_derivatives(a, k) *
                                  velocity_dot_primal_variable_gradient);
                        value -= positivity_preservation_coeff * k2 *
                                 (gauss_shape_functions[c] * velocity_convective_terms[a] *
                                  primal_variable_gradient[k]);

                        // adding stream line dissipation derivatives
                        value += positivity_preservation_coeff_derivatives(c, k) *
                                 k1 * velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient;
                        value += positivity_preservation_coeff *
                                 stream_line_diffusion_coeff_derivatives(c, k) *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient;
                        value += positivity_preservation_coeff * k1 *
                                 gauss_shape_functions[c] * r_shape_derivatives(a, k) *
                                 velocity_dot_primal_variable_gradient;
                        value += positivity_preservation_coeff * k1 *
                                 velocity_convective_terms[a] *
                                 gauss_shape_functions[c] * primal_variable_gradient[k];

                        // putting transposed values
                        rResidualDerivatives(block_size + k, a) -=
                            value * gauss_weights[g];

                        // adding source term derivatives
                        value = 0.0;

                        value += gauss_shape_functions[a] * source_derivatives(c, k);
                        value += tau_derivatives(c, k) *
                                 (velocity_convective_terms[a] +
                                  s * gauss_shape_functions[a]) *
                                 source;
                        value += tau *
                                 (gauss_shape_functions[c] * r_shape_derivatives(a, k) +
                                  s_derivatives(c, k) * gauss_shape_functions[a]) *
                                 source;
                        value += tau *
                                 (velocity_convective_terms[a] +
                                  s * gauss_shape_functions[a]) *
                                 source_derivatives(c, k);

                        // putting transposed values
                        rResidualDerivatives(block_size + k, a) +=
                            value * gauss_weights[g];

                        // adding mass term derivatives
                        value = 0.0;

                        value += tau_derivatives(c, k) *
                                 (velocity_convective_terms[a] +
                                  s * gauss_shape_functions[a]) *
                                 primal_variable_relaxed_rate;
                        value += tau *
                                 (gauss_shape_functions[c] * r_shape_derivatives(a, k) +
                                  s_derivatives(c, k) * gauss_shape_functions[a]) *
                                 primal_variable_relaxed_rate;

                        rResidualDerivatives(block_size + k, a) -=
                            value * gauss_weights[g];
                    }
                }
            }
        }
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
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

}; // Class Element

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionAdjointData>
inline std::istream& operator>>(
    std::istream& rIStream,
    StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, TConvectionDiffusionReactionAdjointData>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionAdjointData>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, TConvectionDiffusionReactionAdjointData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}
///@}

} // namespace Kratos.
#endif // KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT defined
