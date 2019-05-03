//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT)
#define KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT

// System includes

// External includes

// Project includes
#include "includes/element.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_modelling_application_variables.h"

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
        return Kratos::make_shared<StabilizedConvectionDiffusionReactionAdjointElement>(
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
        if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

        rLeftHandSideMatrix.clear();
        CalculatePrimalDampingMatrixScalarDerivatives(rLeftHandSideMatrix, rCurrentProcessInfo);
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
        if (rLeftHandSideMatrix.size1() != 0)
            rLeftHandSideMatrix.resize(0, 0, false);
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

        KRATOS_ERROR_IF(this->Id() < 1)
            << "Element found with Id " << this->Id() << std::endl;

        const double domain_size = this->GetGeometry().DomainSize();
        KRATOS_ERROR_IF(domain_size <= 0.0) << "Element " << this->Id() << " has non-positive size "
                                            << domain_size << std::endl;

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
        buffer << "Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Element #" << Id();
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

    virtual void CalculateEffectiveKinematicViscosityDerivatives(
        Vector& rOutput,
        const Variable<double>& rDerivativeVariable,
        const TConvectionDiffusionReactionAdjointData& rCurrentData,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR
            << "Calling base CalculateEffectiveKinematicViscosityDerivatives "
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

    /**
     * pointer to the data related to this element
     */
    DataValueContainer mData;

    /**
     * pointer to the element's properties
     */
    Properties::Pointer mpProperties;

    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{

    void CalculateStabilizationTau(double& tau,
                                   double& element_length,
                                   const array_1d<double, 3>& rVelocity,
                                   const Matrix& rContravariantMetricTensor,
                                   const double reaction,
                                   const double effective_kinematic_viscosity,
                                   const double delta_time)
    {
        unsigned int dim = rContravariantMetricTensor.size2();
        Vector velocity(dim);
        for (unsigned int d = 0; d < dim; ++d)
            velocity[d] = rVelocity[d];
        Vector temp(dim);
        noalias(temp) = prod(rContravariantMetricTensor, velocity);
        const double stab_convection = inner_prod(velocity, temp);
        const double stab_diffusion = std::pow(3.0 * effective_kinematic_viscosity, 2) *
                                      norm_frobenius(rContravariantMetricTensor);
        const double stab_dynamics = std::pow(2.0 / delta_time, 2);
        const double stab_reaction = std::pow(reaction, 2);

        const double velocity_norm = norm_2(rVelocity);
        tau = 1.0 / std::sqrt(stab_dynamics + stab_convection + stab_diffusion + stab_reaction);
        element_length = 2.0 * velocity_norm / std::sqrt(stab_convection);
    }

    void CalculateStabilizationTauScalarDerivatives(Vector& rOutput,
                                                    const double tau,
                                                    const double effective_kinematic_viscosity,
                                                    const double reaction,
                                                    const Matrix& rContravariantMetricTensor,
                                                    const Vector& rEffectiveKinematicViscosityScalarDerivatives,
                                                    const Vector& rReactionScalarDerivatives)
    {
        if (rOutput.size() != rReactionScalarDerivatives.size())
            rOutput.resize(rReactionScalarDerivatives.size());

        noalias(rOutput) = (rEffectiveKinematicViscosityScalarDerivatives *
                                (18 * effective_kinematic_viscosity *
                                 norm_frobenius(rContravariantMetricTensor)) +
                            rReactionScalarDerivatives * (2.0 * reaction)) *
                           (-1.0 * std::pow(tau, 3) / 2.0);
    }

    void CalculateCrossWindDiffusionParameters(double& chi,
                                               double& k1,
                                               double& k2,
                                               const double velocity_magnitude,
                                               const double tau,
                                               const double effective_kinematic_viscosity,
                                               const double reaction,
                                               const double alpha,
                                               const double gamma,
                                               const double delta_time,
                                               const double element_length)
    {
        const double reaction_dynamics = reaction + (1 - alpha) / (gamma * delta_time);

        chi = 2.0 / (std::abs(reaction_dynamics) * element_length + 2.0 * velocity_magnitude);

        const double psi_one =
            this->CalculatePsiOne(velocity_magnitude, tau, reaction_dynamics);
        const double psi_two = this->CalculatePsiTwo(reaction_dynamics, tau, element_length);

        double value = 0.0;

        value = 0.5 * std::abs(psi_one - tau * velocity_magnitude * reaction_dynamics) * element_length;
        value -= (effective_kinematic_viscosity + tau * std::pow(velocity_magnitude, 2));
        value += psi_two;

        k1 = std::max(value, 0.0);

        value = 0.5 * std::abs(psi_one) * element_length;
        value -= effective_kinematic_viscosity;
        value += psi_two;

        k2 = std::max(value, 0.0);
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

    void CalculateAbsoluteScalarGradientScalarDerivative(Vector& rOutput,
                                                         const array_1d<double, 3> rScalarGradient,
                                                         const Matrix& rShapeFunctionDerivatives)
    {
        std::size_t number_of_nodes = rShapeFunctionDerivatives.size1();
        const double scalar_gradient_norm = norm_2(rScalarGradient);

        if (rOutput.size() != number_of_nodes)
            rOutput.resize(number_of_nodes);

        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const Vector& shape_function_gradient = row(rShapeFunctionDerivatives, i_node);
            rOutput[i_node] =
                CalculateScalarProduct(shape_function_gradient, rScalarGradient) / scalar_gradient_norm;
        }
    }

    void CalculateResidualScalarDerivative(Vector& rOutput,
                                           const double scalar_value,
                                           const double reaction,
                                           const array_1d<double, 3>& rVelocity,
                                           const Vector& rReactionScalarDerivatives,
                                           const Vector& rSourceScalarDerivatives,
                                           const Vector& rShapeFunctions,
                                           const Matrix& rShapeFunctionDerivatives)
    {
        std::size_t number_of_nodes = rSourceScalarDerivatives.size();
        if (rOutput.size() != number_of_nodes)
            rOutput.resize(number_of_nodes);

        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const Vector& shape_function_gradient = row(rShapeFunctionDerivatives, i_node);
            double value = 0.0;

            value += CalculateScalarProduct(shape_function_gradient, rVelocity);
            value += scalar_value * rReactionScalarDerivatives[i_node];
            value += reaction * rShapeFunctions[i_node];
            value -= rSourceScalarDerivatives[i_node];

            rOutput[i_node] = value;
        }
    }

    void CalculatePositivityPreservationCoefficientDerivatives(
        Vector& rOutput,
        const double chi,
        const double residual,
        const double scalar_gradient_norm,
        const double velocity_norm_square,
        const Vector& rChiScalarDerivatives,
        const Vector& rAbsoluteResidualScalarDerivatives,
        const Vector& rAbsoluteScalarGradientScalarDerivative)
    {
        const double abs_residual = std::abs(residual);

        if (rOutput.size() != rChiScalarDerivatives.size())
            rOutput.resize(rChiScalarDerivatives.size());

        noalias(rOutput) = rAbsoluteResidualScalarDerivatives *
                           (chi / (velocity_norm_square * scalar_gradient_norm));
        noalias(rOutput) +=
            rChiScalarDerivatives *
            (abs_residual / (velocity_norm_square * scalar_gradient_norm));
        noalias(rOutput) -=
            rAbsoluteScalarGradientScalarDerivative *
            (chi * abs_residual / (std::pow(scalar_gradient_norm, 2) * velocity_norm_square));
    }

    double CalculatePsiOne(const double velocity_norm, const double tau, const double reaction_tilde)
    {
        return velocity_norm + tau * velocity_norm * reaction_tilde;
    }

    double CalculatePsiTwo(const double reaction_tilde, const double tau, const double element_length)
    {
        return (reaction_tilde + tau * reaction_tilde * std::abs(reaction_tilde)) *
               std::pow(element_length, 2) / 6.0;
    }

    void CalculatePsiOneScalarDerivatives(Vector& rOutput,
                                          const double velocity_norm,
                                          const double reaction_tilde,
                                          const double tau,
                                          const Vector& rTauScalarDerivatives,
                                          const Vector& rAbsoluteReactionTildeScalarDerivatives)
    {
        if (rOutput.size() != rTauScalarDerivatives.size())
            rOutput.resize(rTauScalarDerivatives.size());

        const double absolute_reaction_tilde = std::abs(reaction_tilde);

        noalias(rOutput) = rTauScalarDerivatives * (velocity_norm * absolute_reaction_tilde);
        noalias(rOutput) += rAbsoluteReactionTildeScalarDerivatives * (tau * velocity_norm);
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

        if (rOutput.size() != rTauScalarDerivatives.size())
            rOutput.resize(rTauScalarDerivatives.size());

        noalias(rOutput) = rReactionTildeDerivatives;
        noalias(rOutput) +=
            rTauScalarDerivatives * (reaction_tilde * absolute_reaction_tilde);
        noalias(rOutput) += rReactionTildeDerivatives * (tau * absolute_reaction_tilde);
        noalias(rOutput) += rAbsoluteReactionTildeScalarDerivatives * (tau * reaction_tilde);
        noalias(rOutput) = rOutput * (std::pow(element_length, 2) / 6.0);
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
        if (rOutput.size() != rPsiOneScalarDerivatives.size())
            rOutput.resize(rPsiOneScalarDerivatives.size());

        noalias(rOutput) = rPsiOneScalarDerivatives;
        noalias(rOutput) -= rTauScalarDerivatives * (velocity_norm * reaction_tilde);
        noalias(rOutput) -= rReactionTildeScalarDerivatives * (tau * velocity_norm);

        const double coeff = psi_one - tau * velocity_norm * reaction_tilde;
        noalias(rOutput) = rOutput * (0.5 * element_length * (coeff) / std::abs(coeff));

        noalias(rOutput) += rPsiTwoScalarDerivatives;
        noalias(rOutput) -= rEffectiveViscosityScalarDerivatives;
        noalias(rOutput) -= rTauScalarDerivatives * std::pow(velocity_norm, 2);
    }

    void CalculateCrossWindDiffusionCoeffScalarDerivatives(
        Vector& rOutput,
        const double psi_one,
        const double element_length,
        const Vector& rPsiOneScalarDerivatives,
        const Vector& rPsiTwoScalarDerivatives,
        const Vector& rEffectiveKinematicViscosityScalarDerivatives)
    {
        if (rOutput.size() != rPsiOneScalarDerivatives.size())
            rOutput.resize(rPsiOneScalarDerivatives.size());

        noalias(rOutput) = rPsiOneScalarDerivatives *
                           (0.5 * psi_one * element_length / std::abs(psi_one));
        noalias(rOutput) -= rEffectiveKinematicViscosityScalarDerivatives;
        noalias(rOutput) += rPsiTwoScalarDerivatives;
    }

    void CalculateAbsoluteScalarValueScalarDerivatives(Vector& rOutput,
                                                       const double scalar_value,
                                                       const Vector& rScalarValueDerivatives)
    {
        if (rOutput.size() != rScalarValueDerivatives.size())
            rOutput.resize(rScalarValueDerivatives.size());

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

    void CalculatePrimalDampingMatrixScalarDerivatives(MatrixType& rLeftHandSideMatrix,
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

        const double delta_time = rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma = rCurrentProcessInfo[NEWMARK_GAMMA];

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
            Matrix contravariant_metric_tensor(r_parameter_derivatives_g.size1(),
                                               r_parameter_derivatives_g.size2());
            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            const array_1d<double, 3> velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            BoundedVector<double, TNumNodes> velocity_convective_terms;
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            TConvectionDiffusionReactionAdjointData current_data;
            this->CalculateConvectionDiffusionReactionAdjointData(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const Variable<double>& primal_variable = this->GetPrimalVariable();

            const double scalar_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);

            array_1d<double, 3> scalar_gradient;
            this->CalculateGradient(scalar_gradient, primal_variable, r_shape_derivatives);
            BoundedVector<double, TNumNodes> scalar_convective_terms;
            this->GetConvectionOperator(scalar_convective_terms,
                                        scalar_gradient, r_shape_derivatives);

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);
            Vector effective_kinematic_viscosity_derivatives;
            this->CalculateEffectiveKinematicViscosityDerivatives(
                effective_kinematic_viscosity_derivatives, primal_variable,
                current_data, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);
            Vector reaction_derivatives;
            this->CalculateReactionTermDerivatives(
                reaction_derivatives, primal_variable, current_data, rCurrentProcessInfo);

            const double source =
                this->CalculateSourceTerm(current_data, rCurrentProcessInfo);
            Vector source_derivatives;
            this->CalculateSourceTermDerivatives(source_derivatives, primal_variable,
                                                 current_data, rCurrentProcessInfo);

            double tau, element_length;
            this->CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, delta_time);

            Vector tau_derivatives;
            this->CalculateStabilizationTauScalarDerivatives(
                tau_derivatives, tau, effective_kinematic_viscosity, reaction,
                contravariant_metric_tensor,
                effective_kinematic_viscosity_derivatives, reaction_derivatives);

            const double s = std::abs(reaction);
            Vector s_derivatives;
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

                this->CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau, effective_kinematic_viscosity,
                    reaction, bossak_alpha, bossak_gamma, delta_time, element_length);

                positivity_preserving_coeff =
                    residual * chi / (velocity_magnitude_square * scalar_gradient_norm);
            }

            Vector chi_derivatives;
            this->CalculateChiScalarDerivatives(
                chi_derivatives, chi, element_length, bossak_alpha,
                bossak_gamma, delta_time, reaction, reaction_derivatives);

            Vector scalar_gradient_norm_derivative;
            this->CalculateAbsoluteScalarGradientScalarDerivative(
                scalar_gradient_norm_derivative, scalar_gradient, r_shape_derivatives);

            Vector residual_derivatives;
            this->CalculateResidualScalarDerivative(
                residual_derivatives, scalar_value, reaction, velocity, reaction_derivatives,
                source_derivatives, gauss_shape_functions, r_shape_derivatives);

            Vector absolute_residual_derivatives;
            this->CalculateAbsoluteScalarValueScalarDerivatives(
                absolute_residual_derivatives, residual, residual_derivatives);

            Vector positivity_preserving_coeff_derivatives;
            this->CalculatePositivityPreservationCoefficientDerivatives(
                positivity_preserving_coeff_derivatives, chi, residual,
                scalar_gradient_norm, velocity_magnitude_square, chi_derivatives,
                absolute_residual_derivatives, scalar_gradient_norm_derivative);

            const double reaction_tilde =
                reaction + (1 - bossak_alpha) / (bossak_gamma * delta_time);
            Vector absolute_reaction_tilde_derivatives;
            this->CalculateAbsoluteScalarValueScalarDerivatives(
                absolute_reaction_tilde_derivatives, reaction_tilde, reaction_derivatives);

            const double psi_one =
                this->CalculatePsiOne(velocity_magnitude, tau, reaction_tilde);
            Vector psi_one_derivatives;
            this->CalculatePsiOneScalarDerivatives(
                psi_one_derivatives, velocity_magnitude, reaction_tilde, tau,
                tau_derivatives, absolute_reaction_tilde_derivatives);

            const double psi_two =
                this->CalculatePsiTwo(reaction_tilde, tau, element_length);
            Vector psi_two_derivatives;
            this->CalculatePsiTwoScalarDerivatives(
                psi_two_derivatives, element_length, tau, reaction_tilde, tau_derivatives,
                reaction_derivatives, absolute_reaction_tilde_derivatives);

            Vector streamline_diffusion_coeff_derivatives;
            this->CalculateStreamLineDiffusionCoeffScalarDerivatives(
                streamline_diffusion_coeff_derivatives, element_length, tau,
                velocity_magnitude, reaction_tilde, psi_one, psi_two,
                psi_one_derivatives, psi_two_derivatives, tau_derivatives,
                reaction_derivatives, effective_kinematic_viscosity_derivatives);

            Vector crosswind_diffusion_coeff_derivatives;
            this->CalculateCrossWindDiffusionCoeffScalarDerivatives(
                crosswind_diffusion_coeff_derivatives, psi_one, element_length, psi_one_derivatives,
                psi_two_derivatives, effective_kinematic_viscosity_derivatives);

            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    double dNa_dNc = 0.0;
                    for (unsigned int i = 0; i < TDim; i++)
                        dNa_dNc += r_shape_derivatives(a, i) * r_shape_derivatives(c, i);

                    double value = 0.0;

                    // adding derivative of the convective term
                    value += gauss_shape_functions[a] * velocity_convective_terms[c];

                    // adding derivative of the diffusion term
                    value += scalar_convective_terms[a] *
                             effective_kinematic_viscosity_derivatives[c];
                    value += effective_kinematic_viscosity * dNa_dNc;

                    // adding reaction term derivatives
                    value += gauss_shape_functions[a] * reaction_derivatives[c] * scalar_value;
                    value += gauss_shape_functions[a] * reaction *
                             gauss_shape_functions[c];

                    // adding SUPG stabilization derivatives
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             velocity_convective_terms[c];
                    value += tau_derivatives[c] *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             velocity_dot_scalar_gradient;
                    value += tau * s_derivatives[c] * gauss_shape_functions[a] *
                             velocity_dot_scalar_gradient;

                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             reaction * gauss_shape_functions[c];
                    value += tau_derivatives[c] *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             reaction * scalar_value;
                    value += tau * (s_derivatives[c] * gauss_shape_functions[a]) *
                             reaction * scalar_value;
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             reaction_derivatives[c] * scalar_value;

                    // Adding cross wind dissipation derivatives
                    value += positivity_preserving_coeff * k2 * dNa_dNc * velocity_magnitude_square;
                    value += positivity_preserving_coeff_derivatives[c] * k2 *
                             scalar_convective_terms[a] * velocity_magnitude_square;
                    value += positivity_preserving_coeff *
                             crosswind_diffusion_coeff_derivatives[c] *
                             scalar_convective_terms[a] * velocity_magnitude_square;

                    value -= positivity_preserving_coeff * k2 *
                             velocity_convective_terms[a] *
                             velocity_convective_terms[c];
                    value -= positivity_preserving_coeff_derivatives[c] * k2 *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;
                    value -= positivity_preserving_coeff *
                             crosswind_diffusion_coeff_derivatives[c] *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;

                    // Adding stream line dissipation derivatives
                    value += positivity_preserving_coeff * k1 *
                             velocity_convective_terms[a] *
                             velocity_convective_terms[c];
                    value += positivity_preserving_coeff_derivatives[c] * k1 *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;
                    value += positivity_preserving_coeff *
                             streamline_diffusion_coeff_derivatives[c] *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;

                    rLeftHandSideMatrix(a, c) += gauss_weights[g] * value;
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GeometricalObject);
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags);
        rSerializer.save("Data", mData);
        rSerializer.save("Properties", mpProperties);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GeometricalObject);
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags);
        rSerializer.load("Data", mData);
        rSerializer.load("Properties", mpProperties);
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
