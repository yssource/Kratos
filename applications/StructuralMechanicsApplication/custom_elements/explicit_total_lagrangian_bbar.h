// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_EXPLICIT_TOTAL_LAGRANGIAN_BBAR_ELEMENT_H_INCLUDED )
#define  KRATOS_EXPLICIT_TOTAL_LAGRANGIAN_BBAR_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/integration_utilities.h"
#include "structural_mechanics_application_variables.h"
namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

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
 * @class ExplicitTotalLagrangianBbar
 * @ingroup StructuralMechanicsApplication
 * @brief This is an element explicitely designed for explicit time integration
 * @details The element uses a Bbar stabilization
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the element
 */
template< const SizeType TDim, const SizeType TNumNodes>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ExplicitTotalLagrangianBbar
    : public Element
{
public:

    ///@name Type Definitions
    ///@{

    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;

    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;

    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// This is the definition of the node.
    typedef Node<3> NodeType;

    /// The base element type
    typedef Element BaseType;

    static constexpr SizeType StrainSize = TDim == 3 ? 6 : 3;
//     static constexpr SizeType StrainSize = TDim == 3 ? 6 : 4;

    static constexpr SizeType ElementSize = TDim * TNumNodes;

    // Counted pointer of ExplicitTotalLagrangianBbar
    KRATOS_CLASS_POINTER_DEFINITION( ExplicitTotalLagrangianBbar );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    ExplicitTotalLagrangianBbar()
    {};

    // Constructor using an array of nodes
    ExplicitTotalLagrangianBbar( IndexType NewId, GeometryType::Pointer pGeometry ):Element(NewId,pGeometry)
    {};

    // Constructor using an array of nodes with properties
    ExplicitTotalLagrangianBbar( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):Element(NewId,pGeometry,pProperties)
    {};

    // Copy constructor
    ExplicitTotalLagrangianBbar(ExplicitTotalLagrangianBbar const& rOther)
        :BaseType(rOther)
        ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    {};

    // Destructor
    ~ExplicitTotalLagrangianBbar() override
    {};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Called to initialize the element.
     * @warning Must be called before any calculation is done
     */
    void Initialize() override;

    /**
      * @brief This resets the constitutive law
      */
    void ResetConstitutiveLaw() override;

    /**
     * @brief Called at the beginning of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the end of eahc solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

    /**
     * @brief Returns the used integration method
     * @return default integration method of the used Geometry
     */
    IntegrationMethod GetIntegrationMethod() const override
    {
        return mThisIntegrationMethod;
    }

    /**
     * @brief Sets on rValues the nodal displacements
     * @param rValues The values of displacements
     * @param Step The step to be computed
     */
    void GetValuesVector(
        Vector& rValues,
        int Step = 0
        ) override;

    /**
     * @brief Sets on rValues the nodal velocities
     * @param rValues The values of velocities
     * @param Step The step to be computed
     */
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) override;

    /**
     * @brief Sets on rValues the nodal accelerations
     * @param rValues The values of accelerations
     * @param Step The step to be computed
     */
    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
      * @param rRightHandSideVector the elemental right hand side vector
      * @param rCurrentProcessInfo the current process info instance
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable rDestinationVariable (double version)
     * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT IS ALLOWED TO WRITE ON ITS NODES.
     * The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rRHSVector input variable containing the RHS vector to be assembled
     * @param rRHSVariable variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable variable in the database to which the rRHSVector will be assembled
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(
        const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        Variable<double >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable (array_1d<double, 3>) version rDestinationVariable.
     * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT IS ALLOWED TO WRITE ON ITS NODES.
     * The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rRHSVector input variable containing the RHS vector to be assembled
     * @param rRHSVariable variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable variable in the database to which the rRHSVector will be assembled
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        Variable<array_1d<double, 3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix The elemental mass matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental damping matrix
      * @param rDampingMatrix The elemental damping matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a boolean Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<bool>& rVariable,
        std::vector<bool>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a integer Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<int>& rVariable,
        std::vector<int>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a double Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a 3 components array_1d on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a 6 components array_1d on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 6>>& rVariable,
        std::vector<array_1d<double, 6>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a Vector Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a Matrix Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief Set a bool Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<bool>& rVariable,
        std::vector<bool>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief Set a int Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<int>& rVariable,
        std::vector<int>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief Set a double Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief Set a Vector Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief Set a Matrix Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief Set a Constitutive Law Value on the Element
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
     void SetValueOnIntegrationPoints(
         const Variable<ConstitutiveLaw::Pointer>& rVariable,
         std::vector<ConstitutiveLaw::Pointer>& rValues,
         const ProcessInfo& rCurrentProcessInfo
         ) override;

     /**
      * @brief Set a array of 3 compoenents Value on the Element
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
     void SetValueOnIntegrationPoints(
         const Variable<array_1d<double, 3 > >& rVariable,
         std::vector<array_1d<double, 3 > > rValues,
         const ProcessInfo& rCurrentProcessInfo
         ) override;

     /**
      * @brief Set a array of 6 compoenents Value on the Element
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
     void SetValueOnIntegrationPoints(
         const Variable<array_1d<double, 6 > >& rVariable,
         std::vector<array_1d<double, 6 > > rValues,
         const ProcessInfo& rCurrentProcessInfo
         ) override;

    // GetValueOnIntegrationPoints are TEMPORARY until they are removed!!!
    // They will be removed from the derived elements; i.e. the implementation
    // should be in CalculateOnIntegrationPoints!
    // Adding these functions here is bcs GiD calls GetValueOnIntegrationPoints

    /**
     * @brief Get on rVariable a bool Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<bool>& rVariable,
        std::vector<bool>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a int Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<int>& rVariable,
        std::vector<int>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a double Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a 3  components array_1d Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a 6 components array_1d Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 6>>& rVariable,
        std::vector<array_1d<double, 6>>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a Vector Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a Matrix Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @todo To be renamed to CalculateOnIntegrationPoints!!!
     * @brief Get on rVariable Constitutive Law from the element
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo the current process info instance
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) override;

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
        std::stringstream buffer;
        buffer << "Base Solid Element #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Base Solid Element #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    /**
     * @brief Internal variables used in the kinematic calculations
     */
    struct KinematicVariables
    {
        Vector  N;
        BoundedMatrix<double, StrainSize, TDim * TNumNodes>  B;
        array_1d<double, TDim * TNumNodes> Bh;
        double  detF;
        Matrix  F;
        double  detJ0;
        BoundedMatrix<double, TDim, TDim> J0;
        BoundedMatrix<double, TDim, TDim> InvJ0;
        BoundedMatrix<double, TNumNodes, TDim>  DN_DX;

        /**
         * @brief The default constructor
         */
        KinematicVariables()
        {
            detF = 1.0;
            detJ0 = 1.0;
            N = ZeroVector(TNumNodes);
            B = ZeroMatrix(StrainSize, TDim * TNumNodes);
            Bh = ZeroVector(TDim * TNumNodes);
            F = IdentityMatrix(TDim);
            DN_DX = ZeroMatrix(TNumNodes, TDim);
            J0 = ZeroMatrix(TDim, TDim);
            InvJ0 = ZeroMatrix(TDim, TDim);
        }
    };

    /**
     * @brief Internal variables used in the kinematic calculations
     */
    struct ConstitutiveVariables
    {
        Vector StrainVector;
        Vector StressVector;
        Matrix D;

        /**
         * @brief The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         */
        ConstitutiveVariables()
        {
            StrainVector = ZeroVector(StrainSize);
            StressVector = ZeroVector(StrainSize);
            D = ZeroMatrix(StrainSize, StrainSize);
        }
    };

    IntegrationMethod mThisIntegrationMethod; /// Currently selected integration methods

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; /// The vector containing the constitutive laws

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Sets the used integration method
     * @param ThisIntegrationMethod Integration method used
     */
    void SetIntegrationMethod(const IntegrationMethod& ThisIntegrationMethod)
    {
         mThisIntegrationMethod = ThisIntegrationMethod;
    }

    /**
     * @brief Sets the used constitutive laws
     * @param ThisConstitutiveLawVector Constitutive laws used
     */
    void SetConstitutiveLawVector(const std::vector<ConstitutiveLaw::Pointer>& ThisConstitutiveLawVector)
    {
        mConstitutiveLawVector = ThisConstitutiveLawVector;
    }

    /**
     * @brief It initializes the material
     */
    void InitializeMaterial();

    /**
     * @brief This functions updates the kinematics variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param PointNumber The integration point considered
     * @param rIntegrationMethod The integration method considered
     */
    void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod
        );

    /**
     * @brief This functions updates the data structure passed to the CL
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     */
    void SetConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        );

    /**
     * @brief This functions updates the constitutive variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     * @param ThisStressMeasure The stress measure considered
     */
    void CalculateConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure = ConstitutiveLaw::StressMeasure_PK2
        );

    /**
     * @brief This method computes the deformation matrix
     * @param rB The deformation matrix
     * @param rF The deformation gradient
     * @param rDN_DX The shape functions derivatives
     */
    void CalculateB(
        BoundedMatrix<double, StrainSize, TDim * TNumNodes>& rB,
        const Matrix& rF,
        const BoundedMatrix<double, TNumNodes, TDim>& rDN_DX
        );

    /**
     * @brief This functions calculate the derivatives in the reference frame
     * @param rJ0 The jacobian in the reference configuration
     * @param rInvJ0 The inverse of the jacobian in the reference configuration
     * @param rDN_DX The gradient derivative of the shape function
     * @param PointNumber The id of the integration point considered
     * @param ThisIntegrationMethod The integration method considered
     * @return The determinant of the jacobian in the reference configuration
     */
    double CalculateDerivativesOnReferenceConfiguration(
        BoundedMatrix<double, TDim, TDim>& rJ0,
        BoundedMatrix<double, TDim, TDim>& rInvJ0,
        BoundedMatrix<double, TNumNodes, TDim>& rDN_DX,
        const IndexType PointNumber,
        IntegrationMethod ThisIntegrationMethod
        ) const;

    /**
     * @brief This functions calculate the derivatives in the current frame
     * @param rJ The jacobian in the current configuration
     * @param rInvJ The inverse of the jacobian in the current configuration
     * @param rDN_DX The gradient derivative of the shape function
     * @param PointNumber The id of the integration point considered
     * @param ThisIntegrationMethod The integration method considered
     * @return The determinant of the jacobian in the current configuration
     */
    double CalculateDerivativesOnCurrentConfiguration(
        BoundedMatrix<double, TDim, TDim>& rJ,
        BoundedMatrix<double, TDim, TDim>& rInvJ,
        BoundedMatrix<double, TNumNodes, TDim>& rDN_DX,
        const IndexType PointNumber,
        IntegrationMethod ThisIntegrationMethod
        ) const;

    /**
     * @brief This function computes the body force
     * @param IntegrationPoints The array containing the integration points
     * @param PointNumber The id of the integration point considered
     * @return The vector of body forces
     */
    array_1d<double, 3> GetBodyForce(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const IndexType PointNumber
        ) const;

    /**
     * @brief Calculation of the RHS
     * @param rRightHandSideVector The local component of the RHS due to external forces
     * @param rThisKinematicVariables The kinematic variables like shape functions
     * @param rCurrentProcessInfo rCurrentProcessInfo the current process info instance
     * @param rBodyForce The component of external forces due to self weight
     * @param rStressVector The vector containing the stress components
     * @param IntegrationWeight The integration weight of the corresponding Gauss point
     */
    void CalculateAndAddResidualVector(
        VectorType& rRightHandSideVector,
        const KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rCurrentProcessInfo,
        const array_1d<double, 3>& rBodyForce,
        const Vector& rStressVector,
        const double IntegrationWeight
        ) const;

    /**
     * @brief This function add the external force contribution
     * @param rN Shape function
     * @param rCurrentProcessInfo rCurrentProcessInfo the current process info instance
     * @param rBodyForce The component of external forces due to self weight
     * @param rRightHandSideVector The local component of the RHS due to external forces
     * @param IntegrationWeight The integration weight of the corresponding Gauss point
     */
    void CalculateAndAddExtForceContribution(
        const Vector& rN,
        const ProcessInfo& rCurrentProcessInfo,
        const array_1d<double, 3>& rBodyForce,
        VectorType& rRightHandSideVector,
        const double IntegrationWeight
        ) const;

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
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

    /**
     * @brief This method computes directly the lumped mass vector
     * @param rMassMatrix The lumped mass vector
     */
    void CalculateLumpedMassVector(array_1d<double, ElementSize>& rMassVector) const;

    /**
     * @brief This method computes directly the lumped mass matrix
     * @param rMassMatrix The lumped mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateLumpedDampingVector(
        array_1d<double, ElementSize>& rDampingVector,
        const ProcessInfo& rCurrentProcessInfo
        );

    /**
     * @brief This method gets a value directly in the CL
     * @details Avoids code repetition
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @tparam TType The type considered
     */
    template<class TType>
    void GetValueOnConstitutiveLaw(
        const Variable<TType>& rVariable,
        std::vector<TType>& rOutput
        )
    {
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        for ( IndexType point_number = 0; point_number <integration_points.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->GetValue( rVariable,rOutput[point_number]);
        }
    }

    /**
     * @brief This method computes directly in the CL
     * @details Avoids code repetition
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     * @tparam TType The type considered
     */
    template<class TType>
    void CalculateOnConstitutiveLaw(
        const Variable<TType>& rVariable,
        std::vector<TType>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

        KinematicVariables this_kinematic_variables = KinematicVariables();
        ConstitutiveVariables this_constitutive_variables = ConstitutiveVariables();

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
            // Compute element kinematics B, F, DN_DX ...
            this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

            // Compute material reponse
            this->SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

            rOutput[point_number] = mConstitutiveLawVector[point_number]->CalculateValue( Values, rVariable, rOutput[point_number] );
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override;

    void load( Serializer& rSerializer ) override;

}; // class ExplicitTotalLagrangianBbar.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_EXPLICIT_TOTAL_LAGRANGIAN_BBAR_ELEMENT_H_INCLUDED  defined
