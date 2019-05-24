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

#if !defined(KRATOS_EVM_MONOLITHIC_K_EPSILON_VMS_ADJOINT_ELEMENT_H_INCLUDED)
#define KRATOS_EVM_MONOLITHIC_K_EPSILON_VMS_ADJOINT_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/rans_evm_vms_adjoint_element.h"
#include "includes/checks.h"
#include "includes/element.h"
#include "includes/properties.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_modelling_application_variables.h"

#include "custom_elements/evm_k_epsilon/evm_epsilon_adjoint_element.h"
#include "custom_elements/evm_k_epsilon/evm_k_adjoint_element.h"
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_vms_adjoint_element.h"

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

template <unsigned int TDim, unsigned int TNumNodes = TDim + 1>
class EvmMonolithicKEpsilonVMSAdjointElement : public Element
{
    // class ThisExtensions : public AdjointExtensions
    // {
    //     Element* mpElement;

    // public:
    //     explicit ThisExtensions(Element* pElement);

    //     void GetFirstDerivativesVector(std::size_t NodeId,
    //                                    std::vector<IndirectScalar<double>>&
    //                                    rVector, std::size_t Step) override;

    //     void GetSecondDerivativesVector(std::size_t NodeId,
    //                                     std::vector<IndirectScalar<double>>&
    //                                     rVector, std::size_t Step) override;

    //     void GetAuxiliaryVector(std::size_t NodeId,
    //                             std::vector<IndirectScalar<double>>& rVector,
    //                             std::size_t Step) override;

    //     void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables) const override;

    //     void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables) const override;

    //     void GetAuxiliaryVariables(std::vector<VariableData const*>& rVariables) const override;
    // };

public:
    ///@name Type Definitions
    ///@{

    // defining the base type
    typedef Element BaseType;
    // defining the base adjoint base fluid element type
    typedef EvmKEpsilonVMSAdjointElement<TDim, TDim + 1, TDim + 3> AdjointFluidElement;
    // defining the k element type
    typedef EvmKAdjointElement<TDim, TNumNodes, TDim + 3> AdjointKElement;
    // defining the epsilon element type
    typedef EvmEpsilonAdjointElement<TDim, TNumNodes, TDim + 3> AdjointEpsilonElement;

    constexpr static unsigned int TFluidBlockSize = (TDim + 1);

    constexpr static unsigned int TFluidLocalSize = TFluidBlockSize * TNumNodes;

    constexpr static unsigned int TKBlockSize = 1;

    constexpr static unsigned int TKLocalSize = TKBlockSize * TNumNodes;

    constexpr static unsigned int TEpsilonBlockSize = 1;

    constexpr static unsigned int TEpsilonLocalSize = TEpsilonBlockSize * TNumNodes;

    constexpr static unsigned int TElementBlockSize =
        (TFluidBlockSize + TKBlockSize + TEpsilonBlockSize);

    constexpr static unsigned int TElementLocalSize = TElementBlockSize * TNumNodes;

    constexpr static unsigned int TCoordLocalSize = TDim * TNumNodes;

    // variable definitions
    typedef std::size_t IndexType;

    typedef Element::NodeType NodeType;

    typedef Element::NodesArrayType NodesArrayType;

    typedef Element::GeometryType GeometryType;

    typedef Element::PropertiesType PropertiesType;

    typedef Element::VectorType VectorType;

    typedef Element::MatrixType MatrixType;

    ///@name Pointer Definitions
    /// Pointer definition of EvmMonolithicKEpsilonVMSAdjointElement
    KRATOS_CLASS_POINTER_DEFINITION(EvmMonolithicKEpsilonVMSAdjointElement);

    ///@}

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    EvmMonolithicKEpsilonVMSAdjointElement(IndexType NewId = 0)
        : BaseType(NewId)
    {
    }

    /**
     * Constructor using Geometry
     */
    EvmMonolithicKEpsilonVMSAdjointElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    EvmMonolithicKEpsilonVMSAdjointElement(IndexType NewId,
                                           GeometryType::Pointer pGeometry,
                                           PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Destructor
     */
    ~EvmMonolithicKEpsilonVMSAdjointElement() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        // TODO: Implement adjoint extensions
    }

    /**
     * ELEMENTS inherited from this class have to implement next
     * Create and Clone methods: MANDATORY
     */

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<EvmMonolithicKEpsilonVMSAdjointElement>(
            NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<EvmMonolithicKEpsilonVMSAdjointElement>(
            NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<EvmMonolithicKEpsilonVMSAdjointElement>(
            NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
        KRATOS_CATCH("");
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        BaseType::Check(rCurrentProcessInfo);

        AdjointFluidElement adjoint_fluid_element(this->Id(), this->pGetGeometry());
        AdjointKElement adjoint_k_element(this->Id(), this->pGetGeometry());
        AdjointEpsilonElement adjoint_epsilon_element(this->Id(), this->pGetGeometry());

        adjoint_fluid_element.Check(rCurrentProcessInfo);
        adjoint_k_element.Check(rCurrentProcessInfo);
        adjoint_epsilon_element.Check(rCurrentProcessInfo);

        return 0;

        KRATOS_CATCH("");
    }

    /// Returns the adjoint values stored in this element's nodes.
    void GetValuesVector(VectorType& rValues, int Step = 0) override
    {
        if (rValues.size() != TElementLocalSize)
            rValues.resize(TElementLocalSize, false);

        // AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
        // AdjointKElement k_element(this->Id(), this->pGetGeometry());
        // AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

        // fluid_element.SetData(this->Data());
        // k_element.SetData(this->Data());
        // epsilon_element.SetData(this->Data());

        // IndexType i_offset{0};
        // RansCalculationUtilities rans_calculation_utilities;
        // Vector local_vector;

        // fluid_element.GetValuesVector(local_vector, Step);
        // rans_calculation_utilities.PlaceInGlobalVector(rValues, local_vector,
        // i_offset); i_offset += local_vector.size();

        // k_element.GetValuesVector(local_vector, Step);
        // rans_calculation_utilities.PlaceInGlobalVector(rValues, local_vector,
        // i_offset); i_offset += local_vector.size();

        // epsilon_element.GetValuesVector(local_vector, Step);
        // rans_calculation_utilities.PlaceInGlobalVector(rValues, local_vector, i_offset);
    }

    /// Returns the adjoint velocity values stored in this element's nodes.
    void GetFirstDerivativesVector(VectorType& rValues, int Step = 0) override
    {
        if (rValues.size() != TElementLocalSize)
            rValues.resize(TElementLocalSize, false);

        // AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
        // AdjointKElement k_element(this->Id(), this->pGetGeometry());
        // AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

        // fluid_element.SetData(this->Data());
        // k_element.SetData(this->Data());
        // epsilon_element.SetData(this->Data());

        // IndexType i_offset{0};
        // RansCalculationUtilities rans_calculation_utilities;
        // Vector local_vector;

        // fluid_element.GetFirstDerivativesVector(local_vector, Step);
        // rans_calculation_utilities.PlaceInGlobalVector(rValues, local_vector,
        // i_offset); i_offset += local_vector.size();

        // k_element.GetFirstDerivativesVector(local_vector, Step);
        // rans_calculation_utilities.PlaceInGlobalVector(rValues, local_vector,
        // i_offset); i_offset += local_vector.size();

        // epsilon_element.GetFirstDerivativesVector(local_vector, Step);
        // rans_calculation_utilities.PlaceInGlobalVector(rValues, local_vector, i_offset);
    }

    void GetSecondDerivativesVector(VectorType& rValues, int Step) override
    {
        if (rValues.size() != TElementLocalSize)
            rValues.resize(TElementLocalSize, false);

        // AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
        // AdjointKElement k_element(this->Id(), this->pGetGeometry());
        // AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

        // fluid_element.SetData(this->Data());
        // k_element.SetData(this->Data());
        // epsilon_element.SetData(this->Data());

        // IndexType i_offset{0};
        // RansCalculationUtilities rans_calculation_utilities;
        // Vector local_vector;

        // fluid_element.GetSecondDerivativesVector(local_vector, Step);
        // rans_calculation_utilities.PlaceInGlobalVector(rValues, local_vector,
        // i_offset); i_offset += local_vector.size();

        // k_element.GetSecondDerivativesVector(local_vector, Step);
        // rans_calculation_utilities.PlaceInGlobalVector(rValues, local_vector,
        // i_offset); i_offset += local_vector.size();

        // epsilon_element.GetSecondDerivativesVector(local_vector, Step);
        // rans_calculation_utilities.PlaceInGlobalVector(rValues, local_vector, i_offset);
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "EvmMonolithicKEpsilonVMSAdjointElement::"
                        "CalculateLocalSystem method is not implemented.";

        KRATOS_CATCH("");
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
            rLeftHandSideMatrix.size2() != TElementLocalSize)
            rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);

        rLeftHandSideMatrix.clear();
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "EvmMonolithicKEpsilonVMSAdjointElement::"
                        "CalculateRightHandSide method is not implemented.";

        KRATOS_CATCH("");
    }

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
        AdjointKElement k_element(this->Id(), this->pGetGeometry());
        AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

        fluid_element.SetData(this->Data());
        k_element.SetData(this->Data());
        epsilon_element.SetData(this->Data());

        if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
            rLeftHandSideMatrix.size2() != TElementLocalSize)
            rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);

        rLeftHandSideMatrix.clear();

        fluid_element.CalculateFirstDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
        fluid_element.Calculate(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                                rLeftHandSideMatrix, rCurrentProcessInfo);
        fluid_element.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                                rLeftHandSideMatrix, rCurrentProcessInfo);

        k_element.Calculate(RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE,
                            rLeftHandSideMatrix, rCurrentProcessInfo);
        k_element.CalculateFirstDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
        k_element.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                            rLeftHandSideMatrix, rCurrentProcessInfo);

        epsilon_element.Calculate(RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE,
                                  rLeftHandSideMatrix, rCurrentProcessInfo);
        epsilon_element.Calculate(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                                  rLeftHandSideMatrix, rCurrentProcessInfo);
        epsilon_element.CalculateFirstDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
        AdjointKElement k_element(this->Id(), this->pGetGeometry());
        AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

        fluid_element.SetData(this->Data());
        k_element.SetData(this->Data());
        epsilon_element.SetData(this->Data());

        if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
            rLeftHandSideMatrix.size2() != TElementLocalSize)
            rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);

        rLeftHandSideMatrix.clear();

        fluid_element.CalculateSecondDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
        k_element.CalculateSecondDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
        epsilon_element.CalculateSecondDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "EvmMonolithicKEpsilonVMSAdjointElement::"
                        "CalculateMassMatrix method is not implemented.";

        KRATOS_CATCH("")
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "EvmMonolithicKEpsilonVMSAdjointElement::"
                        "CalculateDampingMatrix method is not implemented.";

        KRATOS_CATCH("")
    }

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rSensitivityVariable == SHAPE_SENSITIVITY)
        {
            if (rOutput.size1() != TCoordLocalSize || rOutput.size2() != TElementLocalSize)
                rOutput.resize(TCoordLocalSize, TElementLocalSize, false);

            rOutput.clear();

            AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
            AdjointKElement k_element(this->Id(), this->pGetGeometry());
            AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

            fluid_element.SetData(this->Data());
            k_element.SetData(this->Data());
            epsilon_element.SetData(this->Data());

            fluid_element.CalculateSensitivityMatrix(
                rSensitivityVariable, rOutput, rCurrentProcessInfo);
            k_element.CalculateSensitivityMatrix(
                rSensitivityVariable, rOutput, rCurrentProcessInfo);
            epsilon_element.CalculateSensitivityMatrix(
                rSensitivityVariable, rOutput, rCurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                         << " not supported." << std::endl;
        }

        KRATOS_CATCH("")
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
        std::stringstream buffer;
        buffer << "EvmMonolithicKEpsilonVMSAdjointElement #" << Element::Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EvmMonolithicKEpsilonVMSAdjointElement #" << Element::Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        Element::pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected Operations
    ///@{
    ///@}
private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

    ///@}
};

///@}

} // namespace Kratos

#endif