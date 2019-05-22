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
    typedef EvmKEpsilonVMSAdjointElement<TDim> AdjointFluidElement;
    // defining the k element type
    typedef EvmKAdjointElement<TDim, TNumNodes> AdjointKElement;
    // defining the epsilon element type
    typedef EvmEpsilonAdjointElement<TDim, TNumNodes> AdjointEpsilonElement;

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
        // AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
        // AdjointKElement k_element(this->Id(), this->pGetGeometry());
        // AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

        // VectorType fluid_values;
        // fluid_element.GetValuesVector(fluid_values, Step);

        // VectorType k_values;
        // k_element.GetValuesVector(k_values, Step);

        // VectorType epsilon_values;
        // epsilon_element.GetValuesVector(epsilon_values, Step);

        // const int total_values_size =
        //     fluid_element.size() + k_element.size() + epsilon_element.size();

        // if (rValues.size() != total_values_size)
        //     rValue.resize(total_values_size, false);

        // IndexType current_index{0};
        // for (IndexType i = 0; i < fluid_element.size(); ++i)
        //     rValue[current_index++] = fluid_element[i];
        // for (IndexType i = 0; i < k_element.size(); ++i)
        //     rValue[current_index++] = k_element[i];
        // for (IndexType i = 0; i < epsilon_element.size(); ++i)
        //     rValue[current_index++] = epsilon_element[i];
    }

    /// Returns the adjoint velocity values stored in this element's nodes.
    void GetFirstDerivativesVector(VectorType& rValues, int Step = 0) override
    {
        // AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
        // AdjointKElement k_element(this->Id(), this->pGetGeometry());
        // AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

        // VectorType fluid_values;
        // fluid_element.GetFirstDerivativesVector(fluid_values, Step);

        // VectorType k_values;
        // k_element.GetFirstDerivativesVector(k_values, Step);

        // VectorType epsilon_values;
        // epsilon_element.GetFirstDerivativesVector(epsilon_values, Step);

        // const int total_values_size =
        //     fluid_element.size() + k_element.size() + epsilon_element.size();

        // if (rValues.size() != total_values_size)
        //     rValue.resize(total_values_size, false);

        // IndexType current_index{0};
        // for (IndexType i = 0; i < fluid_element.size(); ++i)
        //     rValue[current_index++] = fluid_element[i];
        // for (IndexType i = 0; i < k_element.size(); ++i)
        //     rValue[current_index++] = k_element[i];
        // for (IndexType i = 0; i < epsilon_element.size(); ++i)
        //     rValue[current_index++] = epsilon_element[i];
    }

    void GetSecondDerivativesVector(VectorType& rValues, int Step) override
    {
        // AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
        // AdjointKElement k_element(this->Id(), this->pGetGeometry());
        // AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

        // VectorType fluid_values;
        // fluid_element.GetSecondDerivativesVector(fluid_values, Step);

        // VectorType k_values;
        // k_element.GetSecondDerivativesVector(k_values, Step);

        // VectorType epsilon_values;
        // epsilon_element.GetSecondDerivativesVector(epsilon_values, Step);

        // const int total_values_size =
        //     fluid_element.size() + k_element.size() + epsilon_element.size();

        // if (rValues.size() != total_values_size)
        //     rValue.resize(total_values_size, false);

        // IndexType current_index{0};
        // for (IndexType i = 0; i < fluid_element.size(); ++i)
        //     rValue[current_index++] = fluid_element[i];
        // for (IndexType i = 0; i < k_element.size(); ++i)
        //     rValue[current_index++] = k_element[i];
        // for (IndexType i = 0; i < epsilon_element.size(); ++i)
        //     rValue[current_index++] = epsilon_element[i];
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_THROW_ERROR(std::runtime_error,
                           "this function is not implemented.", "")

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

        KRATOS_THROW_ERROR(std::runtime_error,
                           "this function is not implemented.", "")

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

        IndexType i_offset{0}, j_offset{0};
        Matrix local_matrix;
        RansCalculationUtilities rans_calculation_utilities;

        fluid_element.CalculateFirstDerivativesLHS(local_matrix, rCurrentProcessInfo);
        rans_calculation_utilities.PlaceInGlobalMatrix(
            rLeftHandSideMatrix, local_matrix, i_offset, j_offset);
        i_offset += local_matrix.size1();

        fluid_element.Calculate(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                                local_matrix, rCurrentProcessInfo);
        rans_calculation_utilities.PlaceInGlobalMatrix(
            rLeftHandSideMatrix, local_matrix, i_offset, j_offset);
        i_offset += local_matrix.size1();

        fluid_element.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                                local_matrix, rCurrentProcessInfo);
        rans_calculation_utilities.PlaceInGlobalMatrix(
            rLeftHandSideMatrix, local_matrix, i_offset, j_offset);

        i_offset = 0;
        j_offset += local_matrix.size2();

        k_element.Calculate(RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE,
                            local_matrix, rCurrentProcessInfo);
        rans_calculation_utilities.PlaceInGlobalMatrix(
            rLeftHandSideMatrix, local_matrix, i_offset, j_offset);
        i_offset += local_matrix.size1();

        k_element.CalculateFirstDerivativesLHS(local_matrix, rCurrentProcessInfo);
        rans_calculation_utilities.PlaceInGlobalMatrix(
            rLeftHandSideMatrix, local_matrix, i_offset, j_offset);
        i_offset += local_matrix.size1();

        k_element.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                            local_matrix, rCurrentProcessInfo);
        rans_calculation_utilities.PlaceInGlobalMatrix(
            rLeftHandSideMatrix, local_matrix, i_offset, j_offset);

        i_offset = 0;
        j_offset += local_matrix.size2();

        epsilon_element.Calculate(RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE,
                                  local_matrix, rCurrentProcessInfo);
        rans_calculation_utilities.PlaceInGlobalMatrix(
            rLeftHandSideMatrix, local_matrix, i_offset, j_offset);
        i_offset += local_matrix.size1();

        epsilon_element.Calculate(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                                  local_matrix, rCurrentProcessInfo);
        rans_calculation_utilities.PlaceInGlobalMatrix(
            rLeftHandSideMatrix, local_matrix, i_offset, j_offset);
        i_offset += local_matrix.size1();

        epsilon_element.CalculateFirstDerivativesLHS(local_matrix, rCurrentProcessInfo);
        rans_calculation_utilities.PlaceInGlobalMatrix(
            rLeftHandSideMatrix, local_matrix, i_offset, j_offset);

        KRATOS_CATCH("");
    }

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const double inv_one_minus_bossak_alpha = 1.0 / (1.0 - rCurrentProcessInfo[BOSSAK_ALPHA]);

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

        IndexType i_offset{0}, j_offset{0};
        Matrix local_matrix;
        RansCalculationUtilities rans_calculation_utilities;

        fluid_element.CalculateSecondDerivativesLHS(local_matrix, rCurrentProcessInfo);
        rans_calculation_utilities.PlaceInGlobalMatrix(
            rLeftHandSideMatrix, local_matrix, i_offset, j_offset);
        i_offset += local_matrix.size1();
        j_offset += local_matrix.size2();

        k_element.CalculateSecondDerivativesLHS(local_matrix, rCurrentProcessInfo);
        noalias(local_matrix) = local_matrix * inv_one_minus_bossak_alpha;
        rans_calculation_utilities.PlaceInGlobalMatrix(
            rLeftHandSideMatrix, local_matrix, i_offset, j_offset);
        i_offset += local_matrix.size1();
        j_offset += local_matrix.size2();

        epsilon_element.CalculateSecondDerivativesLHS(local_matrix, rCurrentProcessInfo);
        noalias(local_matrix) = local_matrix * inv_one_minus_bossak_alpha;
        rans_calculation_utilities.PlaceInGlobalMatrix(
            rLeftHandSideMatrix, local_matrix, i_offset, j_offset);

        KRATOS_CATCH("");
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