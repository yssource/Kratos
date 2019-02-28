//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Riccardo Rossi
//
#include "embedded_incompressible_potential_flow_element.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
Element::Pointer EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_shared<EmbeddedIncompressiblePotentialFlowElement>(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_shared<EmbeddedIncompressiblePotentialFlowElement>(
        NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_shared<EmbeddedIncompressiblePotentialFlowElement>(
        NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties());
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    //NOT EMBEDDED
    BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // TODO: improve speed
    Matrix tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    BaseType::EquationIdVector(rResult, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                                   ProcessInfo& rCurrentProcessInfo)
{
    BaseType::GetDofList(rElementalDofList, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int Dim, int NumNodes>
int EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Generic geometry check
    int out = BaseType::Check(rCurrentProcessInfo);
    if (out != 0)
    {
        return out;
    }

    KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0.0)
        << this->Id() << "Area cannot be less than or equal to 0" << std::endl;

    for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL,this->GetGeometry()[i]);
    }

    return out;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <int Dim, int NumNodes>
std::string EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "EmbeddedIncompressiblePotentialFlowElement #" << this->Id();
    return buffer.str();
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EmbeddedIncompressiblePotentialFlowElement #" << this->Id();
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::PrintData(std::ostream& rOStream) const
{
    this->pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class EmbeddedIncompressiblePotentialFlowElement<2, 3>;

} // namespace Kratos
