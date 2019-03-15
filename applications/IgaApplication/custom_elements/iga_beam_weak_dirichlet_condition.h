/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
//                  Lukas Rauch
*/

#if !defined(KRATOS_IGA_BEAM_WEAK_DIRICHLET_CONDITION_H_INCLUDED)
#define KRATOS_IGA_BEAM_WEAK_DIRICHLET_CONDITION_H_INCLUDED


// System includes
#include "includes/define.h"
#include "includes/element.h"

// External includes
#include <Eigen/Core>
#include "C:\_Masterarbeit\bibliotheken\HyperJet\HyperJet.h"

// Project includes
#include "iga_base_element.h"


namespace Kratos
{

class IgaBeamWeakDirichletCondition
    : public IgaBaseElement<4>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( IgaBeamWeakDirichletCondition );

    using IgaBaseElementType::IgaBaseElementType;

    ~IgaBeamWeakDirichletCondition() override
    {
    };

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo) override;

    void Initialize() override;

    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool ComputeLeftHandSide,
        const bool ComputeRightHandSide) override;

    void PrintInfo(std::ostream& rOStream) const override;
};

} // namespace Kratos

#endif // !defined(KRATOS_IGA_BEAM_WEAK_DIRICHLET_CONDITION_H_INCLUDED)
