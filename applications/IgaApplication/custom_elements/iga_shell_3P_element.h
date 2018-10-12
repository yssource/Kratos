/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

#if !defined(KRATOS_IGA_SHELL_3P_ELEMENT_H_INCLUDED)
#define KRATOS_IGA_SHELL_3P_ELEMENT_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/element.h"

// External includes

// Project includes
#include "iga_base_element.h"


namespace Kratos
{

class IgaShell3PElement
    : public IgaBaseElement<3>
{
private:
    struct Configuration
    {
        Vector3 a1;
        Vector3 a2;
        Vector3 a3;
        
        Vector3 a11;
        Vector3 a12;
        Vector3 a22;
    };

public:
    KRATOS_CLASS_POINTER_DEFINITION( IgaShell3PElement );

    using IgaBaseElementType::IgaBaseElementType;

    ~IgaShell3PElement() override
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

    void PrintInfo(
        std::ostream& rOStream) const override;
};

} // namespace Kratos

#endif // !defined(KRATOS_IGA_SHELL_3P_ELEMENT_H_INCLUDED)
