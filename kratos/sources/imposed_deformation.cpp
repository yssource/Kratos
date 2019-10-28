//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes


// External includes

// Project includes
#include "includes/imposed_deformation.h"

namespace Kratos
{
/// Flags related to the Parameters of the imposed deformation
KRATOS_CREATE_LOCAL_FLAG( ImposedDeformation, IS_INITIALIZED,  0 );

/***********************************************************************************/
/***********************************************************************************/

ImposedDeformation::Pointer ImposedDeformation::Clone() const
{
    KRATOS_ERROR << "Called the virtual function for Clone"<< std::endl;
    return nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

ImposedDeformation::Pointer ImposedDeformation::Create(Kratos::Parameters NewParameters) const
{
    const std::string& r_name = NewParameters["name"].GetString();
    return KratosComponents<ImposedDeformation>::Get(r_name).Clone();
}

/***********************************************************************************/
/***********************************************************************************/

bool ImposedDeformation::Has(const Variable<bool>& rThisVariable)
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool ImposedDeformation::Has(const Variable<int>& rThisVariable)
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool ImposedDeformation::Has(const Variable<double>& rThisVariable)
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool ImposedDeformation::Has(const Variable<Vector>& rThisVariable)
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool ImposedDeformation::Has(const Variable<Matrix>& rThisVariable)
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool ImposedDeformation::Has(const Variable<array_1d<double, 3> >& rThisVariable)
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool ImposedDeformation::Has(const Variable<array_1d<double, 6> >& rThisVariable)
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool& ImposedDeformation::GetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

int& ImposedDeformation::GetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

double& ImposedDeformation::GetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& ImposedDeformation::GetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& ImposedDeformation::GetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3>& ImposedDeformation::GetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<array_1d<double, 3>>& rThisVariable,
    array_1d<double, 3>& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 6 >& ImposedDeformation::GetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<array_1d<double, 6>>& rThisVariable,
    array_1d<double, 6>& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::SetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<bool>& rVariable,
    const bool& Value,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "Called the virtual function for SetValue"<< std::endl;;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::SetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<int>& rVariable,
    const int& Value,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "Called the virtual function for SetValue"<< std::endl;;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::SetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<double>& rVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "Called the virtual function for SetValue"<< std::endl;;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::SetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<Vector>& rVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "Called the virtual function for SetValue"<< std::endl;;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::SetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<Matrix>& rVariable,
    const Matrix& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "Called the virtual function for SetValue"<< std::endl;;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::SetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<array_1d<double, 3>>& rVariable,
    const array_1d<double, 3>& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "Called the virtual function for SetValue"<< std::endl;;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::SetValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Variable<array_1d<double, 6>>& rVariable,
    const array_1d<double, 6>& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "Called the virtual function for SetValue"<< std::endl;;
}

/***********************************************************************************/
/***********************************************************************************/

bool& ImposedDeformation::CalculateValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

int& ImposedDeformation::CalculateValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

double& ImposedDeformation::CalculateValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& ImposedDeformation::CalculateValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& ImposedDeformation::CalculateValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3>& ImposedDeformation::CalculateValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<array_1d<double, 3 > >& rVariable,
    array_1d<double, 3 > & rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 6>& ImposedDeformation::CalculateValue(
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<array_1d<double, 6 > >& rVariable,
    array_1d<double, 6 >& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::Initialize(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    this->Set(IS_INITIALIZED);
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::CalculateResponse (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues,
    const ConstitutiveLaw::StressMeasure& rStressMeasure
    )
{
    switch(rStressMeasure)
    {
        case ConstitutiveLaw::StressMeasure_PK1:
        CalculateResponsePK1(pConstitutiveLaw, rParameterValues);
        break;

    case ConstitutiveLaw::StressMeasure_PK2:
        CalculateResponsePK2(pConstitutiveLaw, rParameterValues);
        break;

    case ConstitutiveLaw::StressMeasure_Kirchhoff:
        CalculateResponseKirchhoff(pConstitutiveLaw, rParameterValues);
        break;

    case ConstitutiveLaw::StressMeasure_Cauchy:
        CalculateResponseCauchy(pConstitutiveLaw, rParameterValues);
        break;

    default:
        KRATOS_ERROR << " Stress Measure not Defined "<< std::endl;;
        break;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::CalculateResponsePK1 (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues
    )
{
    KRATOS_ERROR << "Calling virtual function for CalculateResponsePK1" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::CalculateResponsePK2 (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues
    )
{
    KRATOS_ERROR << "Calling virtual function for CalculateResponsePK2" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::CalculateResponseKirchhoff (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues
    )
{
    KRATOS_ERROR << "Calling virtual function for CalculateResponseKirchhoff" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::CalculateResponseCauchy(
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues
    )
{
    KRATOS_ERROR << "Calling virtual function for CalculateResponseCauchy" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::InitializeResponse (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues,
    const ConstitutiveLaw::StressMeasure& rStressMeasure
    )
{
    switch(rStressMeasure)
    {
    case ConstitutiveLaw::StressMeasure_PK1:
        InitializeResponsePK1(pConstitutiveLaw, rParameterValues);
        break;

    case ConstitutiveLaw::StressMeasure_PK2:
        InitializeResponsePK2(pConstitutiveLaw, rParameterValues);
        break;

    case ConstitutiveLaw::StressMeasure_Kirchhoff:
        InitializeResponseKirchhoff(pConstitutiveLaw, rParameterValues);
        break;

    case ConstitutiveLaw::StressMeasure_Cauchy:
        InitializeResponseCauchy(pConstitutiveLaw, rParameterValues);
        break;

    default:
        KRATOS_ERROR << " Stress Measure not Defined "<< std::endl;;
        break;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::InitializeResponsePK1 (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues
    )
{
    // NOT MANDATORY
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::InitializeResponsePK2 (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues
    )
{
    // NOT MANDATORY
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::InitializeResponseKirchhoff (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues
    )
{
    // NOT MANDATORY
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::InitializeResponseCauchy (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues
    )
{
    // NOT MANDATORY
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::FinalizeResponse (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues,
    const ConstitutiveLaw::StressMeasure& rStressMeasure
    )
{
    switch(rStressMeasure)
    {
    case ConstitutiveLaw::StressMeasure_PK1:
        FinalizeResponsePK1(pConstitutiveLaw, rParameterValues);
        break;

    case ConstitutiveLaw::StressMeasure_PK2:
        FinalizeResponsePK2(pConstitutiveLaw, rParameterValues);
        break;

    case ConstitutiveLaw::StressMeasure_Kirchhoff:
        FinalizeResponseKirchhoff(pConstitutiveLaw, rParameterValues);
        break;

    case ConstitutiveLaw::StressMeasure_Cauchy:
        FinalizeResponseCauchy(pConstitutiveLaw, rParameterValues);
        break;

    default:
        KRATOS_ERROR << " Stress Measure not Defined "<< std::endl;;
        break;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::FinalizeResponsePK1 (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues
    )
{
    // NOT MANDATORY
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::FinalizeResponsePK2(
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues
    )
{
    // NOT MANDATORY
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::FinalizeResponseKirchhoff (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues
    )
{
    // NOT MANDATORY
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::FinalizeResponseCauchy (
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rParameterValues
    )
{
    // NOT MANDATORY
}

/***********************************************************************************/
/***********************************************************************************/

void ImposedDeformation::Reset(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // NOT MANDATORY
}

/***********************************************************************************/
/***********************************************************************************/

int ImposedDeformation::Check(
    const ConstitutiveLaw* pConstitutiveLaw,
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    return 0;
}

} // namespace Kratos.
