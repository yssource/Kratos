// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Agustina Giuliodori
//  Collaborator:
//

// System includes

// Project includes
#include "custom_constitutive/nitinol_pseudo_elasticity_3d.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
NitinolPseudoElasticity3D<TElasticBehaviourLaw>::NitinolPseudoElasticity3D()
    : BaseType()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
NitinolPseudoElasticity3D<TElasticBehaviourLaw>::NitinolPseudoElasticity3D(const NitinolPseudoElasticity3D& rOther)
    : BaseType(rOther)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
ConstitutiveLaw::Pointer NitinolPseudoElasticity3D<TElasticBehaviourLaw>::Clone() const
{
    return Kratos::make_shared<NitinolPseudoElasticity3D>(*this);
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
NitinolPseudoElasticity3D<TElasticBehaviourLaw>::~NitinolPseudoElasticity3D()
{
}

/***********************************************************************************/
/***********************************************************************************/
    
template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::
CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::
CalculateMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::
CalculateMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{








}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::
FinalizeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::
FinalizeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::
FinalizeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::
FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{




}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
Vector& NitinolPseudoElasticity3D<TElasticBehaviourLaw>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
Matrix& NitinolPseudoElasticity3D<TElasticBehaviourLaw>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    // if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
    //     rValue = MathUtils<double>::StressVectorToTensor(this->GetPreviousStressVector());
    // } else {
    //     rValue = BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    // }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
int NitinolPseudoElasticity3D<TElasticBehaviourLaw>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(VoigtSize == this->GetStrainSize()) << "You are combining not compatible constitutive laws" << std::endl;
    return check_base;
}

/***********************************************************************************/
/***********************************************************************************/

template class NitinolPseudoElasticity3D<ElasticIsotropic3D>;

} // namespace Kratos
