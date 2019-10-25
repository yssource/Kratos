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
    Vector& r_integrated_stress_vector = rValues.GetStressVector();
    Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix();
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // Elastic Matrix C0
    Matrix r_constitutive_matrix;
    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        this->CalculatePseudoElasticMatrix(r_constitutive_matrix);
    }

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

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D::CalculatePseudoElasticMatrix(
    Matrix& rConstitutiveMatrix
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    this->CheckClearElasticMatrix(rConstitutiveMatrix);

    const double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
    const double c2 = c1 * ( 1 - NU );
    const double c3 = c1 * NU;
    const double c4 = c1 * 0.5 * ( 1 - 2 * NU );

    rConstitutiveMatrix( 0, 0 ) = c2;
    rConstitutiveMatrix( 0, 1 ) = c3;
    rConstitutiveMatrix( 0, 2 ) = c3;
    rConstitutiveMatrix( 1, 0 ) = c3;
    rConstitutiveMatrix( 1, 1 ) = c2;
    rConstitutiveMatrix( 1, 2 ) = c3;
    rConstitutiveMatrix( 2, 0 ) = c3;
    rConstitutiveMatrix( 2, 1 ) = c3;
    rConstitutiveMatrix( 2, 2 ) = c2;
    rConstitutiveMatrix( 3, 3 ) = c4;
    rConstitutiveMatrix( 4, 4 ) = c4;
    rConstitutiveMatrix( 5, 5 ) = c4;
}


template class NitinolPseudoElasticity3D<ElasticIsotropic3D>;

} // namespace Kratos
