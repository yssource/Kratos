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
#include "custom_utilities/tangent_operator_calculator_utility.h"

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
        this->CalculatePseudoElasticMatrix(r_constitutive_matrix, rValues, mMartensitePercentage);
    }

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculatePseudoElasticMatrix(r_constitutive_matrix, rValues, mMartensitePercentage);

        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }

        // S0 = C:(E-Ep)
        array_1d<double, VoigtSize> predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector - mTransformationStrain);
        array_1d<double, VoigtSize> stress_deviator;
        bool is_loading;
        this->CheckIfLoading(r_constitutive_matrix, r_strain_vector, is_loading);
        const double uniaxial_stress = this->CalculatePseudoDruckerPragerUniaxialStress(predictive_stress_vector, rValues, stress_deviator);
        const double threshold = this->CalculateThreshold(rValues, is_loading);
        const double yield_condition = uniaxial_stress - threshold;
        bool save_internal_vars = false;
        this->IntegrateStressVector(predictive_stress_vector, yield_condition, rValues, save_internal_vars, is_loading, stress_deviator);
        r_integrated_stress_vector = predictive_stress_vector;

        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            this->CalculateTangentTensor(rValues); // this modifies the ConstitutiveMatrix
            const bool check_analytical_tensor = false;
            Matrix analytical_tensor(VoigtSize, VoigtSize);
            if (check_analytical_tensor)
                this->CalculateAnalyticalTensor(rValues, is_loading, stress_deviator, yield_condition, r_constitutive_matrix, analytical_tensor);
            // KRATOS_WATCH(analytical_tensor)
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::CalculateAnalyticalTensor(
    ConstitutiveLaw::Parameters& rValues,
    const bool IsLoading,
    const array_1d<double, VoigtSize>& rDeviator,
    const double YieldCondition,
    const Matrix& rPseudoElasticMatrix,
    Matrix& rAnalyticalTensor
    )
{
    if (rAnalyticalTensor.size1() != VoigtSize)
        rAnalyticalTensor.resize(VoigtSize, VoigtSize);
    array_1d<double, VoigtSize> identity_vector;
    ConstitutiveLawUtilities<VoigtSize>::CalculateFirstVector(identity_vector);

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double Ea = r_material_properties[AUSTENITIC_YOUNG_MODULUS];
    const double Em = r_material_properties[MARTENSITIC_YOUNG_MODULUS];
    const double E  = Ea * Em / (Em + (mMartensitePercentage * (Ea - Em)));
    const double NU = r_material_properties[POISSON_RATIO];
    const double bulk_mod  = E / (3.0 * (1.0 - (2 * NU)));
    const double shear_mod = E / (2.0 * (1.0 + NU));
    const double max_residual_strain = r_material_properties[MAXIMUN_RESIDUAL_STRAIN];
    const double yield_start_forward = r_material_properties[YIELD_STRESS_START_FORWARD_TRANSFORMATION];
    const double yield_compression_start_forward = r_material_properties[YIELD_STRESS_COMPRESSION_START_FORWARD_TRANSFORMATION];
    const double alpha = (yield_compression_start_forward - yield_start_forward) / (yield_start_forward + yield_compression_start_forward);
    const Matrix identity_matrix = IdentityMatrix(VoigtSize);
    Matrix aux_correction(VoigtSize, VoigtSize);

    if (IsLoading) {
        const double yield_end_forward   = r_material_properties[YIELD_STRESS_END_FORWARD_TRANSFORMATION];
        const double slope_forward = (yield_end_forward - yield_start_forward) / max_residual_strain;
        const double transformation_consistency_factor = YieldCondition / (9.0 * std::pow(alpha, 2) * bulk_mod + 3.0 * shear_mod + slope_forward);
        const Vector aux_vector = std::sqrt(1.5) * rDeviator / norm_2(rDeviator) + alpha * identity_vector;
        const Vector flow_vector = aux_vector / norm_2(aux_vector);
        const double aux_D = (9.0 * std::pow(alpha, 2) * bulk_mod + 3.0 * shear_mod + slope_forward);
        noalias(aux_correction) = (identity_matrix - (1.0 / 3.0) * outer_prod(identity_vector, identity_vector)) - outer_prod(flow_vector, flow_vector);
        const Matrix correction = transformation_consistency_factor * 4.0 * std::pow(shear_mod, 2) * std::sqrt(1.5) * (1.0 / norm_2(rDeviator)) * aux_correction;
        const Vector CE_derivF = alpha * bulk_mod * identity_vector + 2.0 * shear_mod * std::sqrt(1.5) * flow_vector;
        noalias(rAnalyticalTensor) = rPseudoElasticMatrix - (1 / aux_D) * outer_prod(CE_derivF, CE_derivF) - correction;
    } else {
        const double yield_end_backward   = r_material_properties[YIELD_STRESS_END_BACKWARDS_TRANSFORMATION];
        const double yield_start_backward   = r_material_properties[YIELD_STRESS_START_BACKWARDS_TRANSFORMATION];
        const double slope_backward = (yield_start_backward - yield_end_backward) / max_residual_strain;
        const double transformation_consistency_factor = YieldCondition / (9.0 * std::pow(alpha, 2) * bulk_mod + 3.0 * shear_mod + slope_backward);
        const Vector flow_vector = mTransformationStrain / norm_2(mTransformationStrain);
        const double aux_D = (9.0 * std::pow(alpha, 2) * bulk_mod + 3.0 * shear_mod + slope_backward);
        noalias(aux_correction) = (identity_matrix - (1.0 / 3.0) * outer_prod(identity_vector, identity_vector)) - outer_prod(flow_vector, flow_vector);
        const Matrix correction = transformation_consistency_factor * 4.0 * std::pow(shear_mod, 2) * std::sqrt(1.5) * (1.0 / norm_2(rDeviator)) * aux_correction;
        const Vector CE_derivF = alpha * bulk_mod * identity_vector + 2.0 * shear_mod * std::sqrt(1.5) * flow_vector;
        noalias(rAnalyticalTensor) = rPseudoElasticMatrix - (1 / aux_D) * outer_prod(CE_derivF, CE_derivF) - correction;
    }
}


/***********************************************************************************/
/***********************************************************************************/
template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        KRATOS_ERROR << "Analytic solution not available" << std::endl;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::IntegrateStressVector(
    array_1d<double, VoigtSize>& rStressVector,
    const double YieldCondition,
    ConstitutiveLaw::Parameters& rValues,
    const bool SaveInternalVars,
    const bool IsLoading,
    const array_1d<double, VoigtSize>& rDeviator
    )
{
    if (IsLoading) {
        if (YieldCondition <= tolerance && mMartensitePercentage <= tolerance) {
            return;
        } else if (YieldCondition > tolerance && mMartensitePercentage >= 0.99) {
            return;
        } else if (YieldCondition > tolerance && mMartensitePercentage >= tolerance && mMartensitePercentage < 0.99) {
            this->ForwardTransformation(YieldCondition, rValues, rStressVector, rDeviator, SaveInternalVars);
        }
    } else { // unloading
        if (YieldCondition >= tolerance && mMartensitePercentage >= 0.99) {
            return;
        } else if (YieldCondition < tolerance && mMartensitePercentage <= tolerance) {
            return;
        } else if (YieldCondition < tolerance && mMartensitePercentage > tolerance && mMartensitePercentage <= 0.99) {
            this->BackwardTransformation(YieldCondition, rValues, rStressVector, rDeviator, SaveInternalVars);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::ForwardTransformation(
    const double YieldCondition,
    ConstitutiveLaw::Parameters& rValues,
    array_1d<double, VoigtSize>& rStressVector,
    const array_1d<double, VoigtSize>& rDeviator,
    const bool SaveInternalVars
    )
{
    array_1d<double, VoigtSize> identity_vector;
    ConstitutiveLawUtilities<VoigtSize>::CalculateFirstVector(identity_vector);
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double Ea = r_material_properties[AUSTENITIC_YOUNG_MODULUS];
    const double Em = r_material_properties[MARTENSITIC_YOUNG_MODULUS];
    const double E  = Ea * Em / (Em + (mMartensitePercentage * (Ea - Em)));
    const double NU = r_material_properties[POISSON_RATIO];
    const double bulk_mod  = E / (3.0 * (1.0 - (2 * NU)));
    const double shear_mod = E / (2.0 * (1.0 + NU));
    const double max_residual_strain = r_material_properties[MAXIMUN_RESIDUAL_STRAIN];
    const double yield_start_forward = r_material_properties[YIELD_STRESS_START_FORWARD_TRANSFORMATION];
    const double yield_end_forward   = r_material_properties[YIELD_STRESS_END_FORWARD_TRANSFORMATION];
    const double slope_forward = (yield_end_forward - yield_start_forward) / max_residual_strain;
    const double yield_compression_start_forward = r_material_properties[YIELD_STRESS_COMPRESSION_START_FORWARD_TRANSFORMATION];
    const double alpha = (yield_compression_start_forward - yield_start_forward) / (yield_start_forward + yield_compression_start_forward);

    const double transformation_consistency_factor = YieldCondition / (9.0 * std::pow(alpha, 2) * bulk_mod + 3.0 * shear_mod + slope_forward);
    const Vector aux_vector = std::sqrt(1.5) * rDeviator / norm_2(rDeviator) + alpha * identity_vector;
    const Vector flow_vector = aux_vector / norm_2(aux_vector);

    double updated_martensite_percentage = mMartensitePercentage + transformation_consistency_factor / max_residual_strain;
    updated_martensite_percentage  = (updated_martensite_percentage >= 1.0) ? 1.0 : updated_martensite_percentage;
    const Vector updated_transformation_strain = mTransformationStrain + transformation_consistency_factor * std::sqrt(1.5) * flow_vector;

    if (SaveInternalVars) {
        mMartensitePercentage = updated_martensite_percentage;
        mTransformationStrain = updated_transformation_strain;
    }

    Matrix updated_pseudo_elastic_matrix;
    this->CalculatePseudoElasticMatrix(updated_pseudo_elastic_matrix, rValues, updated_martensite_percentage);
    rStressVector -= prod(updated_pseudo_elastic_matrix, updated_transformation_strain);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::BackwardTransformation(
    const double YieldCondition,
    ConstitutiveLaw::Parameters& rValues,
    array_1d<double, VoigtSize>& rStressVector,
    const array_1d<double, VoigtSize>& rDeviator,
    const bool SaveInternalVars
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double Ea = r_material_properties[AUSTENITIC_YOUNG_MODULUS];
    const double Em = r_material_properties[MARTENSITIC_YOUNG_MODULUS];
    const double E  = Ea * Em / (Em + (mMartensitePercentage * (Ea - Em)));
    const double NU = r_material_properties[POISSON_RATIO];
    const double bulk_mod  = E / (3.0 * (1.0 - (2 * NU)));
    const double shear_mod = E / (2.0 * (1.0 + NU));
    const double yield_start_forward   = r_material_properties[YIELD_STRESS_START_FORWARD_TRANSFORMATION];
    const double max_residual_strain   = r_material_properties[MAXIMUN_RESIDUAL_STRAIN];
    const double yield_start_backwards = r_material_properties[YIELD_STRESS_START_BACKWARDS_TRANSFORMATION];
    const double yield_end_backwards   = r_material_properties[YIELD_STRESS_END_BACKWARDS_TRANSFORMATION];
    const double slope_backwards = (yield_start_backwards - yield_end_backwards) / max_residual_strain;
    const double yield_compression_start_forward = r_material_properties[YIELD_STRESS_COMPRESSION_START_FORWARD_TRANSFORMATION];
    const double alpha = (yield_compression_start_forward - yield_start_forward) / (yield_start_forward + yield_compression_start_forward);

    const double transformation_consistency_factor = YieldCondition / (9.0 * std::pow(alpha, 2) * bulk_mod + 3.0 * shear_mod + slope_backwards);
    double updated_martensite_percentage = mMartensitePercentage + transformation_consistency_factor / max_residual_strain;
    updated_martensite_percentage  = (updated_martensite_percentage <= tolerance) ? 0.0 : updated_martensite_percentage;
    const Vector flow_vector = mTransformationStrain / norm_2(mTransformationStrain);
    const Vector updated_transformation_strain = mTransformationStrain + transformation_consistency_factor * std::sqrt(1.5) * flow_vector;

    if (SaveInternalVars) {
        mMartensitePercentage = updated_martensite_percentage;
        mTransformationStrain = updated_transformation_strain;
    }

    Matrix updated_pseudo_elastic_matrix;
    this->CalculatePseudoElasticMatrix(updated_pseudo_elastic_matrix, rValues, updated_martensite_percentage);
    rStressVector -= prod(updated_pseudo_elastic_matrix, updated_transformation_strain);
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
        this->CalculatePseudoElasticMatrix(r_constitutive_matrix, rValues, mMartensitePercentage);
    }

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculatePseudoElasticMatrix(r_constitutive_matrix, rValues, mMartensitePercentage);

        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }

        // S0 = C:(E-Ep)
        array_1d<double, VoigtSize> predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector - mTransformationStrain);
        array_1d<double, VoigtSize> stress_deviator;
        bool is_loading;
        this->CheckIfLoading(r_constitutive_matrix, r_strain_vector, is_loading);
        const double uniaxial_stress = this->CalculatePseudoDruckerPragerUniaxialStress(predictive_stress_vector, rValues, stress_deviator);
        const double threshold = this->CalculateThreshold(rValues, is_loading);
        const double yield_condition = uniaxial_stress - threshold;
        bool save_internal_vars = true;
        this->IntegrateStressVector(predictive_stress_vector, yield_condition, rValues, save_internal_vars, is_loading, stress_deviator);
        r_integrated_stress_vector = predictive_stress_vector;
        mPreviousStrain = r_strain_vector;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
bool NitinolPseudoElasticity3D<TElasticBehaviourLaw>::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == MARTENSITE_PERCENTAGE) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
bool NitinolPseudoElasticity3D<TElasticBehaviourLaw>::Has(const Variable<Vector>& rThisVariable)
{
    if (rThisVariable == TRANSFORMATION_STRAIN) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }
    return BaseType::Has(rThisVariable);
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
    if (rThisVariable == TRANSFORMATION_STRAIN) {
        rValue = mTransformationStrain;
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }
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
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
double& NitinolPseudoElasticity3D<TElasticBehaviourLaw>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == MARTENSITE_PERCENTAGE) {
        return mMartensitePercentage;
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }
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
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::CalculatePseudoElasticMatrix(
    Matrix& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues,
    const double MartensitePercentage
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double Ea = r_material_properties[AUSTENITIC_YOUNG_MODULUS];
    const double Em = r_material_properties[MARTENSITIC_YOUNG_MODULUS];
    const double E = Ea * Em / (Em + (MartensitePercentage * (Ea - Em)));
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

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void NitinolPseudoElasticity3D<TElasticBehaviourLaw>::CheckIfLoading(
    const Matrix& rPseudoElasticMatrix, 
    const Vector& rStrainVector,
    bool& rIsLoading
    )
{
    const double aux = norm_2(prod(rPseudoElasticMatrix, rStrainVector)) - norm_2(prod(rPseudoElasticMatrix, mPreviousStrain));
    if (aux < 0.0)
        rIsLoading = false;
    else
        rIsLoading = true; // Todo what if == 0.0
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
double NitinolPseudoElasticity3D<TElasticBehaviourLaw>::CalculatePseudoDruckerPragerUniaxialStress(
    const array_1d<double, VoigtSize>& rStressVector,
    ConstitutiveLaw::Parameters& rValues,
    array_1d<double, VoigtSize>& rDeviator
)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    double I1, J2;
    array_1d<double, VoigtSize> stress_deviator;
    ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant(rStressVector, I1);
    ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rStressVector, I1, rDeviator, J2);

    const double yield_compression_start_forward = r_material_properties[YIELD_STRESS_COMPRESSION_START_FORWARD_TRANSFORMATION];
    const double yield_start_forward = r_material_properties[YIELD_STRESS_START_FORWARD_TRANSFORMATION];
    const double alpha = (yield_compression_start_forward - yield_start_forward) / (yield_start_forward + yield_compression_start_forward);

    return std::sqrt(3.0 * J2) + alpha * I1;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
double NitinolPseudoElasticity3D<TElasticBehaviourLaw>::CalculateThreshold(
    ConstitutiveLaw::Parameters& rValues,
    const bool IsLoading
)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double max_residual_strain = r_material_properties[MAXIMUN_RESIDUAL_STRAIN];
    if (IsLoading) {
        const double yield_start_forward = r_material_properties[YIELD_STRESS_START_FORWARD_TRANSFORMATION];
        const double yield_end_forward   = r_material_properties[YIELD_STRESS_END_FORWARD_TRANSFORMATION];
        const double slope_forward = (yield_end_forward - yield_start_forward) / max_residual_strain;
        return yield_start_forward + slope_forward * max_residual_strain * mMartensitePercentage;
    } else {
        const double yield_start_backwards = r_material_properties[YIELD_STRESS_START_BACKWARDS_TRANSFORMATION];
        const double yield_end_backwards   = r_material_properties[YIELD_STRESS_END_BACKWARDS_TRANSFORMATION];
        const double slope_backward = (yield_start_backwards - yield_end_backwards) / max_residual_strain;
        return yield_end_backwards + slope_backward * max_residual_strain * mMartensitePercentage;
    }
}
/***********************************************************************************/
/***********************************************************************************/

template class NitinolPseudoElasticity3D<ElasticIsotropic3D>;

} // namespace Kratos
