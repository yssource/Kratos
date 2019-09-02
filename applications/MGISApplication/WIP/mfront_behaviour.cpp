// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_constitutive/mfront_behaviour.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos {

  /******************************CONSTRUCTOR******************************************/
  /***********************************************************************************/

  MFrontBehaviourBase::MFrontBehaviourBase(const Hypothesis h,
                                           const std::string& l,
                                           const std::string& b)
      : ConstitutiveLaw() {
    this->behaviour =
        Kratos::make_shared<Behaviour>(mgis::behaviour::load(l, b, h));
    if ((this->behaviour->btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) &&
        (this->behaviour->kinematic == Behaviour::SMALLSTRAINKINEMATIC)) {
      this->strain_measure = StrainMeasure_Infinitesimal;
      this->stress_measure = StressMeasure_Cauchy;
    } else if ((this->behaviour->btype ==
                Behaviour::FINITESTRAINBASEDBEHAVIOUR) &&
               (thi->behaviour->kinematix ==
                Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY)) {
      this->strain_measure = StrainMeasure_Deformation_Gradient;
      this->stress_measure = StressMeasure_Cauchy;
    } else {
      KRATOS_ERROR() << "unsupported behaviour " << std::endl;
    }
    this->data = BehaviourData(*(this->behaviour));
  }

  // // TO BE DONE
  //   MFrontBehaviourBase::MFrontBehaviourBase(const
  //   FiniteStrainBehaviourOptions& o,
  //                                            const Hypothesis h,
  //                                            const std::string& l,
  //                                            const std::string& b)
  //       : ConstitutiveLaw() {
  //     this->behaviour =
  //         Kratos::make_shared<Behaviour>(mgis::behaviour::load(o, l, b, h));
  //     this->data = BehaviourData(*(this->behaviour));
  //   }

  /******************************COPY
   * CONSTRUCTOR*************************************/
  /***********************************************************************************/

  MFrontBehaviourBase::MFrontBehaviourBase(const MFrontBehaviourBase& rOther)
      : ConstitutiveLaw(rOther) = default;

  /********************************CLONE**********************************************/
  /***********************************************************************************/

  ConstitutiveLaw::Pointer MFrontBehaviourBase::Clone() const {
    return Kratos::make_shared<MFrontBehaviourBase>(*this);
  }

  /*******************************DESTRUCTOR******************************************/
  /***********************************************************************************/

  MFrontBehaviourBase::~MFrontBehaviourBase() = default;

  /***********************************************************************************/
  /***********************************************************************************/

  SizeType MFrontBehaviourBase::WorkingSpaceDimension() override {
    return mgis::behaviour::getSpaceDimension(this->behaviour->hypothesis);
  }  // end of MFrontBehaviourBase::WorkingSpaceDimension

  SizeType MFrontBehaviourBase::GetStrainSize() override {
    if (this->strain_measure == StrainMeasure_Deformation_Gradient) {
      return mgis::behaviour::getTensorSize(this->behaviour->hypothesis);
    } else if (this->strain_measure == StrainMeasure_Infinitesimal) {
      return mgis::behaviour::getStensorSize(this->behaviour->hypothesis);
    }
    KRATOS_ERROR() << "unsupported behaviour " << std::endl;
    return 0;
  } // end of MFrontBehaviourBase::GetStrainSize

  StrainMeasure MFrontBehaviourBase::GetStrainMeasure() override {
    return this->strain_measure;
  }

  StressMeasure MFrontBehaviourBase::GetStressMeasure() override {
    return this->stress_measure;
  }

  void MFrontBehaviourBase::CalculateMaterialResponsePK2(
      ConstitutiveLaw::Parameters& rValues) {
    KRATOS_TRY;
    //     // b.- Get Values to compute the constitutive law:
    //     Flags& r_constitutive_law_options = rValues.GetOptions();
    //
    //     Vector& r_strain_vector = rValues.GetStrainVector();
    //
    //     // NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN
    //     // MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    //     if (r_constitutive_law_options.IsNot(
    //             ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //       CalculateCauchyGreenStrain(rValues, r_strain_vector);
    //     }
    //
    //     if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
    //       Vector& r_stress_vector = rValues.GetStressVector();
    //       CalculatePK2Stress(r_strain_vector, r_stress_vector, rValues);
    //     }
    //
    //     if (r_constitutive_law_options.Is(
    //             ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
    //       Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
    //       CalculateElasticMatrix(r_constitutive_matrix, rValues);
    //     }
    //
    KRATOS_CATCH("");
  }

  /***********************************************************************************/
  /***********************************************************************************/

  // NOTE: Since we are in the hypothesis of small strains we can use the same
  // function for everything

  void MFrontBehaviourBase::CalculateMaterialResponsePK1(
      ConstitutiveLaw::Parameters& rValues) {

  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MFrontBehaviourBase::CalculateMaterialResponseKirchhoff(
      ConstitutiveLaw::Parameters& rValues) {
//    CalculateMaterialResponsePK2(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MFrontBehaviourBase::CalculateMaterialResponseCauchy(
      ConstitutiveLaw::Parameters& rValues) {
//    CalculateMaterialResponsePK2(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MFrontBehaviourBase::InitializeMaterialResponsePK1(
      ConstitutiveLaw::Parameters& rValues) {
    //     // Small deformation so we can call the Cauchy method
    //     InitializeMaterialResponseCauchy(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MFrontBehaviourBase::InitializeMaterialResponsePK2(
      ConstitutiveLaw::Parameters& rValues) {
    //     // Small deformation so we can call the Cauchy method
    //     InitializeMaterialResponseCauchy(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MFrontBehaviourBase::InitializeMaterialResponseCauchy(
      ConstitutiveLaw::Parameters& rValues) {
    // TODO: Add if necessary
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MFrontBehaviourBase::InitializeMaterialResponseKirchhoff(
      ConstitutiveLaw::Parameters& rValues) {
    // Small deformation so we can call the Cauchy method
//    InitializeMaterialResponseCauchy(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MFrontBehaviourBase::FinalizeMaterialResponsePK1(
      ConstitutiveLaw::Parameters& rValues) {
    // Small deformation so we can call the Cauchy method
//    FinalizeMaterialResponseCauchy(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MFrontBehaviourBase::FinalizeMaterialResponsePK2(
      ConstitutiveLaw::Parameters& rValues) {
    // Small deformation so we can call the Cauchy method
//    FinalizeMaterialResponseCauchy(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MFrontBehaviourBase::FinalizeMaterialResponseCauchy(
      ConstitutiveLaw::Parameters& rValues) {
    // TODO: Add if necessary
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MFrontBehaviourBase::FinalizeMaterialResponseKirchhoff(
      ConstitutiveLaw::Parameters& rValues) {
    // Small deformation so we can call the Cauchy method
//    FinalizeMaterialResponseCauchy(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  double& MFrontBehaviourBase::CalculateValue(
      ConstitutiveLaw::Parameters& rParameterValues,
      const Variable<double>& rThisVariable,
      double& rValue) {
    //     Vector& r_strain_vector = rParameterValues.GetStrainVector();
    //     Vector& r_stress_vector = rParameterValues.GetStressVector();
    //
    //     if (rThisVariable == STRAIN_ENERGY) {
    //       this->CalculateCauchyGreenStrain(rParameterValues,
    //       r_strain_vector);
    //       this->CalculatePK2Stress(r_strain_vector, r_stress_vector,
    //                                rParameterValues);
    //
    //       rValue = 0.5 * inner_prod(r_strain_vector,
    //                                 r_stress_vector);  // Strain energy =
    //                                 0.5*E:C:E
    //     }
    //
    return (rValue);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  Vector& MFrontBehaviourBase::CalculateValue(
      ConstitutiveLaw::Parameters& rParameterValues,
      const Variable<Vector>& rThisVariable,
      Vector& rValue) {
    if (rThisVariable == STRAIN ||
        rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {
      this->CalculateCauchyGreenStrain(rParameterValues, rValue);
    } else if (rThisVariable == STRESSES ||
               rThisVariable == CAUCHY_STRESS_VECTOR ||
               rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
               rThisVariable == PK2_STRESS_VECTOR) {
      // Get Values to compute the constitutive law:
      Flags& r_flags = rParameterValues.GetOptions();

      // Previous flags saved
      const bool flag_const_tensor =
          r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
      const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

      r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
      r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

      // We compute the stress
      MFrontBehaviourBase::CalculateMaterialResponseCauchy(rParameterValues);
      rValue = rParameterValues.GetStressVector();

      // Previous flags restored
      r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR,
                  flag_const_tensor);
      r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    }

    return (rValue);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  Matrix& MFrontBehaviourBase::CalculateValue(
      ConstitutiveLaw::Parameters& rParameterValues,
      const Variable<Matrix>& rThisVariable,
      Matrix& rValue) {
    //     if (rThisVariable == CONSTITUTIVE_MATRIX ||
    //         rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
    //         rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
    //       this->CalculateElasticMatrix(rValue, rParameterValues);
    //     }

    return (rValue);
  }

  //*************************CONSTITUTIVE LAW GENERAL FEATURES
  //*************************
  /***********************************************************************************/

  void MFrontBehaviourBase::GetLawFeatures(Features& rFeatures) {
    // Set the type of law
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    if (this->behaviour->symmetry == Behaviour::ISOTROPIC) {
      rFeatures.mOptions.Set(ISOTROPIC);
    } else if (this->behaviour->symmetry == Behaviour::ORTHOTROPIC) {
      rFeatures.mOptions.Set(ORTHOTROPIC);
    }
    // Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(this->GetStrainMeasure());
    // Set the strain size
    rFeatures.mStrainSize = this->GetStrainSize();
    // Set the spacedimension
    rFeatures.mSpaceDimension = this->WorkingSpaceDimension;
  }

  /***********************************************************************************/
  /***********************************************************************************/

  int MFrontBehaviourBase::Check(const Properties& rMaterialProperties,
                             const GeometryType& rElementGeometry,
                             const ProcessInfo& rCurrentProcessInfo) {

    return 0;
  }

}  // Namespace Kratos
