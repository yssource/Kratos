//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/borja_hencky_cam_clay_plane_strain_2D_law.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   BorjaHenckyCamClayPlasticPlaneStrain2DLaw::BorjaHenckyCamClayPlasticPlaneStrain2DLaw()
      : NonLinearHenckyElasticPlasticPlaneStrain2DLaw()
   {
      mpHardeningLaw   = HardeningLaw::Pointer( new CamClayHardeningLaw() );
      mpYieldCriterion = YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
      mpFlowRule       = FlowRule::Pointer( new BorjaCamClayExplicitFlowRule(mpYieldCriterion) );
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   BorjaHenckyCamClayPlasticPlaneStrain2DLaw::BorjaHenckyCamClayPlasticPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
   {
      mpHardeningLaw    =  pHardeningLaw;
      mpYieldCriterion  =  YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
      mpFlowRule        =  pFlowRule;
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   BorjaHenckyCamClayPlasticPlaneStrain2DLaw::BorjaHenckyCamClayPlasticPlaneStrain2DLaw(const BorjaHenckyCamClayPlasticPlaneStrain2DLaw& rOther)
      : NonLinearHenckyElasticPlasticPlaneStrain2DLaw(rOther)
   {

   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveLaw::Pointer BorjaHenckyCamClayPlasticPlaneStrain2DLaw::Clone() const
   {
      BorjaHenckyCamClayPlasticPlaneStrain2DLaw::Pointer p_clone(new BorjaHenckyCamClayPlasticPlaneStrain2DLaw(*this));
      return p_clone;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   BorjaHenckyCamClayPlasticPlaneStrain2DLaw::~BorjaHenckyCamClayPlasticPlaneStrain2DLaw()
   {
   }


   double& BorjaHenckyCamClayPlasticPlaneStrain2DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue)
   {
      /*      if ( rThisVariable == PLASTIC_STRESS_LIKE)
              {
              rValue = this->GetValue( PRECONSOLIDATION, rValue);
              }
              else if ( rThisVariable == PLASTIC_STRESS_LIKE)
              {
              rValue = this->GetValue( VOLUMETRIC_PLASTIC, rValue);
              }*/
      if ( rThisVariable == PENALTY_PARAMETER)
      {
      }
      else if ( rThisVariable == M_MODULUS ) {
         double Swelling = mpYieldCriterion->GetHardeningLaw().GetProperties()[ SWELLING_SLOPE ];
         double MeanStress ;
         MeanStress = this->GetValue( STRESS_INV_P, MeanStress);
         double K = MeanStress/Swelling;

         double Alpha = mpYieldCriterion->GetHardeningLaw().GetProperties()[ ALPHA_SHEAR ];
         double G = mpYieldCriterion->GetHardeningLaw().GetProperties()[ INITIAL_SHEAR_MODULUS ];
         G += Alpha*MeanStress;  // this modulus is approximated (not the real one)

         rValue = K + 4.0*G / 3.0;
      }
      else if ( ( rThisVariable == YOUNG_MODULUS) || ( rThisVariable == EQUIVALENT_YOUNG_MODULUS ) || ( rThisVariable == SHEAR_MODULUS ) || (rThisVariable == BULK_MODULUS ) )
      {
         double Swelling = mpYieldCriterion->GetHardeningLaw().GetProperties()[ SWELLING_SLOPE ];
         double MeanStress ;
         MeanStress = this->GetValue( STRESS_INV_P, MeanStress);
         double K = MeanStress/Swelling;
         if ( rThisVariable == BULK_MODULUS)
         {
            rValue = K;
            return rValue;
         }

         double Alpha = mpYieldCriterion->GetHardeningLaw().GetProperties()[ ALPHA_SHEAR ];
         double G = mpYieldCriterion->GetHardeningLaw().GetProperties()[ INITIAL_SHEAR_MODULUS ];
         G += Alpha*MeanStress;  // this modulus is approximated (not the real one)

         if ( rThisVariable == SHEAR_MODULUS)
         {
            rValue = G;
            return rValue;
         }

         rValue = 9.0*K*G / ( 3.0*K + G);

      }
      else
      {
         rValue = NonLinearHenckyElasticPlasticPlaneStrain2DLaw::GetValue( rThisVariable, rValue);
      }

      return rValue;
   }

   void BorjaHenckyCamClayPlasticPlaneStrain2DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
   {
      /*      if ( rThisVariable == PLASTIC_STRESS_LIKE )
              {
              double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
              const double OtherSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];
              double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
              double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
              ReferencePressure *= OCR;

              double PreconsolidationStress = rValue;
              double VolumetricPlasticDef = - (OtherSlope - SwellingSlope) * std::log(-PreconsolidationStress / ReferencePressure);
              FlowRule::RadialReturnVariables ReturnMappingVariables;
              ReturnMappingVariables.DeltaGamma = VolumetricPlasticDef;
              mpFlowRule->UpdateInternalVariables( ReturnMappingVariables);
              }
              else if ( rThisVariable == PLASTIC_STRAIN_LIKE)
              {
              FlowRule::RadialReturnVariables ReturnMappingVariables;
              ReturnMappingVariables.DeltaGamma = rValue;
              mpFlowRule->UpdateInternalVariables( ReturnMappingVariables );
              }*/
      if ( rThisVariable == PENALTY_PARAMETER)
      {
      }
      else 
      {
         NonLinearHenckyElasticPlasticPlaneStrain2DLaw::SetValue( rThisVariable, rValue, rCurrentProcessInfo );
      }
   }

   void BorjaHenckyCamClayPlasticPlaneStrain2DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
   {
         NonLinearHenckyElasticPlasticPlaneStrain2DLaw::SetValue( rThisVariable, rValue, rCurrentProcessInfo );

   }


   int BorjaHenckyCamClayPlasticPlaneStrain2DLaw::Check( const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
   {

      return 0;
   }


} // Namespace Kratos
