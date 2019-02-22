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
#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "custom_constitutive/borja_hencky_cam_clay_axisym_2D_law.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   BorjaHenckyCamClayPlasticAxisym2DLaw::BorjaHenckyCamClayPlasticAxisym2DLaw()
      : NonLinearHenckyElasticPlasticAxisym2DLaw()
   {
      mpHardeningLaw   = HardeningLaw::Pointer( new CamClayHardeningLaw() );
      mpYieldCriterion = YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
      mpFlowRule       = FlowRule::Pointer( new BorjaCamClayExplicitFlowRule(mpYieldCriterion) );
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   BorjaHenckyCamClayPlasticAxisym2DLaw::BorjaHenckyCamClayPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
   {
      mpHardeningLaw    =  pHardeningLaw;
      mpYieldCriterion  =  YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
      mpFlowRule        =  pFlowRule;
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   BorjaHenckyCamClayPlasticAxisym2DLaw::BorjaHenckyCamClayPlasticAxisym2DLaw(const BorjaHenckyCamClayPlasticAxisym2DLaw& rOther)
      : NonLinearHenckyElasticPlasticAxisym2DLaw(rOther)
   {

   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveLaw::Pointer BorjaHenckyCamClayPlasticAxisym2DLaw::Clone() const
   {
      BorjaHenckyCamClayPlasticAxisym2DLaw::Pointer p_clone(new BorjaHenckyCamClayPlasticAxisym2DLaw(*this));
      return p_clone;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   BorjaHenckyCamClayPlasticAxisym2DLaw::~BorjaHenckyCamClayPlasticAxisym2DLaw()
   {
   }


   //*************************************** GET VALUE *********************************
   double&  BorjaHenckyCamClayPlasticAxisym2DLaw::GetValue(const Variable<double>& rThisVariable, double & rValue)
   {
      if ( rThisVariable == M_MODULUS ) {
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
         //G += Alpha*MeanStress;  // this modulus is approximated (not the real one)
         Matrix ELCG = mElasticLeftCauchyGreen;
         Vector Hencky;
         //Vector Hencky = ConvertCauchyGreenTensorToHenckyStrain( ELCG);
         // instead for calling the function, I just copy here the function. :'P --
         {
            Matrix & rCauchyGreenMatrix = ELCG;
            Matrix EigenVectors;
            Vector EigenValues;
            SolidMechanicsMathUtilities<double>::EigenVectors(rCauchyGreenMatrix, EigenVectors, EigenValues);

            Matrix Aux = ZeroMatrix(3,3);
            for (unsigned int i = 0; i < 3; ++i)
               Aux(i,i) = (std::log(EigenValues(i)))/2.0;

            Aux = prod(Aux, (EigenVectors));
            Aux = prod(trans(EigenVectors), Aux);
            Vector Result = MathUtils<double>::StrainTensorToVector(Aux, 6);
            Hencky = Result;
         }


         double volum = 0.0;
         for (unsigned int i = 0; i < 3; i++)
            volum += Hencky(i);

         double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
         double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
         ReferencePressure /= OCR;
         G += Alpha* ReferencePressure * std::exp( - volum / Swelling) ;

         if ( rThisVariable == SHEAR_MODULUS)
         {
            rValue = G;
            return rValue;
         }

         rValue = 9.0*K*G / ( 3.0*K + G);

      }
      else {
         rValue = NonLinearHenckyElasticPlasticAxisym2DLaw::GetValue( rThisVariable, rValue);
      }

      return rValue;
   }

   void BorjaHenckyCamClayPlasticAxisym2DLaw::SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
   {

         NonLinearHenckyElasticPlasticAxisym2DLaw::SetValue( rThisVariable, rValue, rCurrentProcessInfo );

   }


   int BorjaHenckyCamClayPlasticAxisym2DLaw::Check( const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
{

   return 0;
}


} // Namespace Kratos
