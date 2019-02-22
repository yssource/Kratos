//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/borja_hencky_casm_axisym_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

namespace Kratos
{

	//******************************CONSTRUCTOR*******************************************
	//************************************************************************************

	BorjaHenckyCasmPlasticAxisym2DLaw::BorjaHenckyCasmPlasticAxisym2DLaw()
		: NonLinearHenckyElasticPlasticAxisym2DLaw()
	{
		mpHardeningLaw   = HardeningLaw::Pointer( new CasmHardeningLaw() );
		mpYieldCriterion = YieldCriterion::Pointer( new CasmYieldCriterion(mpHardeningLaw) );
		mpFlowRule       = FlowRule::Pointer( new BorjaCasmExplicitFlowRule(mpYieldCriterion) );
std::cout<<"   CASM 2D axisym constructed"<<std::endl;
	}

	//******************************CONSTRUCTOR*******************************************
	//************************************************************************************

	BorjaHenckyCasmPlasticAxisym2DLaw::BorjaHenckyCasmPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
	{
		mpHardeningLaw    =  pHardeningLaw;
		mpYieldCriterion  =  YieldCriterion::Pointer( new CasmYieldCriterion(mpHardeningLaw) );
		mpFlowRule        =  pFlowRule;
std::cout<<"   CASM 2D axisym constructed"<<std::endl;
	}

	//******************************COPY CONSTRUCTOR**************************************
	//************************************************************************************

	BorjaHenckyCasmPlasticAxisym2DLaw::BorjaHenckyCasmPlasticAxisym2DLaw(const BorjaHenckyCasmPlasticAxisym2DLaw& rOther)
		: NonLinearHenckyElasticPlasticAxisym2DLaw(rOther)
	{
	}

	//********************************CLONE***********************************************
	//************************************************************************************

	ConstitutiveLaw::Pointer BorjaHenckyCasmPlasticAxisym2DLaw::Clone() const
	{
		BorjaHenckyCasmPlasticAxisym2DLaw::Pointer p_clone(new BorjaHenckyCasmPlasticAxisym2DLaw(*this));
		return p_clone;
	}

	//*******************************DESTRUCTOR*******************************************
	//************************************************************************************

	BorjaHenckyCasmPlasticAxisym2DLaw::~BorjaHenckyCasmPlasticAxisym2DLaw()
	{
	}


	//*************************************** GET VALUE *********************************
	double&  BorjaHenckyCasmPlasticAxisym2DLaw::GetValue(const Variable<double>& rThisVariable, double & rValue)
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
		else if ( rThisVariable == CRITICAL_STATE_M )
        {
            const double InitialShearM = mpYieldCriterion->GetHardeningLaw().GetProperties()[ CRITICAL_STATE_LINE ];

            //stress invariants & invariants derivatives
            double LodeAngle;
            LodeAngle = this->GetValue( STRESS_INV_THETA, LodeAngle);
            LodeAngle *= 3.14159265359/180.0;
        
            //calculate third invariant effect 
            double ThirdInvariantEffect = 1.0;
            ThirdInvariantEffect = mpYieldCriterion->EvaluateThirdInvariantEffectMC( -LodeAngle);

            rValue = InitialShearM/ThirdInvariantEffect;
            std::cout<<ThirdInvariantEffect<<std::endl<<std::endl;

        }
		else if ( (rThisVariable==VOLUMETRIC_PLASTIC) || (rThisVariable==INCR_SHEAR_PLASTIC) || (rThisVariable==PLASTIC_STRAIN) ||
		 (rThisVariable==INCR_VOL_PLASTIC) || (rThisVariable == PRECONSOLIDATION) )
		{
			const PlasticVariablesType& InternalVariables = mpFlowRule->GetPlasticVariables();
			//volumetric plastic strain
			if ( rThisVariable==VOLUMETRIC_PLASTIC) {
				rValue = InternalVariables.EquivalentPlasticStrain; 
			}
			//plastic shear strain
			else if (rThisVariable==PLASTIC_STRAIN) {
				rValue = InternalVariables.PlasticShearStrain;
			}
			//incremental plastic shear strain
			else if (rThisVariable==INCR_SHEAR_PLASTIC) {
				rValue = InternalVariables.DeltaPlasticShearStrain;
			}
			//incremental volumetric plastic strain
			else if (rThisVariable==INCR_VOL_PLASTIC) {
				rValue = InternalVariables.DeltaEqPlasticStrain;
			}
			//preconsolidation pressure
			else if (rThisVariable==PRECONSOLIDATION) {
				rValue = InternalVariables.PreconsolidationPressure;
			}
			


			
		}

		else {
			rValue = NonLinearHenckyElasticPlasticAxisym2DLaw::GetValue( rThisVariable, rValue);
		}

		return rValue;
	}

	void BorjaHenckyCasmPlasticAxisym2DLaw::SetPlasticVariables ( const double& rInitialPreconPressure, const double& rInitialBonding) //TODO
	{
		mpFlowRule->SetPlasticVariables(rInitialPreconPressure, rInitialBonding);
	}

	const double BorjaHenckyCasmPlasticAxisym2DLaw::GetPreconPressure() //INHERT
	{
		return mpFlowRule->GetPlasticVariables().PreconsolidationPressure;
	}

	void BorjaHenckyCasmPlasticAxisym2DLaw::SetValue(const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) //INHERIT
	{
		if ( rThisVariable == PRECONSOLIDATION)
		{
			mpFlowRule->SetPreconsolidation(rValue);
		}
		else
		{
			NonLinearHenckyElasticPlasticAxisym2DLaw::SetValue( rThisVariable, rValue, rCurrentProcessInfo );
		}
	}





	void BorjaHenckyCasmPlasticAxisym2DLaw::SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
	{

			NonLinearHenckyElasticPlasticAxisym2DLaw::SetValue( rThisVariable, rValue, rCurrentProcessInfo );
	}
   
	int BorjaHenckyCasmPlasticAxisym2DLaw::Check( const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
	{
		return 0;
	}


} // Namespace Kratos
