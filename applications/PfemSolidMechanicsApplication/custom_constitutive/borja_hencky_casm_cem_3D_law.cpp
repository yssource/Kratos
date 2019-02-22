//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                 October 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "custom_constitutive/borja_hencky_casm_cem_3D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

	//******************************CONSTRUCTOR*******************************************
  	//************************************************************************************

	BorjaHenckyCasmCemPlastic3DLaw::BorjaHenckyCasmCemPlastic3DLaw()
		: NonLinearHenckyElasticPlastic3DLaw()
  	{
		mpHardeningLaw   = HardeningLaw::Pointer( new CasmCemHardeningLaw() );
		mpYieldCriterion = YieldCriterion::Pointer( new CasmCemYieldCriterion(mpHardeningLaw) );
		mpFlowRule       = FlowRule::Pointer( new BorjaCasmCemExplicitFlowRule(mpYieldCriterion) );
std::cout<<"   CASM-CEM 3D constructed"<<std::endl;
	}

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

	BorjaHenckyCasmCemPlastic3DLaw::BorjaHenckyCasmCemPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
	{
		mpHardeningLaw    =  pHardeningLaw;
		mpYieldCriterion  =  YieldCriterion::Pointer( new CasmCemYieldCriterion(mpHardeningLaw) );
		mpFlowRule        =  pFlowRule;
std::cout<<"   CASM-CEM 3D constructed"<<std::endl;
	}

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  BorjaHenckyCasmCemPlastic3DLaw::BorjaHenckyCasmCemPlastic3DLaw(const BorjaHenckyCasmCemPlastic3DLaw& rOther)
		: NonLinearHenckyElasticPlastic3DLaw(rOther)
  {
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer BorjaHenckyCasmCemPlastic3DLaw::Clone() const
	{
		BorjaHenckyCasmCemPlastic3DLaw::Pointer p_clone(new BorjaHenckyCasmCemPlastic3DLaw(*this));
		return p_clone;
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  BorjaHenckyCasmCemPlastic3DLaw::~BorjaHenckyCasmCemPlastic3DLaw()
  {
  }


  //******************************** GET VALUE ********************************
	double& BorjaHenckyCasmCemPlastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue)
  {
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
		else if ( ( rThisVariable == YOUNG_MODULUS) || ( rThisVariable == SHEAR_MODULUS ) || (rThisVariable == BULK_MODULUS ) )
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
            //std::cout<<ThirdInvariantEffect<<std::endl<<std::endl;

        }
		else if ( (rThisVariable==VOLUMETRIC_PLASTIC) || (rThisVariable==INCR_SHEAR_PLASTIC) || (rThisVariable==PLASTIC_STRAIN) ||
		 (rThisVariable==INCR_VOL_PLASTIC) || (rThisVariable == PRECONSOLIDATION) || (rThisVariable == BONDING) )
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
			//bonding
			else {
				rValue = InternalVariables.Bonding;
			}
		}
		else
		{
			rValue = NonLinearHenckyElasticPlastic3DLaw::GetValue( rThisVariable, rValue);
		}

		return rValue;
  }
  
  void BorjaHenckyCasmCemPlastic3DLaw::SetPlasticVariables ( const double& rInitialPreconPressure, const double& rInitialBonding)
	{
		mpFlowRule->SetPlasticVariables(rInitialPreconPressure, rInitialBonding);
	}

	const double BorjaHenckyCasmCemPlastic3DLaw::GetBonding()
	{
		return mpFlowRule->GetPlasticVariables().Bonding;
	}
	
	const double BorjaHenckyCasmCemPlastic3DLaw::GetPreconPressure()
	{
		return mpFlowRule->GetPlasticVariables().PreconsolidationPressure;
	}
    const double BorjaHenckyCasmCemPlastic3DLaw::GetCriticalStateM()
    {
        double MM;
        MM = this->GetValue(CRITICAL_STATE_M, MM);
        return MM;
    }

  void BorjaHenckyCasmCemPlastic3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
  {
    if ( rThisVariable == PENALTY_PARAMETER) {
    }
    else if ( rThisVariable == PRECONSOLIDATION ) {
    	mpFlowRule->SetPreconsolidation(rValue);
    }
    else if ( rThisVariable == BONDING ) {
    	mpFlowRule->SetBonding(rValue);
    }
    else {
    	NonLinearHenckyElasticPlastic3DLaw::SetValue( rThisVariable, rValue, rCurrentProcessInfo );
    }
  }

	void BorjaHenckyCasmCemPlastic3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
  	{
			NonLinearHenckyElasticPlastic3DLaw::SetValue( rThisVariable, rValue, rCurrentProcessInfo );
	}

	int BorjaHenckyCasmCemPlastic3DLaw::Check( const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
	{
		return 0;
	}


} // Namespace Kratos
