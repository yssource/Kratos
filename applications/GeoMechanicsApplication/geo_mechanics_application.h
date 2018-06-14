//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


#if !defined(KRATOS_GEO_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_GEO_MECHANICS_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

// Application includes
#include "geo_mechanics_application_variables.h"

#include "custom_conditions/U_Pw_force_condition.hpp"
#include "custom_conditions/U_Pw_face_load_condition.hpp"
#include "custom_conditions/U_Pw_normal_face_load_condition.hpp"
#include "custom_conditions/U_Pw_normal_flux_condition.hpp"
#include "custom_conditions/U_Pw_normal_flux_FIC_condition.hpp"
#include "custom_conditions/U_Pw_face_load_interface_condition.hpp"
#include "custom_conditions/U_Pw_normal_flux_interface_condition.hpp"
#include "custom_conditions/line_load_2D_diff_order_condition.hpp"
#include "custom_conditions/line_normal_load_2D_diff_order_condition.hpp"
#include "custom_conditions/line_normal_fluid_flux_2D_diff_order_condition.hpp"
#include "custom_conditions/surface_load_3D_diff_order_condition.hpp"
#include "custom_conditions/surface_normal_load_3D_diff_order_condition.hpp"
#include "custom_conditions/surface_normal_fluid_flux_3D_diff_order_condition.hpp"

#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_elements/U_Pw_small_strain_interface_element.hpp"
#include "custom_elements/U_Pw_small_strain_link_interface_element.hpp"
#include "custom_elements/U_Pw_small_strain_FIC_element.hpp"
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"

#ifdef VG_NON_LOCAL

#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"
#include "custom_constitutive/bilinear_cohesive_2D_law.hpp"

#include "custom_constitutive/custom_flow_rules/local_damage_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/nonlocal_damage_flow_rule.hpp"

#include "custom_constitutive/simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/simo_ju_local_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/simo_ju_local_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/simo_ju_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/simo_ju_nonlocal_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/modified_mises_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/modified_mises_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/modified_mises_nonlocal_damage_plane_stress_2D_law.hpp"

#endif


namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KratosGeoMechanicsApplication : public KratosApplication {
public:
	///@name Type Definitions
	///@{


	/// Pointer definition of KratosGeoMechanicsApplication
	KRATOS_CLASS_POINTER_DEFINITION(KratosGeoMechanicsApplication);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	KratosGeoMechanicsApplication();

	/// Destructor.
	virtual ~KratosGeoMechanicsApplication(){}


	///@}
	///@name Operators
	///@{


	///@}
	///@name Operations
	///@{

	virtual void Register() override;



	///@}
	///@name Access
	///@{


	///@}
	///@name Inquiry
	///@{


	///@}
	///@name Input and output
	///@{

	/// Turn back information as a string.
	virtual std::string Info() const override
	{
		return "KratosGeoMechanicsApplication";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const override
	{
		rOStream << Info();
		PrintData(rOStream);
	}

	///// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const override
	{
  		KRATOS_WATCH("in KratosGeoMechanicsApplication");
  		KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );

		rOStream << "Variables:" << std::endl;
		KratosComponents<VariableData>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Elements:" << std::endl;
		KratosComponents<Element>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Conditions:" << std::endl;
		KratosComponents<Condition>().PrintData(rOStream);
    }


	///@}
	///@name Friends
	///@{


	///@}

protected:
	///@name Protected static Member Variables
	///@{


	///@}
	///@name Protected member Variables
	///@{


	///@}
	///@name Protected Operators
	///@{


	///@}
	///@name Protected Operations
	///@{


	///@}
	///@name Protected  Access
	///@{


	///@}
	///@name Protected Inquiry
	///@{


	///@}
	///@name Protected LifeCycle
	///@{


	///@}

private:
	///@name Static Member Variables
	///@{

	// static const ApplicationCondition  msApplicationCondition;

	///@}
	///@name Member Variables
	///@{

	// const Elem2D   mElem2D;
	// const Elem3D   mElem3D;

	///@}
	///@name Private Operators
	///@{


	///@}
	///@name Private Operations
	///@{


	///@}
	///@name Private  Access
	///@{


	///@}
	///@name Private Inquiry
	///@{


	///@}
	///@name Un accessible methods
	///@{

	const UPwSmallStrainElement<2,3> mUPwSmallStrainElement2D3N;
	const UPwSmallStrainElement<2,4> mUPwSmallStrainElement2D4N;
	const UPwSmallStrainElement<3,4> mUPwSmallStrainElement3D4N;
	const UPwSmallStrainElement<3,8> mUPwSmallStrainElement3D8N;

	const UPwSmallStrainInterfaceElement<2,4> mUPwSmallStrainInterfaceElement2D4N;
	const UPwSmallStrainInterfaceElement<3,6> mUPwSmallStrainInterfaceElement3D6N;
	const UPwSmallStrainInterfaceElement<3,8> mUPwSmallStrainInterfaceElement3D8N;

	const UPwSmallStrainLinkInterfaceElement<2,4> mUPwSmallStrainLinkInterfaceElement2D4N;
	const UPwSmallStrainLinkInterfaceElement<3,6> mUPwSmallStrainLinkInterfaceElement3D6N;
	const UPwSmallStrainLinkInterfaceElement<3,8> mUPwSmallStrainLinkInterfaceElement3D8N;

	const UPwSmallStrainFICElement<2,3> mUPwSmallStrainFICElement2D3N;
	const UPwSmallStrainFICElement<2,4> mUPwSmallStrainFICElement2D4N;
	const UPwSmallStrainFICElement<3,4> mUPwSmallStrainFICElement3D4N;
	const UPwSmallStrainFICElement<3,8> mUPwSmallStrainFICElement3D8N;

	const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D6N;
	const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D8N;
	const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D9N;
	const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement3D10N;
	const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement3D20N;
	const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement3D27N;


	const UPwForceCondition<2,1> mUPwForceCondition2D1N;
	const UPwForceCondition<3,1> mUPwForceCondition3D1N;
	const UPwFaceLoadCondition<2,2> mUPwFaceLoadCondition2D2N;
	const UPwFaceLoadCondition<3,3> mUPwFaceLoadCondition3D3N;
	const UPwFaceLoadCondition<3,4> mUPwFaceLoadCondition3D4N;
	const UPwNormalFaceLoadCondition<2,2> mUPwNormalFaceLoadCondition2D2N;
	const UPwNormalFaceLoadCondition<3,3> mUPwNormalFaceLoadCondition3D3N;
	const UPwNormalFaceLoadCondition<3,4> mUPwNormalFaceLoadCondition3D4N;
	const UPwNormalFluxCondition<2,2> mUPwNormalFluxCondition2D2N;
	const UPwNormalFluxCondition<3,3> mUPwNormalFluxCondition3D3N;
	const UPwNormalFluxCondition<3,4> mUPwNormalFluxCondition3D4N;

	const UPwFaceLoadInterfaceCondition<2,2> mUPwFaceLoadInterfaceCondition2D2N;
	const UPwFaceLoadInterfaceCondition<3,4> mUPwFaceLoadInterfaceCondition3D4N;
	const UPwNormalFluxInterfaceCondition<2,2> mUPwNormalFluxInterfaceCondition2D2N;
	const UPwNormalFluxInterfaceCondition<3,4> mUPwNormalFluxInterfaceCondition3D4N;

	const UPwNormalFluxFICCondition<2,2> mUPwNormalFluxFICCondition2D2N;
	const UPwNormalFluxFICCondition<3,3> mUPwNormalFluxFICCondition3D3N;
	const UPwNormalFluxFICCondition<3,4> mUPwNormalFluxFICCondition3D4N;

	const LineLoad2DDiffOrderCondition mLineLoadDiffOrderCondition2D3N;
	const LineNormalLoad2DDiffOrderCondition mLineNormalLoadDiffOrderCondition2D3N;
	const LineNormalFluidFlux2DDiffOrderCondition mLineNormalFluidFluxDiffOrderCondition2D3N;
	const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D6N;
	const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D8N;
	const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D9N;
	const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D6N;
	const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D8N;
	const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D9N;
	const SurfaceNormalFluidFlux3DDiffOrderCondition mSurfaceNormalFluidFluxDiffOrderCondition3D6N;
    const SurfaceNormalFluidFlux3DDiffOrderCondition mSurfaceNormalFluidFluxDiffOrderCondition3D8N;
    const SurfaceNormalFluidFlux3DDiffOrderCondition mSurfaceNormalFluidFluxDiffOrderCondition3D9N;




	/// Assignment operator.
	KratosGeoMechanicsApplication& operator=(KratosGeoMechanicsApplication const& rOther);

	/// Copy constructor.
	KratosGeoMechanicsApplication(KratosGeoMechanicsApplication const& rOther);


	///@}

}; // Class KratosGeoMechanicsApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_GEO_MECHANICS_APPLICATION_H_INCLUDED  defined
