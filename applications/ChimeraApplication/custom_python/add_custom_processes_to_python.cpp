//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//                   Navaneeth K Narayanan
//
// ==============================================================================

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/custom_hole_cutting_process.h"
#include "custom_processes/custom_extract_variables_process.h"
#include "custom_processes/custom_apply_chimera_using_mpc_process.h"
#include "custom_processes/custom_calculate_signed_distance_process.h"
//#include "custom_processes/calculate_signed_distance_to_2d_skin_process.h"
#include "custom_processes/calculate_signed_distance_to_2d_condition_skin_process.h"
#include "custom_processes/calculate_chimera_signed_distance_to_3d_condition_skin_process.h"
#include "custom_processes/apply_multi_point_constraints_process.h"
namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython()
{
	using namespace boost::python;

	/*
	 *  Custom holeCuttingProcess
	 */
	class_<CustomHoleCuttingProcess,bases<Process> >("CustomHoleCuttingProcess", init<>())
		.def("ExtractMeshAtCentroidDistance", &CustomHoleCuttingProcess::ExtractMeshAtCentroidDistance)
		.def("ExtractMeshBetweenLimits", &CustomHoleCuttingProcess::ExtractMeshBetweenLimits);
		

	/*
	 * CustomCalculateSignedDistanceProcess
	 */ 
	class_<CustomCalculateSignedDistanceProcess<2>,bases<Process> >("CustomCalculateSignedDistanceProcess2D", init<>())
			.def("ExtractDistance", &CustomCalculateSignedDistanceProcess<2>::ExtractDistance);

	class_<CustomCalculateSignedDistanceProcess<3>,bases<Process> >("CustomCalculateSignedDistanceProcess3D", init<>())
			.def("ExtractDistance", &CustomCalculateSignedDistanceProcess<3>::ExtractDistance);		

	/*
	 * CustomExtractVariablesProcess
	 */
	class_<CustomExtractVariablesProcess,bases<Process> >("CustomExtractVariablesProcess", init<>())
			.def("ExtractVariable", &CustomExtractVariablesProcess::ExtractVariable< array_1d<double, 3> >)
			.def("ExtractVariable", &CustomExtractVariablesProcess::ExtractVariable<double>);			


    
	/*
	 * CustomApplyMpcConstraintProcessforChimera for 2d and 3d
	 */
	class_<CustomApplyChimeraUsingMpcProcess<2>,bases<Process> >("CustomApplyChimeraUsingMpcProcess2d", init<ModelPart&,ModelPart&,ModelPart&,double>())
			.def("ApplyMpcConstraint", &CustomApplyChimeraUsingMpcProcess<2>::ApplyMpcConstraint)
			.def("ApplyChimeraUsingMpc2d", &CustomApplyChimeraUsingMpcProcess<2>::ApplyChimeraUsingMpc)
			.def("SetOverlapDistance",&CustomApplyChimeraUsingMpcProcess<2>::SetOverlapDistance)
			.def("CalculateNodalAreaAndNodalMass",&CustomApplyChimeraUsingMpcProcess<2>::CalculateNodalAreaAndNodalMass);
			

	class_<CustomApplyChimeraUsingMpcProcess<3>,bases<Process> >("CustomApplyChimeraUsingMpcProcess3d", init<ModelPart&,ModelPart&,ModelPart&,double>())
			.def("ApplyMpcConstraint", &CustomApplyChimeraUsingMpcProcess<3>::ApplyMpcConstraint)		
			.def("ApplyChimeraUsingMpc3d", &CustomApplyChimeraUsingMpcProcess<3>::ApplyChimeraUsingMpc)
			.def("SetOverlapDistance",&CustomApplyChimeraUsingMpcProcess<3>::SetOverlapDistance)
			.def("CalculateNodalAreaAndNodalMass",&CustomApplyChimeraUsingMpcProcess<3>::CalculateNodalAreaAndNodalMass);
			

	/*
	 * Calculate_signed_distance_2d
	 */

	//class_<CalculateSignedDistanceTo2DSkinProcess,bases<Process> >("CalculateSignedDistanceTo2DSkinProcess", init<ModelPart&,ModelPart&>())
		//.def("Execute", &CalculateSignedDistanceTo2DSkinProcess::Execute);	
	class_<CalculateChimeraSignedDistanceTo3DConditionSkinProcess,bases<Process> >("CalculateChimeraSignedDistanceTo3DConditionSkinProcess", init<ModelPart&,ModelPart&>())
		.def("Execute", &CalculateChimeraSignedDistanceTo3DConditionSkinProcess::Execute);	
			
}



} // namespace Python.

} // Namespace Kratos

