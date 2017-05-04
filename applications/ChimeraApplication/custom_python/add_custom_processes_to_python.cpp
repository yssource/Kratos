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
#include "custom_processes/custom_whole_cutting_process.h"
#include "custom_processes/custom_extract_variables_process.h"
#include "custom_processes/custom_apply_mpc_constraint_2d_process.h"
#include "custom_processes/custom_apply_mpc_constraint_3d_process.h"

namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython()
{
	using namespace boost::python;

	/*
	 *  Custom WholeCuttingProcess
	 */
	class_<CustomWholeCuttingProcess,bases<Process> >("CustomWholeCuttingProcess", init<>())
		.def("ExtractSurfaceMeshAtDistance", &CustomWholeCuttingProcess::ExtractSurfaceMeshAtDistance)
		.def("ExtractVolumeMeshBetweenLimits", &CustomWholeCuttingProcess::ExtractVolumeMeshBetweenLimits);

	/*
	 * CustomExtractVariablesProcess
	 */
	class_<CustomExtractVariablesProcess,bases<Process> >("CustomExtractVariablesProcess", init<>())
			.def("ExtractVariable", &CustomExtractVariablesProcess::ExtractVariable< array_1d<double, 3> >)
			.def("ExtractVariable", &CustomExtractVariablesProcess::ExtractVariable<double>);

	/*
	 * CustomApplyMpcConstraintProcessforChimera
	 */
	class_<CustomApplyMpcConstraint2dProcess,bases<Process> >("CustomApplyMpcConstraint2dProcess", init<ModelPart&>())
			.def("ApplyMpcConstraint2d", &CustomApplyMpcConstraint2dProcess::ApplyMpcConstraint2d);

	class_<CustomApplyMpcConstraint3dProcess,bases<Process> >("CustomApplyMpcConstraint3dProcess", init<>())
			.def("ApplyMpcConstraint3d", &CustomApplyMpcConstraint3dProcess::ApplyMpcConstraint3d);


			
}



} // namespace Python.

} // Namespace Kratos

