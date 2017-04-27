/*
==============================================================================
KratosChimeraApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//
//   Project Name:        Kratos
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


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
	class_<CustomApplyMpcConstraint2dProcess,bases<Process> >("CustomApplyMpcConstraint2dProcess", init<>())
			.def("ApplyMpcConstraint2d", &CustomApplyMpcConstraint2dProcess::ApplyMpcConstraint2d);

	class_<CustomApplyMpcConstraint3dProcess,bases<Process> >("CustomApplyMpcConstraint3dProcess", init<>())
			.def("ApplyMpcConstraint3d", &CustomApplyMpcConstraint3dProcess::ApplyMpcConstraint3d);


			
}



} // namespace Python.

} // Namespace Kratos

