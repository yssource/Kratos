//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/define.h"

// rans modelling includes
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utils.h"

namespace Kratos
{
namespace Python
{
void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<RansVariableUtils, VariableUtils>(m, "RansVariableUtils")
        .def(py::init<>())
        .def("FixScalarVariableDofs", &RansVariableUtils::FixScalarVariableDofs)
        .def("ClipScalarVariable", &RansVariableUtils::ClipScalarVariable)
        .def("GetNumberOfNegativeScalarValueNodes", &RansVariableUtils::GetNumberOfNegativeScalarValueNodes)
        .def("GetMinimumScalarValue", &RansVariableUtils::GetMinimumScalarValue)
        .def("GetMaximumScalarValue", &RansVariableUtils::GetMaximumScalarValue)
        .def("GetScalarVariableDifferenceNormSquare",
             &RansVariableUtils::GetScalarVariableDifferenceNormSquare)
        .def("GetScalarVariableSolutionNormSquare", &RansVariableUtils::GetScalarVariableSolutionNormSquare)
        .def("CopyNodalSolutionStepVariablesList",
             &RansVariableUtils::CopyNodalSolutionStepVariablesList);

    py::class_<EvmKepsilonModelUtilities>(m, "EvmKepsilonModelUtilities")
        .def(py::init<>())
        .def("UpdateBoundaryConditions", &EvmKepsilonModelUtilities::UpdateBoundaryConditions)
        .def("AssignInitialValues", &EvmKepsilonModelUtilities::AssignInitialValues);

    py::class_<RansCalculationUtilities>(m, "RansCalculationUtilities")
        .def(py::init<>())
        .def("CreateConnectivityPreservingModelPart",
             &RansCalculationUtilities::CreateConnectivityPreservingModelPart);
}

} // namespace Python.
} // Namespace Kratos
