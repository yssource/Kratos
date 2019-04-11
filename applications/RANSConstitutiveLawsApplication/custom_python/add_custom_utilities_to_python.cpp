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
#include "custom_utilities/rans_variable_utils.h"

namespace Kratos
{
namespace Python
{
void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<RansVariableUtils>(m, "RansVariableUtils")
        .def(py::init<>())
        .def("AddToHistoricalNodeScalarVariable",
             &RansVariableUtils::AddToHistoricalNodeScalarVariable<Variable<double>>)
        .def("GetNumberOfNegativeScalarValueNodes", &RansVariableUtils::GetNumberOfNegativeScalarValueNodes)
        .def("GetMinimumScalarValue", &RansVariableUtils::GetMinimumScalarValue)
        .def("GetMaximumScalarValue", &RansVariableUtils::GetMaximumScalarValue)
        .def("GetScalarVariableIncreaseNormSquare", &RansVariableUtils::GetScalarVariableIncreaseNormSquare)
        .def("GetScalarVariableSolutionNormSquare",
             &RansVariableUtils::GetScalarVariableSolutionNormSquare)
        .def("FixScalarVariableDofs", &RansVariableUtils::FixScalarVariableDofs)
        ;

    py::class_<EvmKepsilonModelUtilities>(m, "EvmKepsilonModelUtilities")
        .def(py::init<>())
        .def("CalculateTurbulentViscosity", &EvmKepsilonModelUtilities::CalculateTurbulentViscosityForModelPart)
        .def("UpdateBoundaryConditions", &EvmKepsilonModelUtilities::UpdateBoundaryConditions)
        .def("AssignInitialValues", &EvmKepsilonModelUtilities::AssignInitialValues)
        ;
}

} // namespace Python.
} // Namespace Kratos
