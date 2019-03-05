//
//   Project Name:        KratosPlasmaDynamicsApplication $
//   Last Modified by:    $Author:    Marc Chung To Sang $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;
}

}  // namespace Python.
} // Namespace Kratos
