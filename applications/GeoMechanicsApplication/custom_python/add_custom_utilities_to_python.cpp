//   
//   Project Name:        KratosGeoMechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/kratos_parameters.h"

#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/interface_element_utilities.hpp"


namespace Kratos
{
	
namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m) 
{
    using namespace pybind11;

#ifdef VG_NON_LOCAL    

    class_< FracturePropagation3DUtilities > 
    (m, "FracturePropagation3DUtilities")
    .def(init<>())
    .def("CheckFracturePropagation",&FracturePropagation3DUtilities::CheckFracturePropagation)
    .def("MappingModelParts",&FracturePropagation3DUtilities::MappingModelParts);
    
    class_< FracturePropagation2DUtilities >
    (m, "FracturePropagation2DUtilities")
    .def(init<>())
    .def("CheckFracturePropagation",&FracturePropagation2DUtilities::CheckFracturePropagation)
    .def("MappingModelParts",&FracturePropagation2DUtilities::MappingModelParts);

#endif

}

}  // namespace Python.
} // Namespace Kratos
