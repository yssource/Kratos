//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "add_imposed_deformation_to_python.h"
#include "includes/define_python.h"
#include "includes/imposed_deformation.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

typedef ImposedDeformation ImposedDeformationBaseType;
template<class TVariableType> bool ImposedDeformationHas(ImposedDeformation& rThisImposedDeformation, TVariableType const& rThisVariable) { return rThisImposedDeformation.Has(rThisVariable); }

//dirty trick. give back a copy instead of a reference
template<class TDataType> const TDataType ImposedDeformationGetValue(ImposedDeformation& rThisImposedDeformation, const ConstitutiveLaw* pConstitutiveLaw, const Variable<TDataType >& rThisVariable, TDataType& value )
{
    TDataType tmp = rThisImposedDeformation.GetValue(pConstitutiveLaw, rThisVariable, value);
    return tmp;
}

// Function to export CalculateValue(...).
// Returns a copy instead of a reference, as GetValue does.
template<class TDataType> const TDataType ImposedDeformationCalculateValue(
    ImposedDeformation& rThisImposedDeformation,
    const ConstitutiveLaw* pConstitutiveLaw,
    ConstitutiveLaw::Parameters& rValues,
    const Variable<TDataType >& rThisVariable,
    TDataType& rValue)
{
    rThisImposedDeformation.CalculateValue(pConstitutiveLaw, rValues, rThisVariable, rValue);
    return rValue;
}

template<class TDataType> void ImposedDeformationSetValue(ImposedDeformation& rThisImposedDeformation, const ConstitutiveLaw* pConstitutiveLaw, const Variable<TDataType>& rThisVariable, TDataType& value, const ProcessInfo& rCurrentProcessInfo)
{ rThisImposedDeformation.SetValue(pConstitutiveLaw, rThisVariable, value, rCurrentProcessInfo); }

void NewInterfaceCalculateResponse(ImposedDeformation& rThisImposedDeformation, const ConstitutiveLaw* pConstitutiveLaw, ConstitutiveLaw::Parameters& rValues,const ConstitutiveLaw::StressMeasure& rStressMeasure)
{rThisImposedDeformation.CalculateResponse (pConstitutiveLaw, rValues,rStressMeasure);}

void AddImposedDeformationToPython(pybind11::module& m)
{
    py::class_< ImposedDeformation, ImposedDeformation::Pointer, Flags>(m,"ImposedDeformation")
    .def(py::init<>() )
    .def("Create",&ImposedDeformation::Create)
    .def("Clone",&ImposedDeformation::Clone)
    .def("Has", &ImposedDeformationHas< Variable<bool> >)
    .def("Has", &ImposedDeformationHas< Variable<int> >)
    .def("Has", &ImposedDeformationHas< Variable<double> >)
    .def("Has", &ImposedDeformationHas< Variable<array_1d<double,3> > >)
    .def("Has", &ImposedDeformationHas< Variable<Vector> >)
    .def("Has", &ImposedDeformationHas< Variable<Matrix> >)
    .def("GetValue", &ImposedDeformationGetValue<bool> )
    .def("GetValue", &ImposedDeformationGetValue<int> )
    .def("GetValue", &ImposedDeformationGetValue<double> )
    .def("GetValue", &ImposedDeformationGetValue<array_1d<double,3>  >)
    .def("GetValue", &ImposedDeformationGetValue<Vector >)
    .def("GetValue", &ImposedDeformationGetValue<Matrix >)
    .def("SetValue", &ImposedDeformationSetValue<int> )
    .def("SetValue", &ImposedDeformationSetValue<double> )
    .def("SetValue", &ImposedDeformationSetValue<array_1d<double,3>  >)
    .def("SetValue", &ImposedDeformationSetValue<Vector >)
    .def("SetValue", &ImposedDeformationSetValue<Matrix >)
    .def("CalculateValue", &ImposedDeformationCalculateValue<bool> )
    .def("CalculateValue", &ImposedDeformationCalculateValue<int> )
    .def("CalculateValue", &ImposedDeformationCalculateValue<double> )
    .def("CalculateValue", &ImposedDeformationCalculateValue<array_1d<double,3>  >)
    .def("CalculateValue", &ImposedDeformationCalculateValue<Vector >)
    .def("CalculateValue", &ImposedDeformationCalculateValue<Matrix >)
    .def("CalculateResponse",&NewInterfaceCalculateResponse)
    .def("CalculateResponsePK1",&ImposedDeformation::CalculateResponsePK1)
    .def("CalculateResponsePK2",&ImposedDeformation::CalculateResponsePK2)
    .def("CalculateResponseKirchhoff",&ImposedDeformation::CalculateResponseKirchhoff)
    .def("CalculateResponseCauchy",&ImposedDeformation::CalculateResponseCauchy)
    .def("InitializeResponse",&ImposedDeformation::InitializeResponse)
    .def("InitializeResponsePK1",&ImposedDeformation::InitializeResponsePK1)
    .def("InitializeResponsePK2",&ImposedDeformation::InitializeResponsePK2)
    .def("InitializeResponseKirchhoff",&ImposedDeformation::InitializeResponseKirchhoff)
    .def("InitializeResponseCauchy",&ImposedDeformation::InitializeResponseCauchy)
    .def("FinalizeResponse",&ImposedDeformation::FinalizeResponse)
    .def("FinalizeResponsePK1",&ImposedDeformation::FinalizeResponsePK1)
    .def("FinalizeResponsePK2",&ImposedDeformation::FinalizeResponsePK2)
    .def("FinalizeResponseKirchhoff",&ImposedDeformation::FinalizeResponseKirchhoff)
    .def("FinalizeResponseCauchy",&ImposedDeformation::FinalizeResponseCauchy)
    .def("Initialize",&ImposedDeformation::Initialize)
    .def("Reset",&ImposedDeformation::Reset)
    .def("Check",&ImposedDeformation::Check)
    ;

}
}  // namespace Python.
}  // namespace Kratos.
