// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "custom_python/add_custom_schemes_to_python.h"
#include "custom_schemes/adjoint_bossak_scheme.h"
#include "custom_schemes/convergence_criteria/adjoint_vel_pr_criteria.h"

namespace Kratos
{
namespace Python
{
using namespace pybind11;

void AddCustomSchemesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;

    class_<
        AdjointBossakScheme<SparseSpaceType, LocalSpaceType>,
        typename AdjointBossakScheme<SparseSpaceType, LocalSpaceType>::Pointer,
        SchemeType>(m,"AdjointBossakScheme")
        .def(init<Parameters&, ResponseFunction::Pointer>())
        ;

	// Convergence criteria
    class_< AdjointVelPrCriteria< SparseSpaceType, LocalSpaceType >,
            typename AdjointVelPrCriteria< SparseSpaceType, LocalSpaceType >::Pointer,
            ConvergenceCriteria< SparseSpaceType, LocalSpaceType  >>
            (m,"AdjointVelPrCriteria")
            .def(init< double, double>())
            .def("SetEchoLevel",&AdjointVelPrCriteria<SparseSpaceType, LocalSpaceType >::SetEchoLevel)
            ;        
}

} // namespace Python

} // namespace Kratos
