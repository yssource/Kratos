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
#include "processes/process.h"
#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

//Processes
#include "custom_processes/apply_multi_point_constraints_process_chimera.h"
#include "custom_utilities/vtk_output.hpp"
#include "custom_utilities/quadtree_binary_cell.h"
#include "custom_utilities/quadtree_binary.h"
//#include "custom_utilities/multipoint_constraint_data.hpp"

namespace Kratos
{

namespace Python
{


  void  AddCustomUtilitiesToPython()
  {
    using namespace boost::python;


      //typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
      //typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
      //typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


          /// Processes
      class_<ApplyMultipointConstraintsProcessChimera, boost::noncopyable, bases<Process>>("ApplyMultipointConstraintsProcessChimera", init<ModelPart&>())
      .def(init< ModelPart&, Parameters& >())
      .def("AddMasterSlaveRelation", &ApplyMultipointConstraintsProcessChimera::AddMasterSlaveRelation)
      .def("PrintData", &ApplyMultipointConstraintsProcessChimera::PrintData);

      class_<VtkOutput, boost::noncopyable>("VtkOutput", init< ModelPart&, std::string, Parameters >())
      .def("PrintOutput", &VtkOutput::PrintOutput)
      .def("PrintOutput", &VtkOutput::PrintOutputSubModelPart);      


  }





}  // namespace Python.

} // Namespace Kratos
