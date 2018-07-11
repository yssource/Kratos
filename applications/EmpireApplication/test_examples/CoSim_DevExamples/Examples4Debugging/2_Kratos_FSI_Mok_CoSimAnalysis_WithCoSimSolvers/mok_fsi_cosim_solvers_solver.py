from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructuralMechanics

import fsi_utilities # here auxiliary functions e.g. for relaxation are declared

# Importing the base class
from co_simulation_base_solver import CoSimulationBaseSolver

# Other imports
import co_simulation_tools as cosim_tools
from co_simulation_tools import csprint, red, green, cyan, bold, magenta

def CreateSolver(cosim_solver_settings, level):
    return MokFSISolverWithCoSimSolvers(cosim_solver_settings, level)

class MokFSISolverWithCoSimSolvers(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings, level):
        super(MokFSISolverWithCoSimSolvers, self).__init__(cosim_solver_settings, level)

        if not len(self.cosim_solver_settings["solvers"]) == 2:
            raise Exception("Exactly two solvers have to be specified for the GaussSeidelStrongCouplingSolver!")

        self.solver_names = []
        self.solvers = {}

        import python_solvers_wrapper_co_simulation as solvers_wrapper

        for solver_name in self.cosim_solver_settings["solvers"]:
            self.solver_names.append(solver_name)
            self.cosim_solver_settings["solvers"][solver_name]["name"] = solver_name # adding the name such that the solver can identify itself
            self.solvers[solver_name] = solvers_wrapper.CreateSolver(
                self.cosim_solver_settings["solvers"][solver_name], self.lvl-1) # -1 to have solver prints on same lvl

    def Initialize(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].Initialize()

        self.fluid_model_part = self.solvers["fluid"].model["MainModelPart"]
        self.structural_model_part = self.solvers["structure"].model["Structure"]

        # FSI parameters
        self.max_iter = 10    # number of inner iterations (set to 1 for explicit coupling)
        self.interface_epsilon = 1e-5  # interface residual (only needed for implicit coupling)
        self.relaxation_coefficient = 0.125 # initial value

        # -------------------
        # ----- Mapping -----
        # -------------------
        with open("ProjectParametersCFD.json", 'r') as parameter_file:
            mapper_params = KratosMultiphysics.Parameters(parameter_file.read())

        project_parameters_mapper_1 = mapper_params["mapper_settings"][0]
        project_parameters_mapper_2 = mapper_params["mapper_settings"][1]
        project_parameters_mapper_3 = mapper_params["mapper_settings"][2]

        self.mapper_1 = KratosMapping.MapperFactory.CreateMapper(self.structural_model_part, self.fluid_model_part, project_parameters_mapper_1)
        self.mapper_2 = KratosMapping.MapperFactory.CreateMapper(self.structural_model_part, self.fluid_model_part, project_parameters_mapper_2)
        self.mapper_3 = KratosMapping.MapperFactory.CreateMapper(self.structural_model_part, self.fluid_model_part, project_parameters_mapper_3)

    def Finalize(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].Finalize()

    def AdvanceInTime(self, current_time):
        self.time = self.solvers[self.solver_names[0]].AdvanceInTime(current_time)
        for solver_name in self.solver_names[1:]:
            time_other_solver = self.solvers[solver_name].AdvanceInTime(current_time)
            if abs(self.time - time_other_solver) > 1e-12:
                raise Exception("Solver time mismatch") # SubStepping not implemented!

        return self.time

    def Predict(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].Predict()

    def InitializeSolutionStep(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].InitializeSolutionStep()

        self.residual = 1
        self.old_displacements = fsi_utilities.GetDisplacements(self.structural_model_part.GetSubModelPart("GENERIC_Beam").Nodes, 2)

    def FinalizeSolutionStep(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].FinalizeSolutionStep()

    def OutputSolutionStep(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].OutputSolutionStep()

    def SolveSolutionStep(self):
        for k in range(self.max_iter):
            # Apply Dirichlet B.C.'s from structural solver to mesh solver
            DisplacementToMesh(self.mapper_1)
            DisplacementToMesh(self.mapper_2)
            DisplacementToMesh(self.mapper_3)

            # Mesh and Fluid are currently solved independently, since the ALE solver does not copy the mesh velocity
            # Solve Mesh
            self.solvers["fluid"].SolveSolutionStep()

            # Apply Neumann B.C.'s from fluid solver to structural solver
            NeumannToStructure(self.mapper_1, KratosMapping.Mapper.SWAP_SIGN | KratosMapping.Mapper.CONSERVATIVE)
            NeumannToStructure(self.mapper_2, KratosMapping.Mapper.SWAP_SIGN | KratosMapping.Mapper.ADD_VALUES | KratosMapping.Mapper.CONSERVATIVE)

            # # Solver Structure
            self.solvers["structure"].SolveSolutionStep()

            # Convergence Checking (only for implicit coupling)
            if self.max_iter > 1:
                displacements = fsi_utilities.GetDisplacements(self.structural_model_part.GetSubModelPart("GENERIC_Beam").Nodes, 2)

                # Compute Residual
                old_residual = self.residual
                self.residual = fsi_utilities.CalculateResidual(displacements,self.old_displacements)

                if (fsi_utilities.Norm(self.residual) <= self.interface_epsilon):
                    fsi_utilities.SetDisplacements(displacements,self.structural_model_part.GetSubModelPart("GENERIC_Beam").Nodes, 2)
                    print("******************************************************")
                    print("************ CONVERGENCE AT INTERFACE ACHIEVED *******")
                    print("******************************************************")
                    break # TODO check if this works bcs it is nested
                else:
                    self.relaxation_coefficient = fsi_utilities.ComputeAitkenRelaxation(self.relaxation_coefficient, self.residual, old_residual, k)
                    relaxed_displacements = fsi_utilities.CalculateRelaxation(self.relaxation_coefficient, self.old_displacements, self.residual)
                    self.old_displacements = relaxed_displacements
                    fsi_utilities.SetDisplacements(relaxed_displacements, self.structural_model_part.GetSubModelPart("GENERIC_Beam").Nodes, 2)

                if (k+1 >= self.max_iter):
                    print("######################################################")
                    print("##### CONVERGENCE AT INTERFACE WAS NOT ACHIEVED ######")
                    print("######################################################")

                print("==========================================================")
                print("COUPLING RESIDUAL = ", fsi_utilities.Norm(self.residual))
                print("COUPLING ITERATION = ", k+1, "/", self.max_iter)
                print("RELAXATION COEFFICIENT = ", self.relaxation_coefficient)
                print("==========================================================")


        # for k in range(self.num_coupling_iterations):
        #     csprint(self.lvl, cyan("Coupling iteration: ")+bold(str(k+1)+" / " + str(self.num_coupling_iterations)))
        #     for solver_name in self.solver_names:
        #         solver = self.solvers[solver_name]
        #         self.__SynchronizeInputData(solver, solver_name)
        #         solver.SolveSolutionStep()
        #         self.__SynchronizeOutputData(solver, solver_name)

        #     if self.convergence_criteria.IsConverged():
        #         csprint(self.lvl, green("##### CONVERGENCE AT INTERFACE WAS ACHIEVED #####"))
        #         break
        #     else:
        #         self.convergence_accelerator.ComputeUpdate()

        #     if k+1 >= self.num_coupling_iterations:
        #         csprint(self.lvl, red("XXXXX CONVERGENCE AT INTERFACE WAS NOT ACHIEVED XXXXX"))


def NeumannToStructure(mapper, flag):
    mapper.InverseMap(KratosStructuralMechanics.POINT_LOAD, KratosMultiphysics.REACTION, flag)

def DisplacementToMesh(mapper):
    mapper.Map(KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.MESH_DISPLACEMENT)