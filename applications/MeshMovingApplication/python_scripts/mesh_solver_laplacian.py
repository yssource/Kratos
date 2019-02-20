from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving

# Import baseclass
from KratosMultiphysics.MeshMovingApplication.mesh_solver_base import MeshSolverBase

import KratosMultiphysics.MeshMovingApplication.auxiliar_methods_solvers as aux_methods


def CreateSolver(mesh_model_part, custom_settings):
    return MeshSolverLaplacian(mesh_model_part, custom_settings)


class MeshSolverLaplacian(MeshSolverBase):
    def __init__(self, mesh_model_part, custom_settings):
        if custom_settings.Has("buffer_size"):
            buffer_size = custom_settings["buffer_size"].GetInt()
            if buffer_size < 2:
                raise Exception("A buffer_size of at least 2 is required!")
        else: # overwritting baseclass-default
            custom_settings.AddEmptyValue("buffer_size").SetInt(2)
        super(MeshSolverLaplacian, self).__init__(mesh_model_part, custom_settings)
        print("::[MeshSolverLaplacian]:: Construction finished")

    def Initialize(self):
        super(MeshSolverLaplacian, self).Initialize()

        # Doing the copy here because the ale-fluid-solver might create new mesh-solvers
        # in "Initialize" if it is doing mesh-motion on subdomains only
        self.laplacian_elements_part = aux_methods.CreateMeshMotionModelPart(
            self.mesh_model_part,
            "LaplacianMeshMovingElement"
        )

        self.mesh_motion_solving_strategies = []
        for i in range(self.settings["domain_size"].GetInt()):
            self.mesh_motion_solving_strategies.append(self._CreateSolvingStrategy())

        for solving_strategy in self.mesh_motion_solving_strategies:
            solving_strategy.Initialize()

        self.print_on_rank_zero("::[MeshSolverLaplacian]:: Finished initialization.")

    def SolveSolutionStep(self):
        KM.VariableUtils().UpdateCurrentToInitialConfiguration(model_part.Nodes) # ALL nodes!

        for i_dir, solving_strategy in enumerate(self.mesh_motion_solving_strategies):
            self.mesh_model_part.ProcessInfo[KM.LAPLACIAN_DIRECTION] = i_dir+1
            solving_strategy.Solve() # calling "Solve" is sufficient for this mesh-solver

        self.MoveMesh()

    def SetEchoLevel(self, level):
        for solving_strategy in self.mesh_motion_solving_strategies:
            solving_strategy.SetEchoLevel(level)

    def Clear(self):
        for solving_strategy in self.mesh_motion_solving_strategies:
            solving_strategy.Clear()

    def _CreateSolvingStrategy(self):
        linear_solver = self.get_linear_solver()
        time_order = self.settings["time_order"].GetInt()
        reform_dofs_each_step = self.settings["reform_dofs_each_step"].GetBool()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        calculate_mesh_velocities = self.settings["calculate_mesh_velocities"].GetBool()
        echo_level = self.settings["echo_level"].GetInt()
        solving_strategy = KratosMeshMoving.LaplacianMeshMovingStrategy(self.mesh_model_part,
                                                            linear_solver,
                                                            time_order,
                                                            reform_dofs_each_step,
                                                            compute_reactions,
                                                             calculate_mesh_velocities,
                                                            echo_level)

        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        builder_and_solver = self.create_builder_and_solver()
        return KM.ResidualBasedLinearStrategy(self.laplacian_elements_part,
                                                              mechanical_scheme,
                                                              linear_solver,
                                                              builder_and_solver,
                                                              self.settings["compute_reactions"].GetBool(),
                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                              False,
                                                              False)



