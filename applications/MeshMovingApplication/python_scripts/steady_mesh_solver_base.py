from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MeshMovingApplication as KMM

# Other imports
from KratosMultiphysics.MeshMovingApplication.mesh_solver_base import MeshSolverBase


class SteadyMeshSolverBase(MeshSolverBase):
    """The base class for mesh motion solvers.

    This class defines the user interface to mesh motion solvers.

    Derived classes must override the function _create_mesh_motion_solving_strategy()
    to customize the mesh motion algorithm. The mesh motion solving strategy and linear
    solver should always be retrieved using the getter functions. Only the
    member variables listed below should be accessed directly.

    Public member variables:
    settings -- Kratos parameters for general mesh motion settings.
    mesh_model_part -- the mesh motion model part.
    """
    def __init__(self, model, custom_settings):
        super(MeshSolverBase,self).__init__(model, custom_settings)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type"           : "steady_mesh_solver_structural_similarity",
            "buffer_size"           : 1,
            "echo_level"            : 0,
            "domain_size"           : -1,
            "model_part_name"       : "",
            "model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name"
            },
            "mesh_motion_linear_solver_settings" : {
                "solver_type" : "amgcl",
                "smoother_type":"ilu0",
                "krylov_type": "gmres",
                "coarsening_type": "aggregation",
                "max_iteration": 200,
                "provide_coordinates": false,
                "gmres_krylov_space_dimension": 100,
                "verbosity" : 0,
                "tolerance": 1e-7,
                "scaling": false,
                "block_size": 1,
                "use_block_matrices_if_possible" : true,
                "coarse_enough" : 5000
            },
            "reform_dofs_each_step"     : false,
            "compute_reactions"         : false
        }""")

        self.settings.ValidateAndAssignDefaults(default_settings)

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        if self.model.HasModelPart(model_part_name):
            self.mesh_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.mesh_model_part = model.CreateModelPart(model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size == -1:
            raise Exception('Please provide the domain size as the "domain_size" (int) parameter!')

        self.mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        KratosMultiphysics.Logger.PrintInfo("::[SteadyMeshSolverBase]:: Construction finished")

    #### Public user interface functions ####

    def AddVariables(self):
        self.mesh_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)
        self.mesh_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_REACTION)
        KratosMultiphysics.Logger.PrintInfo("::[SteadyMeshSolverBase]:: Variables ADDED.")

    def AdvanceInTime(self, current_time):
        self.mesh_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        return 0.0

    def Predict(self):
        self.get_mesh_motion_solving_strategy().Predict() #needed?

    def MoveMesh(self):
        # move local and ghost nodes
        KMM.MoveMesh(self.mesh_model_part.Nodes)

    def PrepareModelPart(self):
        self.mesh_model_part.SetBufferSize(1)
