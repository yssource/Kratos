# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from .algorithm_base import OptimizationAlgorithm
from . import mapper_factory
from . import data_logger_factory
from .custom_timer import Timer
from .custom_variable_utilities import WriteDictionaryDataOnNodalVariable

import numpy as np

# ==============================================================================
class AlgorithmSteepestDescentWithProjection(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = Parameters("""
        {
            "name"                  : "steepest_descent_with_projection",
            "max_iterations"        : 100,
            "relative_tolerance"    : 1e-3,
            "gradient_tolerance"    : 1e-5,
            "update_mapping_matrix" : true,
            "fix_boundaries"        : [],
            "line_search" : {
                "line_search_type"           : "manual_stepping",
                "normalize_search_direction" : true,
                "step_size"                  : 1.0,
                "approximation_tolerance"    : 0.1,
                "increase_factor"            : 1.1,
                "max_increase_factor"        : 10.0
            }
        }""")
        self.algorithm_settings =  optimization_settings["optimization_algorithm"]
        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

        self.optimization_settings = optimization_settings
        self.mapper_settings = optimization_settings["design_variables"]["filter"]

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.design_surface = None
        self.mapper = None
        self.data_logger = None
        self.optimization_utilities = None

        self.objectives = optimization_settings["objectives"]
        self.constraints = optimization_settings["constraints"]

        self.previos_objective_value = None

        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()
        self.gradient_tolerance = self.algorithm_settings["gradient_tolerance"].GetDouble()
        self.update_mapping_matrix = self.algorithm_settings["update_mapping_matrix"].GetBool()
        self.line_search_type = self.algorithm_settings["line_search"]["line_search_type"].GetString()
        self.normalize_search_direction = self.algorithm_settings["line_search"]["normalize_search_direction"].GetBool()
        self.approximation_tolerance = self.algorithm_settings["line_search"]["approximation_tolerance"].GetDouble()
        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.increase_factor = self.algorithm_settings["line_search"]["increase_factor"].GetDouble()
        self.max_step_size = self.step_size*self.algorithm_settings["line_search"]["max_increase_factor"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.SEARCH_DIRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.VECTOR_VARIABLE)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.VECTOR_VARIABLE_MAPPED)
        self.optimization_model_part.AddNodalSolutionStepVariable(KM.NODAL_VAUX)

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Steepest descent algorithm only supports one objective function!")
        if self.constraints.size() > 0:
            raise RuntimeError("Steepest descent algorithm does not allow for any constraints!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.model_part_controller.ImportOptimizationModelPart()
        self.model_part_controller.InitializeMeshController()

        self.analyzer.InitializeBeforeOptimizationLoop()

        self.design_surface = self.model_part_controller.GetDesignSurface()

        self.mapper = mapper_factory.CreateMapper(self.design_surface, self.design_surface, self.mapper_settings)
        self.mapper.Initialize()

        self.data_logger = data_logger_factory.CreateDataLogger(self.model_part_controller, self.communicator, self.optimization_settings)
        self.data_logger.InitializeDataLogging()

        self.optimization_utilities = OptimizationUtilities(self.design_surface, self.optimization_settings)

        # Identify fixed design areas (geometric constraints)
        VariableUtils().SetFlag(KM.BOUNDARY, False, self.optimization_model_part.Nodes)
        for itr in range(self.algorithm_settings["fix_boundaries"].size()):
            sub_model_part_name = self.algorithm_settings["fix_boundaries"][itr].GetString()
            node_set = self.optimization_model_part.GetSubModelPart(sub_model_part_name).Nodes
            VariableUtils().SetFlag(KM.BOUNDARY, True, node_set)

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        for self.optimization_iteration in range(1,self.max_iterations):
            print("\n>===================================================================")
            print("> ",timer.GetTimeStamp(),": Starting optimization iteration ",self.optimization_iteration)
            print(">===================================================================\n")

            timer.StartNewLap()

            self.__initializeNewShape()

            self.__analyzeShape()

            if self.line_search_type == "adaptive_stepping" and self.optimization_iteration > 1:
                self.__adjustStepSize()

            self.__computeShapeUpdate()

            self.__logCurrentOptimizationStep()

            print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
            print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

            if self.__isAlgorithmConverged():
                break
            else:
                self.__determineAbsoluteChanges()

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.data_logger.FinalizeDataLogging()
        self.analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeNewShape(self):
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)
        self.model_part_controller.UpdateMeshAccordingInputVariable(KSO.SHAPE_UPDATE)
        self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyzeShape(self):
        self.communicator.initializeCommunication()
        self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
        self.communicator.requestGradientOf(self.objectives[0]["identifier"].GetString())

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.design_surface, self.optimization_iteration, self.communicator)

        objGradientDict = self.communicator.getStandardizedGradient(self.objectives[0]["identifier"].GetString())
        WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, KSO.DF1DX)

        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ComputeUnitSurfaceNormals()
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(KSO.DF1DX)

        self.model_part_controller.DampNodalVariableIfSpecified(KSO.DF1DX)

    # --------------------------------------------------------------------------
    def __adjustStepSize(self):
        current_a = self.step_size

        # Compare actual and estimated improvement using linear information
        dfda1 = 0.0
        for node in self.design_surface.Nodes:
            s1 = node.GetSolutionStepValue(KSO.SEARCH_DIRECTION)
            dfds1 = node.GetSolutionStepValue(KSO.DF1DX_MAPPED)
            dfda1 = dfda1 + s1[0]*dfds1[0] + s1[1]*dfds1[1] + s1[2]*dfds1[2]

        f2 = self.communicator.getStandardizedValue(self.objectives[0]["identifier"].GetString())
        f1 = self.previos_objective_value

        df_actual = f2 - f1
        df_estimated = current_a*dfda1

        # Adjust step size if necessary
        if f2 < f1:
            estimation_error = (df_actual-df_estimated)/df_actual

            # Increase step size if linear approximation matches the actual improvement within a specified tolerance
            if estimation_error < self.approximation_tolerance:
                new_a = min(current_a*self.increase_factor, self.max_step_size)

            # Leave step size unchanged if a nonliner change in the objective is observed but still a descent direction is obtained
            else:
                new_a = current_a
        else:
            # Search approximation of optimal step using interpolation
            a = current_a
            corrected_step_size = - 0.5 * dfda1 * a**2 / (f2 - f1 - dfda1 * a )

            # Starting from the new design, and assuming an opposite gradient direction, the step size to the approximated optimum behaves reciprocal
            new_a = current_a-corrected_step_size

        self.step_size = new_a

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):
        if self.update_mapping_matrix:
            self.mapper.Update()
        self.mapper.InverseMap(KSO.DF1DX, KSO.DF1DX_MAPPED)





        # Perform projection
        new_timer = Timer()
        new_timer.StartTimer()

        num_design_nodes = self.design_surface.NumberOfNodes()
        dJds = np.zeros(num_design_nodes*3)
        num_boundary_nodes = 0

        for itr, node in enumerate(self.design_surface.Nodes):
            temp_vec = node.GetSolutionStepValue(KSO.DF1DX_MAPPED)
            dJds[3*itr+0] = temp_vec[0]
            dJds[3*itr+1] = temp_vec[1]
            dJds[3*itr+2] = temp_vec[2]

            if node.Is(KM.BOUNDARY):
                num_boundary_nodes += 1

        Cm = np.zeros((num_boundary_nodes*3,num_design_nodes*3))

        boundar_node_index = 0

        for row_itr, node in enumerate(self.design_surface.Nodes):
            if node.Is(KM.BOUNDARY):
                node.SetSolutionStepValue(KSO.VECTOR_VARIABLE,[1,1,1])
                self.mapper.InverseMap(KSO.VECTOR_VARIABLE, KSO.VECTOR_VARIABLE_MAPPED)

                node.SetSolutionStepValue(KSO.VECTOR_VARIABLE,[0,0,0])

                for coll_itr, coll_node in enumerate(self.design_surface.Nodes):
                    another_temp_vec = coll_node.GetSolutionStepValue(KSO.VECTOR_VARIABLE_MAPPED)

                    Cm[3*boundar_node_index+0,3*coll_itr+0] = another_temp_vec[0]
                    Cm[3*boundar_node_index+1,3*coll_itr+1] = another_temp_vec[1]
                    Cm[3*boundar_node_index+2,3*coll_itr+2] = another_temp_vec[2]

                boundar_node_index += 1

        print("\n> Time needed for preparation of projection vectors and matrices = ", new_timer.GetLapTime(), "s")

        new_timer.StartNewLap()

        Cm_transpose = Cm.transpose()

        print("\n> Time needed for computation of transpose = ", new_timer.GetLapTime(), "s")

        new_timer.StartNewLap()

        lambda_fac = np.linalg.solve(np.dot(Cm,Cm_transpose), np.dot(Cm,dJds))
        # lambda_fac = np.linalg.lstsq(Cm_transpose, dJds)[0]

        # plot lambda
        boundary_node_itr = 0
        for node in self.design_surface.Nodes:
            if node.Is(KM.BOUNDARY):
                node.SetSolutionStepValue(KM.NODAL_VAUX,lambda_fac[(3*boundary_node_itr):(3*boundary_node_itr+3)].tolist())
                boundary_node_itr += 1

        print("\n> Time needed for solution of projection equation = ", new_timer.GetLapTime(), "s")

        new_timer.StartNewLap()

        p = -1 * ( dJds - np.dot(Cm_transpose, lambda_fac) )

        print("\n> Time needed for compuation of p = ", new_timer.GetLapTime(), "s")

        max_norm = np.linalg.norm(p,np.inf)
        if self.normalize_search_direction:
            p = p/max_norm

        # Store computed search direction (inverse of descent direction)
        for itr, node in enumerate(self.design_surface.Nodes):
            temp_vec = [p[3*itr+0], p[3*itr+1], p[3*itr+2]]
            node.SetSolutionStepValue(KSO.SEARCH_DIRECTION,temp_vec)





        # self.optimization_utilities.ComputeSearchDirectionSteepestDescent()
        self.optimization_utilities.ComputeControlPointUpdate(self.step_size)

        self.mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
        self.model_part_controller.DampNodalVariableIfSpecified(KSO.SHAPE_UPDATE)

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        self.previos_objective_value = self.communicator.getStandardizedValue(self.objectives[0]["identifier"].GetString())
        self.norm_obj_gradient = self.optimization_utilities.ComputeL2NormOfNodalVariable(KSO.DF1DX_MAPPED)

        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.step_size
        additional_values_to_log["norm_obj_gradient"] = self.norm_obj_gradient
        self.data_logger.LogCurrentValues(self.optimization_iteration, additional_values_to_log)
        self.data_logger.LogCurrentDesign(self.optimization_iteration)

    # --------------------------------------------------------------------------
    def __isAlgorithmConverged(self):

        if self.optimization_iteration > 1 :
            # Check if maximum iterations were reached
            if self.optimization_iteration == self.max_iterations:
                print("\n> Maximal iterations of optimization problem reached!")
                return True

            # Check gradient norm
            if self.optimization_iteration == 2:
                self.initial_norm_obj_gradient = self.norm_obj_gradient
            else:
                if self.norm_obj_gradient < self.gradient_tolerance*self.initial_norm_obj_gradient:
                    print("\n> Optimization problem converged as gradient norm reached specified tolerance of ",self.gradient_tolerance)
                    return True

            # Check for relative tolerance
            relativeChangeOfObjectiveValue = self.data_logger.GetValue("rel_change_obj", self.optimization_iteration)
            if abs(relativeChangeOfObjectiveValue) < self.relative_tolerance:
                print("\n> Optimization problem converged within a relative objective tolerance of ",self.relative_tolerance,"%.")
                return True

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges(self):
        self.optimization_utilities.AddFirstVariableToSecondVariable(KSO.CONTROL_POINT_UPDATE, KSO.CONTROL_POINT_CHANGE)
        self.optimization_utilities.AddFirstVariableToSecondVariable(KSO.SHAPE_UPDATE, KSO.SHAPE_CHANGE)

# ==============================================================================
