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
from algorithm_base import OptimizationAlgorithm
import mapper_factory
import data_logger_factory
from custom_timer import Timer
from custom_variable_utilities import WriteDictionaryDataOnNodalVariable

import numpy as np
from scipy.optimize import minimize

# ==============================================================================
class AlgorithmConjugateGradient(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = Parameters("""
        {
            "name"                  : "conjugate_gradient",
            "max_iterations"        : 100,
            "relative_tolerance"    : 1e-3,
            "update_mapping_matrix" : true,
            "fix_boundaries"        : [],
            "gradient_tolerance"    : 1e-6,
            "line_search" : {
                "line_search_type"           : "manual_stepping",
                "normalize_search_direction" : true,
                "step_size"                  : 1.0
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

        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()
        self.update_mapping_matrix = self.algorithm_settings["update_mapping_matrix"].GetBool()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(SEARCH_DIRECTION)

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
        VariableUtils().SetFlag(BOUNDARY, False, self.optimization_model_part.Nodes)

        radius = self.mapper_settings["filter_radius"].GetDouble()
        search_based_functions = SearchBasedFunctions(self.design_surface)

        for itr in range(self.algorithm_settings["fix_boundaries"].size()):
            sub_model_part_name = self.algorithm_settings["fix_boundaries"][itr].GetString()
            node_set = self.optimization_model_part.GetSubModelPart(sub_model_part_name).Nodes
            search_based_functions.FlagNodesInRadius(node_set, BOUNDARY, radius)

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        self.num_nodes = self.design_surface.NumberOfNodes()
        self.num_dv = 3*self.num_nodes

        X0 = np.zeros(self.num_dv)

        self.optimization_iteration = 1
        self.num_function_calls = 0

        def function_wrapper(X):

            self.num_function_calls = self.num_function_calls+1

            for counter, node in enumerate(self.design_surface.Nodes):
                if node.Is(BOUNDARY):
                    # Enforce geometric constraint
                    node.SetSolutionStepValue(CONTROL_POINT_CHANGE,[0.0, 0.0, 0.0])
                else:
                    [dx,dy,dz] = X[(3*counter):(3*counter+3)]
                    node.SetSolutionStepValue(CONTROL_POINT_CHANGE,[dx,dy,dz])

            self.mapper.Map(CONTROL_POINT_CHANGE, SHAPE_CHANGE)

            self.model_part_controller.DampNodalVariableIfSpecified(SHAPE_CHANGE)

            # Initialize new shape
            for counter, node in enumerate(self.design_surface.Nodes):
                [dx,dy,dz] = node.GetSolutionStepValue(SHAPE_CHANGE)

                # Update Coordinates
                node.X = node.X0 + dx
                node.Y = node.Y0 + dy
                node.Z = node.Z0 + dz

                # Update Reference mesh
                node.X0 = node.X
                node.Y0 = node.Y
                node.Z0 = node.Z

            # Analyze shape
            self.communicator.initializeCommunication()
            self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
            self.analyzer.AnalyzeDesignAndReportToCommunicator(self.design_surface, self.num_function_calls, self.communicator)

            value = self.communicator.getStandardizedValue(self.objectives[0]["identifier"].GetString())

            # Log data
            additional_values_to_log = {}
            additional_values_to_log["norm_obj_gradient"] = 0.0
            additional_values_to_log["outer_iter"] = self.optimization_iteration

            self.data_logger.LogCurrentValues(self.num_function_calls, additional_values_to_log)
            self.data_logger.LogCurrentDesign(self.num_function_calls)

            # Reset mesh udpate
            for counter, node in enumerate(self.design_surface.Nodes):
                [dx,dy,dz] = node.GetSolutionStepValue(SHAPE_CHANGE)

                node.X = node.X - dx
                node.Y = node.Y - dy
                node.Z = node.Z - dz
                node.X0 = node.X
                node.Y0 = node.Y
                node.Z0 = node.Z

            return value

        def gradient_wrapper(X):

            for counter, node in enumerate(self.design_surface.Nodes):
                if node.Is(BOUNDARY):
                    # Enforce geometric constraint
                    node.SetSolutionStepValue(CONTROL_POINT_CHANGE,[0.0, 0.0, 0.0])
                else:
                    [dx,dy,dz] = X[(3*counter):(3*counter+3)]
                    node.SetSolutionStepValue(CONTROL_POINT_CHANGE,[dx,dy,dz])

            self.mapper.Map(CONTROL_POINT_CHANGE, SHAPE_CHANGE)

            self.model_part_controller.DampNodalVariableIfSpecified(SHAPE_CHANGE)

            # Initialize new shape
            for counter, node in enumerate(self.design_surface.Nodes):
                [dx,dy,dz] = node.GetSolutionStepValue(SHAPE_CHANGE)

                # Update Coordinates
                node.X = node.X0 + dx
                node.Y = node.Y0 + dy
                node.Z = node.Z0 + dz

                # Update Reference mesh
                node.X0 = node.X
                node.Y0 = node.Y
                node.Z0 = node.Z

            # Analyze shape
            self.communicator.initializeCommunication()
            self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
            self.communicator.requestGradientOf(self.objectives[0]["identifier"].GetString())
            self.analyzer.AnalyzeDesignAndReportToCommunicator(self.design_surface, self.num_function_calls, self.communicator)

            value = self.communicator.getStandardizedValue(self.objectives[0]["identifier"].GetString())

            gradient_from_kratos = self.communicator.getStandardizedGradient(self.objectives[0]["identifier"].GetString())
            WriteDictionaryDataOnNodalVariable(gradient_from_kratos, self.optimization_model_part, DF1DX)

            # Mapping
            self.model_part_controller.DampNodalVariableIfSpecified(DF1DX)
            self.mapper.InverseMap(DF1DX, DF1DX_MAPPED)

            gradient = np.zeros_like(X)
            for counter, node in enumerate(self.design_surface.Nodes):
                gradient[(3*counter):(3*counter+3)] = node.GetSolutionStepValue(DF1DX_MAPPED)

            norm_obj_gradient = self.optimization_utilities.ComputeL2NormOfNodalVariable(DF1DX_MAPPED)

            additional_values_to_log = {}
            additional_values_to_log["norm_obj_gradient"] = norm_obj_gradient
            additional_values_to_log["outer_iter"] = self.optimization_iteration

            self.data_logger.LogCurrentValues(self.num_function_calls, additional_values_to_log)
            self.data_logger.LogCurrentDesign(1000+self.num_function_calls)

            # Reset mesh udpate
            for counter, node in enumerate(self.design_surface.Nodes):
                [dx,dy,dz] = node.GetSolutionStepValue(SHAPE_CHANGE)

                node.X = node.X - dx
                node.Y = node.Y - dy
                node.Z = node.Z - dz
                node.X0 = node.X
                node.Y0 = node.Y
                node.Z0 = node.Z

            # Return value
            return gradient

        def log_function(X):
            if self.update_mapping_matrix:
                self.mapper.Update()

            self.data_logger.LogCurrentDesign(3000+self.optimization_iteration)

            self.optimization_iteration = self.optimization_iteration+1
            print("##################################################################### Starting opt step: ",self.optimization_iteration)

        # Run optimization
        res = minimize(function_wrapper, X0, method='CG', jac=gradient_wrapper, callback=log_function, options={'disp': True})

        print(res)

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.data_logger.FinalizeDataLogging()
        self.analyzer.FinalizeAfterOptimizationLoop()

# ==============================================================================
