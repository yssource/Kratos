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
    def __init__(self, OptimizationSettings, Analyzer, Communicator, ModelPartController):
        default_algorithm_settings = Parameters("""
        {
            "name"               : "steepest_descent",
            "max_iterations"     : 100,
            "relative_tolerance" : 1e-3,
            "line_search" : {
                "line_search_type"           : "manual_stepping",
                "normalize_search_direction" : true,
                "step_size"                  : 1.0
            }
        }""")
        self.algorithm_settings =  OptimizationSettings["optimization_algorithm"]
        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

        self.Analyzer = Analyzer
        self.Communicator = Communicator
        self.ModelPartController = ModelPartController

        self.objectives = OptimizationSettings["objectives"]
        self.constraints = OptimizationSettings["constraints"]

        self.OptimizationModelPart = ModelPartController.GetOptimizationModelPart()
        self.DesignSurface = ModelPartController.GetDesignSurface()

        self.Mapper = mapper_factory.CreateMapper(self.DesignSurface, self.DesignSurface, OptimizationSettings["design_variables"]["filter"])
        self.DataLogger = data_logger_factory.CreateDataLogger(ModelPartController, Communicator, OptimizationSettings)


    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Steepest descent algorithm only supports one objective function!")
        if self.constraints.size() > 0:
            raise RuntimeError("Steepest descent algorithm does not allow for any constraints!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.only_obj = self.objectives[0]

        self.maxIterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relativeTolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.ModelPartController.InitializeMeshController()
        self.Mapper.Initialize()
        self.Analyzer.InitializeBeforeOptimizationLoop()
        self.DataLogger.InitializeDataLogging()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        self.num_nodes = self.DesignSurface.NumberOfNodes()
        self.num_dv = 3*self.num_nodes

        X0 = np.zeros(self.num_dv)

        self.optimization_iteration = 0
        self.num_function_calls = 0
        self.num_gradient_calls = 0

        def function_wrapper(X):

            self.num_function_calls = self.num_function_calls+1

            # Mapping
            for counter, node in enumerate(self.DesignSurface.Nodes):
                [dx,dy,dz] = X[(3*counter):(3*counter+3)]
                node.SetSolutionStepValue(CONTROL_POINT_CHANGE,[dx,dy,dz])

            self.Mapper.Map(CONTROL_POINT_CHANGE, SHAPE_CHANGE)

            # Initialize new shape
            for counter, node in enumerate(self.DesignSurface.Nodes):
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
            self.Communicator.initializeCommunication()
            self.Communicator.requestValueOf(self.only_obj["identifier"].GetString())
            self.Communicator.requestGradientOf(self.only_obj["identifier"].GetString())
            self.Analyzer.AnalyzeDesignAndReportToCommunicator(self.DesignSurface, self.optimization_iteration, self.Communicator)

            value = self.Communicator.getStandardizedValue(self.only_obj["identifier"].GetString())

            gradient_from_kratos = self.Communicator.getStandardizedGradient(self.only_obj["identifier"].GetString())
            WriteDictionaryDataOnNodalVariable(gradient_from_kratos, self.OptimizationModelPart, DF1DX)

            # Mapping
            self.Mapper.Update()
            self.Mapper.InverseMap(DF1DX, DF1DX_MAPPED)

            gradient = np.zeros_like(X)
            for counter, node in enumerate(self.DesignSurface.Nodes):
                gradient[(3*counter):(3*counter+3)] = node.GetSolutionStepValue(DF1DX_MAPPED)

            self.DataLogger.LogCurrentValues(self.num_function_calls, {})
            self.DataLogger.LogCurrentDesign(self.num_function_calls)

            # Reset mesh udpate
            for counter, node in enumerate(self.DesignSurface.Nodes):
                [dx,dy,dz] = node.GetSolutionStepValue(SHAPE_CHANGE)

                node.X = node.X - dx
                node.Y = node.Y - dy
                node.Z = node.Z - dz
                node.X0 = node.X
                node.Y0 = node.Y
                node.Z0 = node.Z

            # Return value
            return value, gradient

        def log_function(X):
            self.optimization_iteration = self.optimization_iteration+1

            # self.DataLogger.LogCurrentValues(self.optimization_iteration, {})
            # self.DataLogger.LogCurrentDesign(self.optimization_iteration)

        # Run optimization
        res = minimize(function_wrapper, X0, method='CG', jac=True, callback=log_function, options={"disp": True})

        print(res)






        # for self.optimization_iteration in range(1,self.maxIterations):
        #     print("\n>===================================================================")
        #     print("> ",timer.GetTimeStamp(),": Starting optimization iteration ",self.optimization_iteration)
        #     print(">===================================================================\n")

        #     timer.StartNewLap()

        #     # Initialize new shape
        #     self.ModelPartController.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
        #     self.ModelPartController.SetReferenceMeshToMesh()

        #     # Analyze shape
        #     self.Communicator.initializeCommunication()
        #     self.Communicator.requestValueOf(self.only_obj["identifier"].GetString())
        #     self.Communicator.requestGradientOf(self.only_obj["identifier"].GetString())

        #     self.Analyzer.AnalyzeDesignAndReportToCommunicator(self.DesignSurface, self.optimization_iteration, self.Communicator)

        #     objGradientDict = self.Communicator.getStandardizedGradient(self.only_obj["identifier"].GetString())
        #     WriteDictionaryDataOnNodalVariable(objGradientDict, self.OptimizationModelPart, DF1DX)

        #     if self.only_obj["project_gradient_on_surface_normals"].GetBool():
        #         self.ModelPartController.ComputeUnitSurfaceNormals()
        #         self.ModelPartController.ProjectNodalVariableOnUnitSurfaceNormals(DF1DX)

        #     self.ModelPartController.DampNodalVariableIfSpecified(DF1DX)

        #     # Compute shape update
        #     self.Mapper.Update()
        #     self.Mapper.InverseMap(DF1DX, DF1DX_MAPPED)

        #     self.OptimizationUtilities.ComputeSearchDirectionSteepestDescent()
        #     self.OptimizationUtilities.ComputeControlPointUpdate()

        #     self.Mapper.Map(CONTROL_POINT_UPDATE, SHAPE_UPDATE)

        #     self.ModelPartController.DampNodalVariableIfSpecified(SHAPE_UPDATE)

        #     # Log current optimization step
        #     additional_values_to_log = {}
        #     additional_values_to_log["step_size"] = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        #     self.DataLogger.LogCurrentValues(self.optimization_iteration, additional_values_to_log)
        #     self.DataLogger.LogCurrentDesign(self.optimization_iteration)

        #     print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
        #     print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

        #     # Determine absolute changes
        #     self.OptimizationUtilities.AddFirstVariableToSecondVariable(CONTROL_POINT_UPDATE, CONTROL_POINT_CHANGE)
        #     self.OptimizationUtilities.AddFirstVariableToSecondVariable(SHAPE_UPDATE, SHAPE_CHANGE)

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.DataLogger.FinalizeDataLogging()
        self.Analyzer.FinalizeAfterOptimizationLoop()

# ==============================================================================
