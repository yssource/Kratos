from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# For time measures
import time as timer

# ======================================================================================================================================
# Model part & solver
# ======================================================================================================================================

parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

#defining the model_part
main_model_part = ModelPart(ProjectParameters["optimization_settings"]["design_variables"]["optimization_model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["optimization_settings"]["design_variables"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available in kratos core
Model = {ProjectParameters["optimization_settings"]["design_variables"]["optimization_model_part_name"].GetString() : main_model_part}

# Create an optimizer
# Note that internally variables related to the optimizer are added to the model part
optimizerFactory = __import__("optimizer_factory")
optimizer = optimizerFactory.CreateOptimizer( main_model_part, ProjectParameters["optimization_settings"] )

# ======================================================================================================================================
# Analyzer
# ======================================================================================================================================

class ExternalAnalyzer( (__import__("analyzer_base")).analyzerBaseClass ):

    # --------------------------------------------------------------------------
    def __init__(self):
        self.target_Y = 5

    # --------------------------------------------------------------------------
    def analyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):

        # Calculation of value of strain energy
        if communicator.isRequestingValueOf("distance_to_xz"):

            value = 0.0
            for node in currentDesign.Nodes:
                local_deviation = node.Y - self.target_Y
                value = value+local_deviation**2

            communicator.reportValue("distance_to_xz", value)

        # Calculation of gradient of strain energy
        if communicator.isRequestingGradientOf("distance_to_xz"):

            gradientOnDesignSurface = {}
            for node in currentDesign.Nodes:
                local_deviation = node.Y - self.target_Y
                gradientOnDesignSurface[node.Id] = [0,2*local_deviation,0]


            communicator.reportGradient("distance_to_xz", gradientOnDesignSurface)

    # --------------------------------------------------------------------------

# ======================================================================================================================================
# Optimization
# ======================================================================================================================================

optimizer.importAnalyzer( ExternalAnalyzer() )
optimizer.optimize()

# ======================================================================================================================================