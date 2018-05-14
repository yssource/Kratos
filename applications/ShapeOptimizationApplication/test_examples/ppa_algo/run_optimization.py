from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# For time measures
import time as timer

import os
import shutil

# ======================================================================================================================================
# Model part & solver
# ======================================================================================================================================

# Read parameters


jsonFileName = sys.argv[1]
with open(jsonFileName,'r') as parameter_file:
    ProjectParameters = Parameters(parameter_file.read())

# change cwd, and clean temp folder
shutil.copy(jsonFileName,"temp/param.json")
os.chdir("temp")
cwd = os.getcwd()
files = os.listdir(cwd)
for f in files:
    if f.endswith((".bin",".time",".lst")):
        os.remove(os.path.join(cwd, f))

# Defining the model_part
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

# Create an optimizer
# Note that internally variables related to the optimizer are added to the model part
optimizerFactory = __import__("optimizer_factory")
optimizer = optimizerFactory.CreateOptimizer(ProjectParameters["optimization_settings"], main_model_part)

# Create solver for all response functions specified in the optimization settings
# Note that internally variables related to the individual functions are added to the model part
responseFunctionFactory = __import__("response_function_factory")
listOfResponseFunctions = responseFunctionFactory.CreateListOfResponseFunctions(ProjectParameters["optimization_settings"], main_model_part)

# Create structural solver
# Note that internally variables related to the individual functions are added to the model part
csm_analysis = __import__("structural_mechanics_analysis").StructuralMechanicsAnalysis(ProjectParameters, main_model_part)

# ======================================================================================================================================
# Analyzer
# ======================================================================================================================================

class kratosCSMAnalyzer( (__import__("analyzer_base")).analyzerBaseClass ):

    # --------------------------------------------------------------------------
    def initializeBeforeOptimizationLoop( self ):
        csm_analysis.Initialize()

    # --------------------------------------------------------------------------
    def analyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):

        responseIds = list(communicator.list_of_requests.keys())   # ["strain_energy","mass","eigenfrequency","sphere","translation", several rosenbrock ...]
        solveStructureIds = ["strain_energy","eigenfrequency"]

        noPAsArgumentIds = ["strain_energy","eigenfrequency","mass"]
        if hasattr(communicator,"p"): # projected position_algo
            p = communicator.p

        for responseId in solveStructureIds:
            if communicator.isRequestingValueOf(responseId):
                csm_analysis.InitializeTimeStep()
                csm_analysis.SolveTimeStep()
                csm_analysis.FinalizeTimeStep()
                break

        for responseId in responseIds:
            if communicator.isRequestingValueOf(responseId):
                if responseId in noPAsArgumentIds:
                    listOfResponseFunctions[responseId].CalculateValue()
                else:
                    listOfResponseFunctions[responseId].CalculateValue(p)
                communicator.reportValue(responseId, listOfResponseFunctions[responseId].GetValue())

            if communicator.isRequestingGradientOf(responseId):
                if responseId in noPAsArgumentIds:
                    listOfResponseFunctions[responseId].CalculateGradient()
                else:
                    listOfResponseFunctions[responseId].CalculateGradient(p)
                communicator.reportGradient(responseId, listOfResponseFunctions[responseId].GetGradient())


    # --------------------------------------------------------------------------
    def finalizeAfterOptimizationLoop( self ):
        csm_analysis.Finalize()

    # --------------------------------------------------------------------------

structureAnalyzer = kratosCSMAnalyzer()

# ======================================================================================================================================
# Optimization
# ======================================================================================================================================

optimizer.importAnalyzer( structureAnalyzer )
optimizer.optimize()

# ======================================================================================================================================