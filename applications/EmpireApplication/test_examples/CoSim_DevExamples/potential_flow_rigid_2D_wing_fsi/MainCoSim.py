from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import sys
import trace



# Importing the Kratos Library => only needed to get the CoSim-Python-Scripts on the path!
import KratosMultiphysics
import KratosMultiphysics.EmpireApplication

# Importing the base class
from co_simulation_steady_analysis import CoSimulationSteadyAnalysis
import json

#parameter_file_name = "project_parameters_cosim_pure_SDoF.json"
#parameter_file_name = "project_parameters_cosim_pure_fluid.json"
parameter_file_name = "project_parameters_cosim_potential_flow_wing_fsi.json"

with open(parameter_file_name, 'r') as parameter_file:
    cosim_parameters = json.load(parameter_file)

simulation = CoSimulationSteadyAnalysis(cosim_parameters)
simulation.Run()

#TODO: remove tracer before merging
##### tracing (for debugging)
# create a Trace object, telling it what to ignore, and whether to
# do tracing or line-counting or both.
#tracer = trace.Trace(ignoredirs=[sys.prefix, sys.exec_prefix], trace=1, count=0)

# run the new command using the given tracer
#tracer.run('main()')

# make a report, placing output in the current directory
#r = tracer.results()
#r.write_results(show_missing=True, coverdir="/home/mlukash/software/Kratos/applications/EmpireApplication/test_examples/CoSim_DevExamples/potential_flow_rigid_2D_wing_fsi/tracing")
