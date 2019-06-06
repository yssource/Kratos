from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.EmpireApplication

# Importing the base class
from co_simulation_analysis import CoSimulationAnalysis
import json

# read parameter file(s)
parameter_file_name = "project_parameters_cosim_oscilating_cylinder_fsi.json"
with open(parameter_file_name, 'r') as parameter_file:
    cosim_parameters = json.load(parameter_file)


simulation = CoSimulationAnalysis(cosim_parameters)

simulation.Run()
