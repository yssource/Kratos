from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping

def Create(io_settings, solvers, solver_name, cosim_solver_details, level):
    return CoSimulationBaseIO(io_settings, solvers, solver_names, cosim_solver_details, level)

class CoSimulationBaseIO(object):
    def __init__(self, settings, solvers, solver_name, cosim_solver_details, level):
        self.settings = settings
        self.solvers = solvers
        self.solver_name = solver_name
        self.cosim_solver_details = cosim_solver_details
        self.lvl = level
        self.echo_level = 0
        if "echo_level" in self.settings:
            self.echo_level = self.settings["echo_level"]

    def ImportData(self, data_name, from_client):
        pass
    def ImportMesh(self, mesh_name, from_client):
        pass

    def ExportData(self, data_name, to_client):
        pass
    def ExportMesh(self, mesh_name, to_client):
        pass

    def MakeDataAvailable(self, data_name, to_client):
        pass
    def MakeMeshAvailable(self, mesh_name, to_client):
        pass
