from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping

class CoSimulationBaseIO(object):
    def __init__(self, solvers, solver_name, cosim_solver_details, level):
        self.solvers = solvers
        self.solver_name = solver_name
        self.cosim_solver_details = cosim_solver_details
        self.lvl = level
        self.echo_level = 0

    def ImportData(self, data_name, from_client):
        pass
    def ImportMesh(self, mesh_name, from_client):
        pass

    def ExportData(self, data_name, to_client):
        pass
    def ExportMesh(self, mesh_name, to_client):
        pass

    def SetEchoLevel(self, level):
        self.echo_level = level

    def PrintInfo(self):
        print("IO does not yet implement PrintInfo!")

    def Check(self):
        print("IO does not yet implement Check!")
