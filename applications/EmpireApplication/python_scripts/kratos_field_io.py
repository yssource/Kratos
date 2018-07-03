from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping

def Create(io_settings):
    return KratosFieldIO(io_settings)

class KratosFieldIO(object):
    def __init__(self, settings):
        pass

    def ImportData(self, data_name, from_client):
        print("IN IMPORTDATA", data_name, from_client)










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