from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

import numpy as np

def Create(io_settings):
    return KratosSignalIO(io_settings)

class KratosSignalIO(object):
    def __init__(self, settings):
        pass

    def ImportData(self, data_name, geometry_name, from_client):
        # TODO check if var in ModelPart!
        model_part = from_client.model[geometry_name]
        kratos_var = KratosMultiphysics.KratosGlobals.GetVariable(data_name)

        if type(kratos_var) == KratosMultiphysics.DoubleVariable or type(kratos_var) == KratosMultiphysics.Array1DComponentVariable:
            data = ExtractData(model_part, kratos_var)
        elif type(kratos_var) == KratosMultiphysics.Array1DVariable3:
            domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            if not domain_size in [1,2,3]:
                raise Exception("DOMAIN_SIZE has to be 1, 2 or 3!")
            num_nodes = model_part.NumberOfNodes()
            data = np.zeros(num_nodes*domain_size)
            ext = ["_X", "_Y", "_Z"]
            for i in range(domain_size):
                component_var = KratosMultiphysics.KratosGlobals.GetVariable(kratos_var.Name()+ext[i])
                range_begin = i*num_nodes
                range_end = (i+1)*num_nodes
                data[range_begin:range_end] = ExtractData(model_part, component_var)
        else:
            err_msg  = 'Type of variable "' + kratos_var.Name() + '" is not valid\n'
            err_msg += 'It can only be double, component or array3d!'
            raise Exception(err_msg)

        return data

    def ImportMesh(self, MeshName, FromClient):
        pass

    def ExportData(self, DataName, ToClient):
        pass
    def ExportMesh(self, MeshName, ToClient):
        pass

    def MakeDataAvailable(self, DataName, ToClient):
        pass
    def MakeMeshAvailable(self, MeshName, ToClient):
        pass

def ExtractData(model_part, kratos_var):
    num_nodes = model_part.NumberOfNodes()
    data = np.zeros(num_nodes)
    idx = 0
    for node in model_part.Nodes:
        data[idx] = node.GetSolutionStepValue(kratos_var)
        idx += 1
    return data