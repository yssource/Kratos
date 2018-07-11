from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Importing the base class
from co_simulation_base_io import CoSimulationBaseIO

# Other imports
import numpy as np

def Create(io_settings, solvers, solver_name, cosim_solver_details, level):
    return KratosSignalIO(io_settings, solvers, solver_name, cosim_solver_details, level)

class KratosSignalIO(CoSimulationBaseIO):

    def ImportData(self, data_name, from_client):
        # TODO check if var in ModelPart!
        data_definition = from_client.GetDataDefinition(data_name)
        geometry_name = data_definition["geometry_name"]
        var_name = data_definition["data_identifier"]

        model_part = from_client.model[geometry_name]
        kratos_var = KratosMultiphysics.KratosGlobals.GetVariable(var_name)

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

    def ExportData(self, data_name, to_client, data):
        # TODO check if var in ModelPart!
        data_definition = to_client.GetDataDefinition(data_name)
        geometry_name = data_definition["geometry_name"]
        var_name = data_definition["data_identifier"]

        model_part = to_client.model[geometry_name]
        kratos_var = KratosMultiphysics.KratosGlobals.GetVariable(var_name)

        if type(kratos_var) == KratosMultiphysics.DoubleVariable or type(kratos_var) == KratosMultiphysics.Array1DComponentVariable:
            SetData(model_part, kratos_var, data)
        elif type(kratos_var) == KratosMultiphysics.Array1DVariable3:
            domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            if not domain_size in [1,2,3]:
                raise Exception("DOMAIN_SIZE has to be 1, 2 or 3!")
            num_nodes = model_part.NumberOfNodes()
            if data.size != num_nodes*domain_size:
                raise Exception("Size of data does not match number of nodes x domain size!")
            ext = ["_X", "_Y", "_Z"]
            for i in range(domain_size):
                component_var = KratosMultiphysics.KratosGlobals.GetVariable(kratos_var.Name()+ext[i])
                range_begin = i*num_nodes
                range_end = (i+1)*num_nodes
                SetData(model_part, component_var, data[range_begin:range_end])
        else:
            err_msg  = 'Type of variable "' + kratos_var.Name() + '" is not valid\n'
            err_msg += 'It can only be double, component or array3d!'
            raise Exception(err_msg)

    def ExportMesh(self, MeshName, ToClient):
        pass

    def MakeDataAvailable(self, DataName, ToClient):
        pass
    def MakeMeshAvailable(self, MeshName, ToClient):
        pass

def ExtractData(model_part, kratos_var):
    num_nodes = model_part.NumberOfNodes()
    data = np.zeros(num_nodes)
    for idx, node in zip(range(num_nodes), model_part.Nodes):
        data[idx] = node.GetSolutionStepValue(kratos_var)
    return data

def SetData(model_part, kratos_var, data):
    num_nodes = model_part.NumberOfNodes()
    if data.size != num_nodes:
        raise Exception("Size of data does not match number of nodes!")
    for idx, node in zip(range(num_nodes), model_part.Nodes):
        node.SetSolutionStepValue(kratos_var, data[idx])