from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping

# Importing the base class
from co_sim_base_io import CoSimulationBaseIO

def Create(io_settings, solvers, solver_name, cosim_solver_details, level):
    return KratosFieldIO(io_settings, solvers, solver_name, cosim_solver_details, level)

class KratosFieldIO(CoSimulationBaseIO):
    def __init__(self, settings, solvers, solver_name, cosim_solver_details, level):
        super(KratosFieldIO, self).__init__(settings, solvers, solver_name, cosim_solver_details, level)
        # Setting up the Mappers
        print(settings)
        print(cosim_solver_details)
        # KratosMapping.MapperFactory.CreateMapper()

        input_data_list = self.cosim_solver_details[self.solver_name]["input_data_list"]

        print("I am:", self.solver_name)

        print("\n", input_data_list)

        self.mappers = {}
        self.mapping_options = {}

        for input_data in input_data_list:
            if not "io_settings" in input_data: # skip if no mapping is defined
                continue
            data_name = input_data["data_name"]
            print(data_name)
            from_solver = self.solvers[input_data["from_solver"]]
            to_solver = self.solvers[self.solver_name]

            data_definition_from = from_solver.GetDataDefinition(input_data["data_name"])
            data_definition_to = to_solver.GetDataDefinition(input_data["data_name"])

            from_solver_mesh = from_solver.model[data_definition_from["geometry_name"]]
            to_solver_mesh = to_solver.model[data_definition_to["geometry_name"]]

            mapper_settings = KratosMultiphysics.Parameters("""{
                "mapper_type" : ""
            }""")
            mapper_settings["mapper_type"].SetString(input_data["io_settings"]["mapper_type"])
            self.mappers[data_name] = KratosMapping.MapperFactory.CreateMapper(from_solver_mesh, to_solver_mesh, mapper_settings)

            orig_var = KratosMultiphysics.KratosGlobals.GetVariable(data_definition_from["data_identifier"])
            dest_var = KratosMultiphysics.KratosGlobals.GetVariable(data_definition_to["data_identifier"])

            self.mapping_options[data_name] = [orig_var, dest_var]





        print(self.mappers)


        # err



        pass


    def ImportData(self, data_name, from_client):
        if data_name in self.mappers: # check if this IO-instance has the mapper for the data
            mapper = self.mappers[data_name]
            to_solver = self.solvers[self.solver_name]

            data_definition_from = from_client.GetDataDefinition(data_name)
            data_definition_to = to_solver.GetDataDefinition(data_name)

            orig_var = self.mapping_options[data_name][0]
            dest_var = self.mapping_options[data_name][1]
            # options  = self.mapping_options[data_name][2]

            mapper.Map(orig_var, dest_var)#, options)

        pass
    def ImportMesh(self, mesh_name, from_client):
        pass

    def ExportData(self, data_name, to_client):
        if data_name in self.mappers: # check if this IO-instance has the mapper for the data
            mapper = self.mappers[data_name]
            to_solver = self.solvers[self.solver_name]

            data_definition_from = from_client.GetDataDefinition(data_name)
            data_definition_to = to_solver.GetDataDefinition(data_name)

            orig_var = self.mapping_options[data_name][0]
            dest_var = self.mapping_options[data_name][1]
            # options  = self.mapping_options[data_name][2]

            mapper.InverseMap(orig_var, dest_var)#, options)

    def ExportMesh(self, mesh_name, to_client):
        pass

    def MakeDataAvailable(self, data_name, to_client):
        pass
    def MakeMeshAvailable(self, mesh_name, to_client):
        pass
