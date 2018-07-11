from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping

# Importing the base class
from co_simulation_base_io import CoSimulationBaseIO

# Other imports
from co_simulation_tools import csprint, bold, green, red

def Create(io_settings, solvers, solver_name, cosim_solver_details, level):
    return KratosFieldIO(io_settings, solvers, solver_name, cosim_solver_details, level)

class KratosFieldIO(CoSimulationBaseIO):
    def __init__(self, settings, solvers, solver_name, cosim_solver_details, level):
        super(KratosFieldIO, self).__init__(settings, solvers, solver_name, cosim_solver_details, level)
        # Setting up the Mappers

        input_data_list = self.cosim_solver_details[self.solver_name]["input_data_list"]
        output_data_list = self.cosim_solver_details[self.solver_name]["output_data_list"]

        self.mappers = {}
        self.mapping_options = {}
        self.mapper_geometries_map = {}
        self.mapper_flags = {
            "add_values" : KratosMapping.Mapper.ADD_VALUES,
            "swap_sign" : KratosMapping.Mapper.SWAP_SIGN,
            "conservative" : KratosMapping.Mapper.CONSERVATIVE,
        }

        for input_data in input_data_list:
            if not "io_settings" in input_data: # skip if no IO is defined
                continue
            self.__SetupMapper(input_data, input_data["from_solver"], self.solver_name, False)


        for output_data in output_data_list:
            if not "io_settings" in output_data: # skip if no IO is defined
                continue
            self.__SetupMapper(output_data, self.solver_name, output_data["to_solver"], True)


    def __SetupMapper(self, data_entry, from_solver_name, to_solver_name, inverse_map):
        from_solver = self.solvers[from_solver_name]
        to_solver = self.solvers[to_solver_name]

        data_name = data_entry["data_name"]

        data_definition_from = from_solver.GetDataDefinition(data_entry["data_name"])
        data_definition_to = to_solver.GetDataDefinition(data_entry["data_name"])

        geometry_name_from = data_definition_from["geometry_name"]
        geometry_name_to = data_definition_to["geometry_name"]

        # First we check if a mapper for the current geometries exists already
        # Adding the names of the solver to avoid conflicts if geometry names are the same
        if inverse_map:
            name_origin = to_solver_name + "_" + geometry_name_to
            name_destination = from_solver_name + "_" + geometry_name_from
        else:
            name_origin = from_solver_name + "_" + geometry_name_from
            name_destination = to_solver_name + "_" + geometry_name_to

        mapper_exists_already = False
        if name_origin in self.mapper_geometries_map: # if a mapper mapping from this geometry exists already
            if name_destination in self.mapper_geometries_map[name_origin]:
                self.mappers[data_name] = self.mapper_geometries_map[name_origin][name_destination]
                # TODO check also if the mapper type is the same!
                mapper_exists_already = True

        if not mapper_exists_already:
            from_solver_mesh = from_solver.model[geometry_name_from]
            to_solver_mesh = to_solver.model[geometry_name_to]

            mapper_settings = KratosMultiphysics.Parameters("""{
                "mapper_type" : ""
            }""")
            mapper_settings["mapper_type"].SetString(data_entry["io_settings"]["mapper_type"])
            self.mappers[data_name] = KratosMapping.MapperFactory.CreateMapper(from_solver_mesh, to_solver_mesh, mapper_settings)

            # Adding the mapper to the map, in case the same mapper is needed again later (=> same geometry)
            mapper_name_orig = from_solver_name + "_" + geometry_name_from
            mapper_name_dest = to_solver_name + "_" + geometry_name_to
            if not mapper_name_orig in self.mapper_geometries_map:
                self.mapper_geometries_map[mapper_name_orig] = {}
            self.mapper_geometries_map[mapper_name_orig][mapper_name_dest] = self.mappers[data_name]

        if inverse_map: # Swiching bcs the arguments for mapping don't change
            dest_var = KratosMultiphysics.KratosGlobals.GetVariable(data_definition_from["data_identifier"])
            orig_var = KratosMultiphysics.KratosGlobals.GetVariable(data_definition_to["data_identifier"])
        else:
            orig_var = KratosMultiphysics.KratosGlobals.GetVariable(data_definition_from["data_identifier"])
            dest_var = KratosMultiphysics.KratosGlobals.GetVariable(data_definition_to["data_identifier"])

        mapping_flags = KratosMultiphysics.Flags()
        if "mapper_args" in data_entry["io_settings"]:
            for flag_name in data_entry["io_settings"]["mapper_args"]:
                mapping_flags |= self.mapper_flags[flag_name]

        self.mapping_options[data_name] = [orig_var, dest_var, mapping_flags]

        # Printing information related to mapping
        if self.echo_level > 2:
            if mapper_exists_already:
                info_msg  = bold("Existing mapper used") + ' for solver "' + self.solver_name + '": from "'
            else:
                info_msg  = bold("Mapper created") + ' for solver "' + self.solver_name + '": from "'
            info_msg += from_solver_name + ':' + geometry_name_from + '" to "'
            info_msg += to_solver_name + ':' + geometry_name_to + '"'
            csprint(self.lvl, info_msg)

        if self.echo_level > 3:
            info_msg  = bold("Origin_Variable: ") + orig_var.Name() + " | "
            info_msg += bold("Destination_Variable: ") + dest_var.Name()
            if "mapper_args" in data_entry["io_settings"]:
                info_msg += " | with " + bold("Flags") + ": "
                for flag_name in data_entry["io_settings"]["mapper_args"]:
                    info_msg += flag_name + " "
            csprint(self.lvl, info_msg)


    def ImportData(self, data_name, from_client):
        if data_name in self.mappers: # check if this IO-instance has the mapper for the data
            mapper = self.mappers[data_name]
            to_solver = self.solvers[self.solver_name]

            data_definition_from = from_client.GetDataDefinition(data_name)
            data_definition_to = to_solver.GetDataDefinition(data_name)

            orig_var = self.mapping_options[data_name][0]
            dest_var = self.mapping_options[data_name][1]
            flags    = self.mapping_options[data_name][2]

            if self.echo_level > 3:
                info_msg  = bold("Mapping with: ")
                info_msg += bold("Origin_Variable: ") + orig_var.Name() + " | "
                info_msg += bold("Destination_Variable: ") + dest_var.Name()
                csprint(self.lvl, info_msg)

            mapper.Map(orig_var, dest_var, flags)

    def ImportMesh(self, mesh_name, from_client):
        pass

    def ExportData(self, data_name, to_client):
        if data_name in self.mappers: # check if this IO-instance has the mapper for the data
            mapper = self.mappers[data_name]
            from_solver = self.solvers[self.solver_name]

            data_definition_from = from_solver.GetDataDefinition(data_name)
            data_definition_to = to_client.GetDataDefinition(data_name)

            orig_var = self.mapping_options[data_name][0]
            dest_var = self.mapping_options[data_name][1]
            flags    = self.mapping_options[data_name][2]

            if self.echo_level > 3:
                info_msg  = bold("Inverse-Mapping with: ")
                info_msg += bold("Origin_Variable: ") + orig_var.Name() + " | "
                info_msg += bold("Destination_Variable: ") + dest_var.Name()
                csprint(self.lvl, info_msg)

            mapper.InverseMap(orig_var, dest_var, flags)

    def ExportMesh(self, mesh_name, to_client):
        pass

    def MakeDataAvailable(self, data_name, to_client):
        pass
    def MakeMeshAvailable(self, mesh_name, to_client):
        pass
