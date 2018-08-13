from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping

# Importing the base class
from co_simulation_base_io import CoSimulationBaseIO

# Other imports
import numpy as np
from co_simulation_tools import csprint, bold, green, red

def Create(solvers, solver_name, cosim_solver_details, level):
    return KratosIO(solvers, solver_name, cosim_solver_details, level)

class KratosIO(CoSimulationBaseIO):
    def __init__(self, solvers, solver_name, cosim_solver_details, level):
        super(KratosIO, self).__init__(solvers, solver_name, cosim_solver_details, level)
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

    def ImportData(self, data_settings, from_client):
        data_format = data_settings["data_format"]
        data_name = data_settings["data_name"]

        if data_format == "numpy_array":
            # TODO check if var in ModelPart!
            # In this case the from_client is the solver itself
            data_definition = from_client.GetDataDefinition(data_name)
            geometry_name = data_definition["geometry_name"]
            var_name = data_definition["data_identifier"]
            data_array = data_settings["data_array"]

            model_part = from_client.model[geometry_name]
            kratos_var = KratosMultiphysics.KratosGlobals.GetVariable(var_name)

            if type(kratos_var) == KratosMultiphysics.DoubleVariable or type(kratos_var) == KratosMultiphysics.Array1DComponentVariable:
                SetData(model_part, kratos_var, data_array)
            elif type(kratos_var) == KratosMultiphysics.Array1DVariable3:
                domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                if not domain_size in [1,2,3]:
                    raise Exception("DOMAIN_SIZE has to be 1, 2 or 3!")
                num_nodes = NumberOfNodes(model_part)
                if data_array.size != num_nodes*domain_size:
                    raise Exception("Size of data does not match number of nodes x domain size!")
                ext = ["_X", "_Y", "_Z"]
                for i in range(domain_size):
                    component_var = KratosMultiphysics.KratosGlobals.GetVariable(kratos_var.Name()+ext[i])
                    range_begin = i*num_nodes
                    range_end = (i+1)*num_nodes
                    SetData(model_part, component_var, data_array[range_begin:range_end])
            else:
                err_msg  = 'Type of variable "' + kratos_var.Name() + '" is not valid\n'
                err_msg += 'It can only be double, component or array3d!'
                raise Exception(err_msg)

        elif data_format == "kratos_modelpart":
            ## TODO check why the incoming client is not used!
            if data_name in self.mappers: # check if this IO-instance has the mapper for the data
                mapper = self.mappers[data_name]

                orig_var = self.mapping_options[data_name][0]
                dest_var = self.mapping_options[data_name][1]
                flags    = self.mapping_options[data_name][2]

                if self.echo_level > 3:
                    info_msg  = bold("Mapping with: ")
                    info_msg += bold("Origin_Variable: ") + orig_var.Name() + " | "
                    info_msg += bold("Destination_Variable: ") + dest_var.Name()
                    csprint(self.lvl, info_msg)

                mapper.Map(orig_var, dest_var, flags)

        elif data_format == "single_value":
            # TODO check if var in ModelPart!
            # In this case the from_client is the solver itself
            data_definition = from_client.GetDataDefinition(data_name)
            geometry_name = data_definition["geometry_name"]
            var_name = data_definition["data_identifier"]
            value = data_settings["single_value"]

            model_part = from_client.model[geometry_name]
            kratos_var = KratosMultiphysics.KratosGlobals.GetVariable(var_name)

            if type(kratos_var) == KratosMultiphysics.DoubleVariable or type(kratos_var) == KratosMultiphysics.Array1DComponentVariable:
                for node in Nodes(model_part):
                    node.SetSolutionStepValue(kratos_var, value)
            else:
                err_msg  = 'Type of variable "' + kratos_var.Name() + '" is not valid\n'
                err_msg += 'It can only be double, component!'
                raise Exception(err_msg)

        else:
            raise Exception("The requested data_format is not implemented in KratosIO!")


    def ExportData(self, data_settings, to_client):
        data_format = data_settings["data_format"]
        data_name = data_settings["data_name"]

        if data_format == "numpy_array":
            # TODO check if var in ModelPart!
            # In this case the to_client is the solver itself
            data_definition = to_client.GetDataDefinition(data_name)
            geometry_name = data_definition["geometry_name"]
            var_name = data_definition["data_identifier"]
            data_array = data_settings["data_array"]
            buffer_index = data_settings["buffer_index"]

            model_part = to_client.model[geometry_name]
            kratos_var = KratosMultiphysics.KratosGlobals.GetVariable(var_name)

            if type(kratos_var) == KratosMultiphysics.DoubleVariable or type(kratos_var) == KratosMultiphysics.Array1DComponentVariable:
                required_size = NumberOfNodes(model_part)
                if not data_array.size == required_size:
                    # data_array = np.resize(data_array, (1,required_size))
                    data_array.resize(required_size, refcheck=False)
                ExtractData(model_part, kratos_var, data_array, buffer_index)
            elif type(kratos_var) == KratosMultiphysics.Array1DVariable3:
                domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                if not domain_size in [1,2,3]:
                    raise Exception("DOMAIN_SIZE has to be 1, 2 or 3!")

                num_nodes = NumberOfNodes(model_part)
                required_size = num_nodes * domain_size
                if not data_array.size == required_size:
                    # data_array = np.resize(data_array, (1,required_size))
                    data_array.resize(required_size, refcheck=False)

                ext = ["_X", "_Y", "_Z"]
                for i in range(domain_size):
                    component_var = KratosMultiphysics.KratosGlobals.GetVariable(kratos_var.Name()+ext[i])
                    range_begin = i*num_nodes
                    range_end = (i+1)*num_nodes
                    ExtractData(model_part, component_var, data_array[range_begin:range_end], buffer_index)
            else:
                err_msg  = 'Type of variable "' + kratos_var.Name() + '" is not valid\n'
                err_msg += 'It can only be double, component or array3d!'
                raise Exception(err_msg)

        elif data_format == "kratos_modelpart":
            ## TODO check why the incoming client is not used!
            if data_name in self.mappers: # check if this IO-instance has the mapper for the data => I think this is not needed any more
                mapper = self.mappers[data_name]

                orig_var = self.mapping_options[data_name][0]
                dest_var = self.mapping_options[data_name][1]
                flags    = self.mapping_options[data_name][2]

                if self.echo_level > 3:
                    info_msg  = bold("Inverse-Mapping with: ")
                    info_msg += bold("Origin_Variable: ") + orig_var.Name() + " | "
                    info_msg += bold("Destination_Variable: ") + dest_var.Name()
                    csprint(self.lvl, info_msg)

                mapper.InverseMap(orig_var, dest_var, flags)

        else:
            raise Exception("The requested data_format is not implemented in KratosIO!")

def ExtractData(model_part, kratos_var, data_array, buffer_index):
    for i, node in enumerate(Nodes(model_part)):
        data_array[i] = node.GetSolutionStepValue(kratos_var, buffer_index)

def SetData(model_part, kratos_var, data_array):
    num_nodes = NumberOfNodes(model_part)
    if data_array.size != num_nodes:
        raise Exception("Size of data does not match number of nodes!")
    for i, node in enumerate(Nodes(model_part)):
        node.SetSolutionStepValue(kratos_var, data_array[i])


def Nodes(model_part):
    # Wrapper to avoid long call
    return model_part.GetCommunicator().LocalMesh().Nodes

def NumberOfNodes(model_part):
    return len(Nodes(model_part)) # Mesh does currently not expose NumberOfNodes!
