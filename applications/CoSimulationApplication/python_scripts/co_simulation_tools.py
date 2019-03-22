import warnings
cs_data_structure = None
import math
import numpy as np
## Class contains definition of colors. This is to be used as a struct
#
# Example usage print(bcolors.HEADER + "This is a header in header color" + bcolor.ENDC)
# IMPORTANT : The end of the print statement should always contain bcolor.ENDC
class bcolors:
    HEADER    = '\033[95m'
    BLUE      = '\033[94m'
    GREEN     = '\033[92m'
    MAGENTA   = '\033[96m'
    WARNING   = '\033[93m'
    FAIL      = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'

## PrintInfo : Printing information with a label
#
#  @param label         The label for the print
#  @param args          The arguments to be printed
def PrintInfo(label, *args):
    print(label, " ".join(map(str,args)))

## PrintInfo : Printing a warning with a label
#
#  @param label         The label for the print
#  @param args          The arguments to be printed
def PrintWarning(label, *args):
    print(label, " ".join(map(str,args)))


## ImportDataStructure : Imports the data structure which is specified in the parameters file
#
#  @param parameters_file_name   The JSON file name which contains the settings for the co-simulation
def ImportDataStructure(parameters_file_name):
    import json
    global cs_data_structure
    with open(parameters_file_name,'r') as parameter_file:
        parameters = json.load(parameter_file) # This is for getting the flag for database
        if 'data_structure' in parameters['problem_data'].keys():
            data_structure_name = parameters['problem_data']['data_structure']
        else:
            data_structure_name = 'KratosMultiphysics'

        # Initialize cs_data_structure and import corresponding module
        if cs_data_structure is None:
            if data_structure_name == 'KratosMultiphysics':
                cs_data_structure = __import__('KratosMultiphysics')
            elif data_structure_name == 'pyKratos':
                cs_data_structure = __import__('pyKratos')
            else:
                raise Exception('data_structure needs to be KratosMultiphysics or pyKratos')

    return cs_data_structure


## Class CouplingInterfaceData: Class to hold different properties of the data field contributed in
#                           CoSimulation.
#
class CouplingInterfaceData(object):
    def __init__(self, custom_config, solver):

        default_config = cs_data_structure.Parameters("""
        {
            "name" : "default",
            "dimension" : 0,
            "geometry_name" : "",
            "location_on_mesh":"on_nodes"
        }
        """)
        custom_config.ValidateAndAssignDefaults(default_config)

        self.name = custom_config["name"].GetString()
        self.variable = None
        self.filters = []
        self.solver = solver
        self.dimension = custom_config["dimension"].GetInt()
        self.location_on_mesh = custom_config["location_on_mesh"].GetString()
        self.mesh_name = custom_config["geometry_name"].GetString()
        self.destination_data = None
        self.origin_data = None
        self.mapper_settings = None

    def ApplyFilters(self):
        for filter in self.filters:
            filter.Apply()

    def GetPythonList(self):
        data_mesh = self.solver.model[self.mesh_name]
        data = [0]*len(data_mesh.Nodes)*self.dimension
        data_variable = cs_data_structure.KratosGlobals.GetVariable(self.name)
        node_index = 0
        for node in data_mesh.Nodes:
            data_value = node.GetSolutionStepValue(data_variable,0) #TODO what if non-historical?
            for i in range(self.dimension):
                data[node_index*self.dimension + i] = data_value[i]
            node_index+=1
        return data

    def GetNumpyArray(self):
        return np.asarray(self.GetPythonList(), dtype=np.float64)

    def ApplyUpdateToData(self, update):
        data_mesh = self.solver.model[self.mesh_name]
        data_variable = cs_data_structure.KratosGlobals.GetVariable(self.name)
        node_index = 0
        for node in data_mesh.Nodes: # #TODO: local nodes to also work in MPI?
            updated_value = [0]*self.dimension
            # TODO: aditya the data might also be non-historical => GetValue
            for i in range(self.dimension):
                updated_value[i] = update[node_index*self.dimension + i]

            node.SetSolutionStepValue(data_variable, 0, updated_value)
            node_index += 1


