import KratosMultiphysics as KM
from KratosMultiphysics.MappingApplication.python_mapper import PythonMapper

import os

def CreateMapper(model_part_origin, model_part_destination, mapper_settings):
    return EmpireMapperWrapper(model_part_origin, model_part_destination, mapper_settings)

class EmpireMapperWrapper(PythonMapper):
    """Wrapper for the mappers of Empire with the same Interface as the Mappers in Kratos"""

    mapper_count = 0
    mapper_lib = None

    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        super(EmpireMapperWrapper, self).__init__(model_part_origin, model_part_destination, mapper_settings)
        default_parameters = KM.Parameters("""{
            "mapper_type" : "PLEASE_SPECIFY_MAPPER_TYPE",
            "mapper_lib"  : ""
        }""")

        self.mapper_settings.ValidateAndAssignDefaults(default_parameters)

        EmpireMapperWrapper.mapper_count += 1 # required for identification purposes

        if EmpireMapperWrapper.mapper_lib == None: # load it only once
            KM.Logger.PrintInfo("EmpireMapperWrapper", "Attempting to load mapper lib")
            # first try automatic detection using the environment that is set by Empire => startEMPIRE
            if ('EMPIRE_MAPPER_LIBSO_ON_MACHINE' in os.environ):
                KM.Logger.PrintInfo("EmpireMapperWrapper", "EMPIRE_MAPPER_LIBSO_ON_MACHINE found in environment")

            else:
                KM.Logger.PrintInfo("EmpireMapperWrapper", "EMPIRE_MAPPER_LIBSO_ON_MACHINE NOT found in environment, using manually specified path to load the mapper lib")
                mapper_lib_path = self.mapper_settings["mapper_lib"].GetString()
                if mapper_lib_path == "":
                    raise Exception('The automatic detection of the mapper lib failed, the path to the mapper lib has to be specified with "mapper_lib"')

        empire_mapper_type =  self.mapper_settings["mapper_type"].GetString()[7:] # returns the mapper-type without preceeding "empire_"

        # if empire_mapper_type == "nearest_neighbor":

    def __del__(self):
        EmpireMapperWrapper.mapper_count -= 1
