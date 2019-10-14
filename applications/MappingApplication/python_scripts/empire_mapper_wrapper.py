import KratosMultiphysics as KM
from KratosMultiphysics.MappingApplication.python_mapper import PythonMapper

import os
import ctypes as ctp

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
            KM.Logger.PrintInfo("EmpireMapperWrapper", "Loading mapper lib ...")
            EmpireMapperWrapper.mapper_lib = self.__LoadMapperLib()



        empire_mapper_type =  self.mapper_settings["mapper_type"].GetString()[7:] # returns the mapper-type without preceeding "empire_"

        # if empire_mapper_type == "nearest_neighbor":

    def __LoadMapperLib(self):
        KM.Logger.PrintInfo("EmpireMapperWrapper", "Determining path to mapper lib")
        # first try automatic detection using the environment that is set by Empire => startEMPIRE
        if ('EMPIRE_MAPPER_LIBSO_ON_MACHINE' in os.environ):
            KM.Logger.PrintInfo("EmpireMapperWrapper", "EMPIRE_MAPPER_LIBSO_ON_MACHINE found in environment")
            mapper_lib_path = os.environ['EMPIRE_API_LIBSO_ON_MACHINE']

        else:
            KM.Logger.PrintInfo("EmpireMapperWrapper", "EMPIRE_MAPPER_LIBSO_ON_MACHINE NOT found in environment, using manually specified path to load the mapper lib")
            mapper_lib_path = self.mapper_settings["mapper_lib"].GetString()
            if mapper_lib_path == "":
                raise Exception('The automatic detection of the mapper lib failed, the path to the mapper lib has to be specified with "mapper_lib"')

        KM.Logger.PrintInfo("EmpireMapperWrapper", "Attempting to load the mapper lib")
        # TODO check if still both are needed! (does the mapperlib link to MPI?)
        try: # OpenMPI
            loaded_mapper_lib = ctp.CDLL(mapper_lib_path, ctp.RTLD_GLOBAL)
            KM.Logger.PrintInfo('EmpireMapperWrapper', 'Using standard OpenMPI')
        except: # Intel MPI or OpenMPI compiled with "–disable-dlopen" option
            loaded_mapper_lib = ctp.cdll.LoadLibrary(mapper_lib_path)
            KM.Logger.PrintInfo('EmpireMapperWrapper', 'Using Intel MPI or OpenMPI compiled with "–disable-dlopen" option')

        KM.Logger.PrintInfo("EmpireMapperWrapper", "Successfully loaded the mapper lib")

        return loaded_mapper_lib


    def __del__(self):
        EmpireMapperWrapper.mapper_count -= 1
