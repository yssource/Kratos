import KratosMultiphysics as KM
from KratosMultiphysics.MappingApplication.python_mapper import PythonMapper

import os
import ctypes as ctp

def CreateMapper(model_part_origin, model_part_destination, mapper_settings):
    empire_mapper_type =  mapper_settings["mapper_type"].GetString()[7:] # returns the mapper-type without preceeding "empire_"

    if empire_mapper_type == "nearest_neighbor":
        return EmpireNearestNeighborMapper(model_part_origin, model_part_destination, mapper_settings)

    elif empire_mapper_type == "nearest_element":
        return EmpireNearestElementMapper(model_part_origin, model_part_destination, mapper_settings)

    elif empire_mapper_type == "barycentric":
        return EmpireBaryCentricMapper(model_part_origin, model_part_destination, mapper_settings)

    elif empire_mapper_type == "mortar":
        return EmpireMortarMapper(model_part_origin, model_part_destination, mapper_settings)

    else:
        raise Exception('Mapper "{}" not available in empire!'.format(empire_mapper_type))


class EmpireMapperWrapper(PythonMapper):
    """Wrapper for the mappers of Empire with the same Interface as the Mappers in Kratos
    This wrapper requires the development version of Empire which has the mapperlib exposed separately
    """

    mapper_count = 0
    mapper_lib = None

    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        super(EmpireMapperWrapper, self).__init__(model_part_origin, model_part_destination, mapper_settings)

        EmpireMapperWrapper.mapper_count += 1 # required for identification purposes
        self.__coupling_matrices_are_build = False
        self.__inverse_mapper = None

        if EmpireMapperWrapper.mapper_lib == None: # load it only once
            KM.Logger.PrintInfo("EmpireMapperWrapper", "Loading mapper lib ...")
            EmpireMapperWrapper.mapper_lib = self.__LoadMapperLib()


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
        try:
            try: # OpenMPI
                loaded_mapper_lib = ctp.CDLL(mapper_lib_path, ctp.RTLD_GLOBAL)
                KM.Logger.PrintInfo('EmpireMapperWrapper', 'Using standard OpenMPI')
            except: # Intel MPI or OpenMPI compiled with "–disable-dlopen" option
                loaded_mapper_lib = ctp.cdll.LoadLibrary(mapper_lib_path)
                KM.Logger.PrintInfo('EmpireMapperWrapper', 'Using Intel MPI or OpenMPI compiled with "–disable-dlopen" option')
        except OSError:
            raise Exception("Mapper lib could not be loaded!")

        KM.Logger.PrintInfo("EmpireMapperWrapper", "Successfully loaded the mapper lib")

        return loaded_mapper_lib



    def Map(self, variable_origin, variable_destination, mapper_flags=KM.Flags()):
        if not self.__coupling_matrices_are_build:
            self.__BuildCouplingMatrices()
        raise NotImplementedError('"Map" was not implemented for "{}"'.format(self._ClassName()))

    def InverseMap(self, variable_origin, variable_destination, mapper_flags=KM.Flags()):
        # TODO check if using transpose => conservative

        if self.__inverse_mapper == None:
            self.__CreateInverseMapper()
        self.__inverse_mapper.Map(variable_destination, variable_origin, mapper_flags) # TODO check this!
        raise NotImplementedError('"InverseMap" was not implemented for "{}"'.format(self._ClassName()))

    def UpdateInterface(self):
        self.__BuildCouplingMatrices() # TODO check if this works

    def __BuildCouplingMatrices(self):

        if EmpireMapperWrapper.mapper_lib.hasMapper(ctp.c_char_p(self.name)):
            EmpireMapperWrapper.mapper_lib.buildCouplingMatrices(ctp.c_char_p(self.name))
            self.__coupling_matrices_are_build = True
        else:
            pass # TODO print warning


    def __CreateInverseMapper(self):
        return self.__class__(self.model_part_destination, self.model_part_origin, self.mapper_settings) # TODO check this!


    def __del__(self):
        EmpireMapperWrapper.mapper_count -= 1

        if EmpireMapperWrapper.mapper_count == 0: # last mapper was destoyed
            if self.echo_level > 1:
                KM.Logger.PrintInfo('EmpireMapperWrapper', 'Destroying last instance, deleting all meshes & mappers')
            #  delete everything to make sure nothing is left
            EmpireMapperWrapper.mapper_lib.deleteAllMeshes()
            EmpireMapperWrapper.mapper_lib.deleteAllMappers()

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "mapper_lib"  : ""
        }""")
        this_defaults.AddMissingParameters(super(EmpireMapperWrapper, cls)._GetDefaultSettings())
        return this_defaults


class EmpireNearestNeighborMapper(EmpireMapperWrapper):
    """Wrapper for the nearest neighbor mapper of Empire"""

    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        super(EmpireNearestNeighborMapper, self).__init__(model_part_origin, model_part_destination, mapper_settings)


class EmpireNearestElementMapper(EmpireMapperWrapper):
    """Wrapper for the nearest element mapper of Empire"""

    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        super(EmpireNearestElementMapper, self).__init__(model_part_origin, model_part_destination, mapper_settings)


class EmpireBarycentricMapper(EmpireMapperWrapper):
    """Wrapper for the barycentric mapper of Empire"""

    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        super(EmpireBarycentricMapper, self).__init__(model_part_origin, model_part_destination, mapper_settings)


class EmpireMortarMapper(EmpireMapperWrapper):
    """Wrapper for the mortar mapper of Empire"""

    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        super(EmpireMortarMapper, self).__init__(model_part_origin, model_part_destination, mapper_settings)

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "dual"                 : false,
            "enforce_consistency"  : false,
            "opposite_normals"     : false,
        }""")
        this_defaults.AddMissingParameters(super(EmpireMapperWrapper, cls)._GetDefaultSettings())
        return this_defaults