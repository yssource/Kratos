import KratosMultiphysics as KM

class PythonMapper(object):
    """Baseclass for python based mappers in Kratos
    The inteface matches the C++ version ("custom_mappers/mapper.h")
    """
    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        self.model_part_origin = model_part_origin
        self.model_part_destination = model_part_destination
        self.mapper_settings = mapper_settings # Note: no validation done here, should be done in derived class

        self.mapper_settings.ValidateAndAssignDefaults(self._GetDefaultSettings())
        self.echo_level = self.mapper_settings["echo_level"].GetInt()

    def Map(self, variable_origin, variable_destination, mapper_flags=KM.Flags()):
        raise NotImplementedError('"Map" was not implemented for "{}"'.format(self._ClassName()))

    def InverseMap(self, variable_origin, variable_destination, mapper_flags=KM.Flags()):
        raise NotImplementedError('"InverseMap" was not implemented for "{}"'.format(self._ClassName()))

    def UpdateInterface(self):
        raise NotImplementedError('"UpdateInterface" was not implemented for "{}"'.format(self._ClassName()))

    @classmethod
    def _GetDefaultSettings(cls):
        return KM.Parameters("""{
            "mapper_type" : "",
            "echo_level"  : 0
        }""")

    @classmethod
    def _ClassName(cls):
        return cls.__name__
