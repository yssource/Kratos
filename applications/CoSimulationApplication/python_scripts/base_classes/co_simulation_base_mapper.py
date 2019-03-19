# co simulation imports
import co_simulation_tools as tools
# Other imports
cs_data_structure = tools.cs_data_structure

##
#  IMPORTANT : This is a BASE CLASS
#               Please do not change anything in this class.
#
#  This class is intended to server as the base class for all the coupled solvers.
class CoSimulationBaseMapper(object):
    def __init__(self, origin_geo, destination_geo, settings):
        self.origin_geo = origin_geo
        self.destination_geo = destination_geo

    ## Map :  Maps the origin_geo to destination_geo
    #
    #  @param from_data                 The origin data to map from
    #  @param to_data                   The destination data to map to
    #  @param mapper_flags              Additional flags for the mapper
    def Map(self, from_data, to_data, mapper_flags):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseMapper : Calling Map function !" + tools.bcolors.ENDC)


    ## InverseMap :  Maps data from destination_geo to the origin_geo
    #
    #  @param from_data                 The origin data to map from
    #  @param to_data                   The destination data to map to
    #  @param mapper_flags              Additional flags for the mapper
    def InverseMap(self, from_data, to_data, mapper_flags):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseMapper : Calling InverseMap function !" + tools.bcolors.ENDC)

    ## UpdateOriginGeometry :  Updates the Origin geometry and re computes the mapping.
    #
    #  @param new_origin_geo            The new geometry with will act as origin geometry.
    def UpdateOriginGeometry(self, new_origin_geo):
        self.origin_geo = new_origin_geo

    ## UpdateDestinationGeometry :  Updates the Origin geometry and re computes the mapping.
    #
    #  @param new_destination_geo       The new geometry with will act as destination geometry.
    def UpdateDestinationGeometry(self, new_destination_geo):
        self.destination_geo = new_destination_geo