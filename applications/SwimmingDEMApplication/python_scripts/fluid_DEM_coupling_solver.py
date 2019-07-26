from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics
from python_solver import PythonSolver

# Import applications
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

def CreateSolver(model, custom_settings):
    return FluidDEMSolver(model, custom_settings)

class FluidDEMSolver(FluidSolver):

    def PrepareModelPart(self):
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Replace default elements and conditions
            self._ReplaceElementsAndConditions()
            ## Executes the check and prepare model process
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self.main_model_part.SetBufferSize(self.min_buffer_size)

        KratosMultiphysics.Logger.PrintInfo("FluidSolver", "Model reading finished.")

    ## FluidDEMSolver specific methods.

    def _ReplaceElementsAndConditions(self):
        ## Get number of nodes and domain size
        elem_num_nodes = self._GetElementNumNodes()
        cond_num_nodes = self._GetConditionNumNodes()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ## If there are no elements and/or conditions, default to triangles/tetra meshes to avoid breaking the ReplaceElementsAndConditionsProcess
        ## This only affects the input name (if there are no elements or conditions to replace, nothing is replaced).
        if elem_num_nodes == 0:
            elem_num_nodes = domain_size + 1
        if cond_num_nodes == 0:
            cond_num_nodes = domain_size

        ## Complete the element name
        if (self.element_name is not None):
            new_elem_name = self.element_name + str(int(domain_size)) + "D" + str(int(elem_num_nodes)) + "N"
        else:
            raise Exception("There is no element name. Define the self.element_name string variable in your derived solver.")

        ## Complete the condition name
        if (self.condition_name is not None):
            new_cond_name = self.condition_name + str(int(domain_size)) + "D" + str(int(cond_num_nodes)) + "N"
        else:
            raise Exception("There is no condition name. Define the self.condition_name string variable in your derived solver.")

        ## Set the element and condition names in the Json parameters
        #self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""{}""")
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)

        ## Call the replace elements and conditions process
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

    # def _GetElementNumNodes(self):
    #     if self.main_model_part.NumberOfElements() != 0:
    #         if sys.version_info[0] >= 3: # python3 syntax
    #             element_num_nodes = len(self.main_model_part.Elements.__iter__().__next__().GetNodes())
    #         else: # python2 syntax
    #             element_num_nodes = len(self.main_model_part.Elements.__iter__().next().GetNodes())
    #     else:
    #         element_num_nodes = 0

    #     element_num_nodes = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(element_num_nodes)
    #     return element_num_nodes

    # def _GetConditionNumNodes(self):
    #     if self.main_model_part.NumberOfConditions() != 0:
    #         if sys.version_info[0] >= 3: # python3 syntax
    #             condition_num_nodes = len(self.main_model_part.Conditions.__iter__().__next__().GetNodes())
    #         else: # python2 syntax
    #             condition_num_nodes = len(self.main_model_part.Conditions.__iter__().next().GetNodes())
    #     else:
    #         condition_num_nodes = 0

    #     condition_num_nodes = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(condition_num_nodes)
    #     return condition_num_nodes

    # def _ExecuteCheckAndPrepare(self):
    #     ## Check that the input read has the shape we like
    #     prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
    #     prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
    #     prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])

    #     import check_and_prepare_model_process_fluid
    #     check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()

    # def _ComputeDeltaTime(self):
    #     # Automatic time step computation according to user defined CFL number
    #     if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
    #         delta_time = self.EstimateDeltaTimeUtility.EstimateDt()
    #     # User-defined delta time
    #     else:
    #         delta_time = self.settings["time_stepping"]["time_step"].GetDouble()

    #     return delta_time

    # def _GetAutomaticTimeSteppingUtility(self):
    #     if (self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
    #         EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility2D(self.GetComputingModelPart(),
    #                                                                  self.settings["time_stepping"])
    #     else:
    #         EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility3D(self.GetComputingModelPart(),
    #                                                                  self.settings["time_stepping"])

    #     return EstimateDeltaTimeUtility