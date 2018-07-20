from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.EmpireApplication

# Importing the base class
from co_simulation_solvers.co_simulation_base_solver import CoSimulationBaseSolver

# Other imports
from co_simulation_tools import csprint, yellow


class KratosEmpireSolver(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings, level):
        super(KratosEmpireSolver, self).__init__(cosim_solver_settings, level)
        self.model = KratosMultiphysics.Model()
        self.client_model_part = KratosMultiphysics.ModelPart("ClientModelPart")
        self.model.AddModelPart(self.client_model_part)

        input_file_name = self.cosim_solver_settings["input_file"]
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as parameter_file:
            parameters = json.load(parameter_file)

        self.client_xml_file_name = parameters["client_xml_file_name"]
        self.time_step = parameters["time_step"] # dummy

    def Initialize(self):
        csprint(self.lvl, yellow(self.__class__.__name__ + ":") + " Starting to initialize Empire")
        import empire_wrapper
        csprint(self.lvl, yellow(self.__class__.__name__ + ":") + " Wrapper-Import Successful")
        empire = empire_wrapper.EmpireWrapper()
        csprint(self.lvl, yellow(self.__class__.__name__ + ":") + " Wrapper Created")
        empire.Connect(self.client_xml_file_name)

        empire.ReceiveMesh("client_mesh", self.client_model_part)

    def Finalize(self):
        self.empire.Disconnect()

    def AdvanceInTime(self, current_time):
        self.step_is_repeated = False
        return current_time + self.delta_time # needed for some checks only

    def ImportData(self, data_name, from_client):
        '''This function first receives the data from Empire
        Then it calls the BaseClass, i.e. the IO to do sth with
        the data, e.g. mapping
        '''
        # the following is needed for iterative coupling to tell the client
        # to repeat the same timestep
        if self.step_is_repeated:
            self.empire.SendConvergenceSignal(0) # 0 means that it is NOT converged
        self.step_is_repeated = True

        # receivedata

        super(KratosEmpireSolver, self).ImportData(data_name, from_client)

    def ExportData(self, data_name, to_client):
        '''This function first calls the BaseClass, i.e. the IO to
        do sth with the data, e.g. mapping
        Then it calls Empire to send the data
        '''
        super(KratosEmpireSolver, self).ExportData(data_name, to_client)
        # senddata

    def FinalizeSolutionStep(self):
        '''If this function is called then it means that either convergence is achieved
        or the maximum number of iterations is achieved (in case of iterative coupling)
        In either way and also for explicit coupling, here the connected client is
        informed to advance in time
        '''
        self.empire.SendConvergenceSignal(1) # 1 means that it IS converged