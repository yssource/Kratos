from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid


# Other imports
import co_simulation_tools as cs_tools
import co_simulation_ios.co_simulation_io_factory



def CreateSolver(cosim_solver_settings, level):
    return KratosFluidSolver(cosim_solver_settings, level)

class KratosFluidSolver(object):

    def __init__(self, cosim_solver_settings, level):
        self.cosim_solver_settings = cosim_solver_settings
        self.lvl = level
        self.echo_level = 0
        if "echo_level" in self.cosim_solver_settings:
            self.echo_level = self.cosim_solver_settings["echo_level"]
        self.io_is_initialized = False

        self.model = KratosMultiphysics.Model()

        input_file_name = self.cosim_solver_settings["input_file"]
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        self.project_parameters = KratosMultiphysics.Parameters(parameters)

    def Initialize(self):
        self.parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()

        if (self.parallel_type == "MPI"):
            import KratosMultiphysics.mpi as KratosMPI
            self.is_printing_rank = (KratosMPI.mpi.rank == 0)
        else:
            self.is_printing_rank = True

        self.solver = python_solvers_wrapper_fluid.CreateSolver(self.model, self.project_parameters)
        self.solver.AddVariables()

        self.solver.ImportModelPart()
        self.solver.PrepareModelPart()
        self.solver.AddDofs()

        order_processes_initialization = ["gravity",
                "initial_conditions_process_list",
                "boundary_conditions_process_list",
                "auxiliar_process_list"]
        self._list_of_processes = self._CreateProcesses("processes", order_processes_initialization)
        self._list_of_output_processes = self._CreateProcesses("output_processes", order_processes_initialization)
        self._list_of_processes.extend(self._list_of_output_processes)  # Adding the output processes to the regular processes

        ##here we initialize user-provided processes
        for process in self._list_of_processes:
            process.ExecuteInitialize()

        self.solver.Initialize()
        self.solver.Check()
        for process in self._list_of_processes:
            process.Check()
            process.ExecuteBeforeSolutionLoop()

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        if self.solver.GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.time = self.solver.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()

        ## If the echo level is high enough, print the complete list of settings used to run the simualtion
        if self.is_printing_rank and self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("Fluid Dynamics Analysis", "Analysis -START- ")

    def InitializeIO(self, solvers, io_echo_level):
        solver_name = self.cosim_solver_settings["name"]
        if self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is already initialized!')

        self.io = co_simulation_ios.co_simulation_io_factory.CreateIO(self._GetIOName(),
                                      solvers,
                                      solver_name,
                                      self.lvl)
        self.io.SetEchoLevel(io_echo_level)
        self.io_is_initialized = True

    def Finalize(self):
        for process in self._list_of_processes:
            process.ExecuteFinalize()

        self.solver.Finalize()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("Fluid Dynamics Analysis", "Analysis -END- ")

    def AdvanceInTime(self, current_time):
        new_time = self.solver.AdvanceInTime(current_time)
        self.time = new_time # only needed to print the time correctly
        self.delta_time = new_time - current_time
        return new_time

    def Predict(self):
        self.solver.Predict()

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        for process in self._list_of_processes:
            process.ExecuteInitializeSolutionStep()

        self.solver.InitializeSolutionStep()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("Fluid Dynamics Analysis", "STEP: ",
                                                self.solver.GetComputingModelPart().ProcessInfo[
                                                    KratosMultiphysics.STEP])
            KratosMultiphysics.Logger.PrintInfo("Fluid Dynamics Analysis", "TIME: ", self.time)

    def FinalizeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.
        """
        self.solver.FinalizeSolutionStep()

        for process in self._list_of_processes:
            process.ExecuteFinalizeSolutionStep()

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """
        # first we check if one of the output processes will print output in this step
        # this is done to save computation in case none of them will print
        is_output_step = False
        for output_process in self._list_of_output_processes:
            if output_process.IsOutputStep():
                is_output_step = True
                break

        if is_output_step:  # at least one of the output processes will print output
            for process in self._list_of_processes:
                process.ExecuteBeforeOutputStep()

            for output_process in self._list_of_output_processes:
                if output_process.IsOutputStep():
                    output_process.PrintOutput()

            for process in self._list_of_processes:
                process.ExecuteAfterOutputStep()

    def SolveSolutionStep(self):
        self.solver.SolveSolutionStep()

    def GetBufferSize(self):
        model_part_name = self.project_parameters["solver_settings"]["model_part_name"].GetString()
        return self.model[model_part_name].GetBufferSize()

    def GetDeltaTime(self):
        if not hasattr(self, 'delta_time'):
            raise Exception("DeltaTime can only be querried after it has been computed at least once")
        return self.delta_time

    def SetEchoLevel(self, level):
        self.echo_level = level

    def ImportData(self, data_name, from_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ImportData(data_name, from_client)

    def ImportMesh(self, mesh_name, from_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ImportMesh(mesh_name, from_client)

    def ExportData(self, data_name, to_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ExportData(data_name, to_client)

    def ExportMesh(self, mesh_name, to_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ExportMesh(mesh_name, to_client)

    def GetDataDefinition(self, data_name):
        return self.cosim_solver_settings["data"][data_name]

    def PrintInfo(self):
        cs_tools.solverprint(self.lvl, "KratosSolver", cs_tools.bold(self._Name()))

    ################################### BAUSTELLE ###################################(from analysis_stage)
    ################################### BAUSTELLE ###################################
    def Check(self):
        is_distributed = cs_tools.COSIM_SPACE.IsDistributed()
        if is_distributed and not self.parallel_type == "MPI":
            warning_msg  = 'WARNING: Global "parallel_type" (MPI) is different '
            warning_msg += 'from local one (' + self.parallel_type + ')!'
            cs_tools.solverprint(self.lvl, self._Name(), ": " + cs_tools.red(warning_msg))
        elif not is_distributed and not self.parallel_type == "OpenMP":
            warning_msg  = 'WARNING: Global "parallel_type" (OpenMP) is different '
            warning_msg += 'from local one (' + self.parallel_type + ')!'
            cs_tools.solverprint(self.lvl, self._Name(), ": " + cs_tools.red(warning_msg))

    def _GetIOName(self):
        return "kratos"

    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()

    def _Name(self):
        return self.__class__.__name__

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = []
        # The list of processes will contain a list with each individual process already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: gravity is constructed first. Outlet process might need its information.
        # Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        if parameter_name == "processes":
            processes_block_names = ["gravity", "initial_conditions_process_list", "boundary_conditions_process_list", "auxiliar_process_list"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                info_msg  = "Using the old way to create the processes, this will be removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                KratosMultiphysics.Logger.PrintWarning("FluidDynamicsAnalysis", info_msg)
                from process_factory import KratosProcessFactory
                factory = KratosProcessFactory(self.model)
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        list_of_processes += factory.ConstructListOfProcesses(self.project_parameters[process_name])
            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        raise Exception("Mixing of process initialization is not alowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                info_msg  = "Using the old way to create the gid-output, this will be removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                KratosMultiphysics.Logger.PrintWarning("FluidDynamicsAnalysis", info_msg)

                if self.parallel_type == "OpenMP":
                    from gid_output_process import GiDOutputProcess as OutputProcess
                elif self.parallel_type == "MPI":
                    from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

                gid_output = OutputProcess(self.solver.GetComputingModelPart(),
                                       self.project_parameters["problem_data"]["problem_name"].GetString(),
                                       self.project_parameters["output_configuration"])

                list_of_processes += [gid_output,]
        else:
            raise NameError("wrong parameter name")

        return list_of_processes