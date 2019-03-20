# co simulation imports
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_coupled_solver import CoSimulationBaseCoupledSolver
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Other imports
import os
import time

def Create(custom_settings):
    return GaussSeidelIterativeStrongCouplingSolver(custom_settings)
#Comment Use "Strong" or "Iterative"
#Comment same for weak/loose
#Comment I vote for strong & weak
class GaussSeidelIterativeStrongCouplingSolver(CoSimulationBaseCoupledSolver):
    def __init__(self, custom_settings):
        super(GaussSeidelIterativeStrongCouplingSolver, self).__init__(custom_settings)
        if not self.number_of_participants == 2:
            raise Exception(cs_tools.bcolors.FAIL + "Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")
        self.num_coupling_iterations = self.settings["num_coupling_iterations"].GetInt()

        ### Importing the Participant modules
        self.coupling_sequence_list = self.full_settings["coupled_solver_settings"]["coupling_sequence"]

        ### Making the filters
        self._CreateFilters(self.coupling_sequence_list)

        ### Creating the convergence criterion
        self.convergence_criteria_list = self._CreateConvergenceCriteria(self.settings["convergence_criteria"])

        ### Creating the convergence accelerators
        self.convergence_accelerators_list = self._CreateConvergenceAccelerators(self.settings["convergence_accelerators"])

    def Initialize(self):
        super(GaussSeidelIterativeStrongCouplingSolver, self).Initialize()
        for conv_criteria in self.convergence_criteria_list:
            conv_criteria.Initialize()
        for accelerator in self.convergence_accelerators_list:
            accelerator.Initialize()

    def Finalize(self):
        super(GaussSeidelIterativeStrongCouplingSolver, self).Finalize()
        for conv_criteria in self.convergence_criteria_list:
            conv_criteria.Finalize()
        for accelerator in self.convergence_accelerators_list:
            accelerator.Finalize()


    def InitializeSolutionStep(self):
        super(GaussSeidelIterativeStrongCouplingSolver, self).InitializeSolutionStep()
        for conv_criteria in self.convergence_criteria_list:
            conv_criteria.InitializeSolutionStep()
        for accelerator in self.convergence_accelerators_list:
            accelerator.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(GaussSeidelIterativeStrongCouplingSolver, self).FinalizeSolutionStep()
        for conv_criteria in self.convergence_criteria_list:
            conv_criteria.FinalizeSolutionStep()
        for accelerator in self.convergence_accelerators_list:
            accelerator.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        if self.coupling_started:
            for iteration in range(self.num_coupling_iterations):
                #time.sleep(1)
                if self.echo_level > 0:
                    cs_tools.PrintInfo("\t"+ cs_tools.bcolors.HEADER + str(self._Name()) ,
                                        cs_tools.bcolors.MAGENTA + "Coupling iteration: ", cs_tools.bcolors.BOLD + str(iteration+1) +
                                        " / " + cs_tools.bcolors.BLUE + str(self.num_coupling_iterations) + cs_tools.bcolors.ENDC)

                for conv_criteria in self.convergence_criteria_list:
                    conv_criteria.InitializeCouplingIteration()
                for accelerator in self.convergence_accelerators_list:
                    accelerator.InitializeCouplingIteration()

                for solver_name, solver in self.participating_solvers.items():
                    solver.InitializeCouplingIteration()

                for solver_name, solver in self.participating_solvers.items():
                    self._SynchronizeInputData(solver_name)
                    cs_tools.PrintInfo("\t"+cs_tools.bcolors.GREEN + cs_tools.bcolors.BOLD + "SolveSolutionStep for Solver", solver_name + cs_tools.bcolors.ENDC)
                    solver.SolveSolutionStep()
                    self._SynchronizeOutputData(solver_name)

                for solver_name, solver in self.participating_solvers.items():
                    solver.FinalizeCouplingIteration()

                for accelerator in self.convergence_accelerators_list:
                    accelerator.FinalizeCouplingIteration()
                for conv_criteria in self.convergence_criteria_list:
                    conv_criteria.FinalizeCouplingIteration()

                is_converged = True #Comment I think this would be suitable for list-comprehension
                for conv_criteria in self.convergence_criteria_list:
                    is_converged = is_converged and conv_criteria.IsConverged()
                if is_converged or iteration+1 >= self.num_coupling_iterations:
                    if self.echo_level > 0:
                        if is_converged:
                            cs_tools.PrintInfo(cs_tools.bcolors.GREEN + "### CONVERGENCE WAS ACHIEVED ###" + cs_tools.bcolors.ENDC )
                        if iteration+1 >= self.num_coupling_iterations:
                            cs_tools.PrintWarning("\t"+cs_tools.bcolors.FAIL + "### CONVERGENCE NOT ACHIEVED IN STRONG COUPLING ITERATIONS ###" + cs_tools.bcolors.ENDC)
                    break
                else:
                    for accelerator in self.convergence_accelerators_list:
                        accelerator.ComputeUpdate()
        else:
            for solver_name, solver in self.participating_solvers.items():
                cs_tools.PrintInfo("\t"+cs_tools.bcolors.GREEN + cs_tools.bcolors.BOLD + "SolveSolutionStep for Solver", solver_name + cs_tools.bcolors.ENDC)
                solver.SolveSolutionStep()

    def _Name(self):
        return self.settings['name'].GetString()


    def PrintInfo(self):
        super(GaussSeidelIterativeStrongCouplingSolver, self).PrintInfo()

    ## _CreateFilters : Protected Function to make filter objects list and store in the datafield
    #
    #  @param conv_acc_settings dict: setting of the convergence accelerator to be make
    def _CreateFilters(self, co_simulation_solver_settings): # probably better in some utils file
        import KratosMultiphysics.CoSimulationApplication.co_simulation_convergence_accelerators.co_simulation_convergence_accelerator_factory as factory
        num_solvers = co_simulation_solver_settings.size()
        solver_cosim_details = {}
        for i_solver in range(num_solvers):
            solver_name = co_simulation_solver_settings[i_solver]["name"].GetString()
            solver = self.participating_solvers[solver_name]

            ## for all the input data
            input_data_list = self.coupling_sequence[solver_name]["input_data_list"]
            num_input_data = input_data_list.size()
            for i in range(num_input_data):
                input_data = input_data_list[i]
                filters_list = input_data["filters"]
                destination_data = solver.GetInterfaceData(input_data["destination_data"].GetString())
                for filter in filters_list:
                    accelerator = factory.CreateConvergenceAccelerator(filter, destination_data) ## TODO: should change to filter
                    destination_data.filters.append(accelerator)


            ## for all the output data
            output_data_list = self.coupling_sequence[solver_name]["output_data_list"]
            num_output_data = output_data_list.size()
            for i in range(num_output_data):
                output_data = output_data_list[i]
                filters_list = output_data["filters"]
                origin_data = solver.GetInterfaceData(output_data["origin_data"].GetString())
                for filter in filters_list:
                    accelerator = factory.CreateConvergenceAccelerator(filter, destination_data) ## TODO: should change to filter
                    origin_data.filters.append(accelerator)

    ## _CreateConvergenceCriteria : Private Function to make convergence criteria objects list #Comment protected
    #
    #  @param conv_acc_settings dict: setting of the convergence criteria to be make
    def _CreateConvergenceCriteria(self, conv_criteria_settings): # probably better in some utils file
        conv_criteria = []
        import KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_convergence_criteria as criteria
        num_criteria = conv_criteria_settings.size()
        for i in range(num_criteria):
            criteria_setting = conv_criteria_settings[i]
            solver_name = criteria_setting["solver"].GetString()
            solver = self.participating_solvers[solver_name]
            data_name = criteria_setting["data_name"].GetString()
            data = solver.data_map[data_name]
            criteria = criteria.Create(criteria_setting, data) # Change to use interface data
            conv_criteria.append(criteria)

        return conv_criteria

    ## _CreateConvergenceAccelerators : Protected Function to make convergence accelerator objects list
    #
    #  @param conv_acc_settings dict: setting of the convergence accelerator to be make
    def _CreateConvergenceAccelerators(self, conv_acc_settings): # probably better in some utils file
        conv_accelerators = []
        num_acceleratos = conv_acc_settings.size()
        import KratosMultiphysics.CoSimulationApplication.co_simulation_convergence_accelerators.co_simulation_convergence_accelerator_factory as factory
        for i in range(num_acceleratos):
            accelerator_settings = conv_acc_settings[i]
            solver_name = accelerator_settings["solver"].GetString()
            solver = self.participating_solvers[solver_name]
            data = solver.GetInterfaceData(accelerator_settings["data_name"].GetString())
            accelerator = factory.CreateConvergenceAccelerator(accelerator_settings, data)
            conv_accelerators.append(accelerator)

        return conv_accelerators
