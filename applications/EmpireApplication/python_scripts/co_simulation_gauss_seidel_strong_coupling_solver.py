from __future__ import print_function, absolute_import, division

# Importing the base class
from co_simulation_base_solver import CoSimulationBaseSolver

# Other imports
import co_simulation_convergence_accelerator_factory as convergence_accelerator_factory
import co_simulation_convergence_criteria_factory as convergence_criteria_factory
import co_simulation_tools as cosim_tools
from co_simulation_tools import csprint, red, green, cyan, bold, magenta

def CreateSolver(cosim_solver_settings, level):
    return GaussSeidelStrongCouplingSolver(cosim_solver_settings, level)

class GaussSeidelStrongCouplingSolver(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings, level):
        super(GaussSeidelStrongCouplingSolver, self).__init__(cosim_solver_settings, level)

        if not len(self.cosim_solver_settings["solvers"]) == 2:
            raise Exception("Exactly two solvers have to be specified for the GaussSeidelStrongCouplingSolver!")

        self.solver_names = []
        self.solvers = {}

        import python_solvers_wrapper_co_simulation as solvers_wrapper

        for solver_settings in self.cosim_solver_settings["coupling_loop"]:
            solver_name = solver_settings["name"]
            if solver_name in self.solver_names:
                raise NameError('Solver name "' + solver_name + '" defined twice!')
            self.solver_names.append(solver_name)
            self.cosim_solver_settings["solvers"][solver_name]["name"] = solver_name # adding the name such that the solver can identify itself
            self.solvers[solver_name] = solvers_wrapper.CreateSolver(
                self.cosim_solver_settings["solvers"][solver_name], self.lvl-1) # -1 to have solver prints on same lvl

        self.cosim_solver_details = cosim_tools.GetSolverCoSimulationDetails(
            self.cosim_solver_settings["coupling_loop"])

        # With this setting the coupling can start later
        self.start_coupling_time = 0.0
        if "start_coupling_time" in self.cosim_solver_settings:
            self.start_coupling_time = self.cosim_solver_settings["start_coupling_time"]
        self.coupling_start_printed = False

    def Initialize(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].Initialize()
        for solver_name in self.solver_names:
            self.solvers[solver_name].InitializeIO(self.solvers, self.cosim_solver_details)


        self.num_coupling_iterations = self.cosim_solver_settings["num_coupling_iterations"]
        self.convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(
            self.cosim_solver_settings["convergence_accelerator_settings"], self.solvers, self.cosim_solver_details, self.lvl)
        self.convergence_criteria = convergence_criteria_factory.CreateConvergenceCriteria(
            self.cosim_solver_settings["convergence_criteria_settings"], self.solvers, self.cosim_solver_details, self.lvl)

    def Finalize(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].Finalize()

    def AdvanceInTime(self, current_time):
        self.time = self.solvers[self.solver_names[0]].AdvanceInTime(current_time)
        for solver_name in self.solver_names[1:]:
            time_other_solver = self.solvers[solver_name].AdvanceInTime(current_time)
            if abs(self.time - time_other_solver) > 1e-12:
                raise Exception("Solver time mismatch")

        self.convergence_accelerator.AdvanceInTime()
        self.convergence_criteria.AdvanceInTime()

        if self.start_coupling_time > 0.0 and self.time > self.start_coupling_time and not self.coupling_start_printed:
            csprint(self.lvl, magenta("<< Starting Coupling >>"))
            self.coupling_start_printed = True

        return self.time

    def Predict(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].Predict()

    def InitializeSolutionStep(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].FinalizeSolutionStep()

    def OutputSolutionStep(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].OutputSolutionStep()

    def SolveSolutionStep(self):
        for k in range(self.num_coupling_iterations):
            csprint(self.lvl, cyan("Coupling iteration: ")+bold(str(k+1)+" / " + str(self.num_coupling_iterations)))
            for solver_name in self.solver_names:
                solver = self.solvers[solver_name]
                self.__SynchronizeInputData(solver, solver_name)
                solver.SolveSolutionStep()
                self.__SynchronizeOutputData(solver, solver_name)

            if self.convergence_criteria.IsConverged():
                csprint(self.lvl, green("##### CONVERGENCE AT INTERFACE WAS ACHIEVED #####"))
                break
            # else:
                ## ComputeUpdate(r, x)
                # @param r residual r_k
                # @param x solution x_k
                # Computes the approximated update in each iteration.
                def ComputeUpdate( self, r, x ):
                    return delta_x
                        #     self.convergence_accelerator.ComputeUpdate(...)
            if k+1 >= self.num_coupling_iterations:
                csprint(self.lvl, red("XXXXX CONVERGENCE AT INTERFACE WAS NOT ACHIEVED XXXXX"))

    def __SynchronizeInputData(self, solver, solver_name):
        if (self.time - self.start_coupling_time) > 1e-12:
            input_data_list = self.cosim_solver_details[solver_name]["input_data_list"]
            for input_data in input_data_list:
                from_solver = self.solvers[input_data["from_solver"]]
                solver.ImportData(input_data["data_name"], from_solver)

    def __SynchronizeOutputData(self, solver, solver_name):
        if (self.time - self.start_coupling_time) > 1e-12:
            output_data_list = self.cosim_solver_details[solver_name]["output_data_list"]
            for output_data in output_data_list:
                to_solver = self.solvers[output_data["to_solver"]]
                solver.ExportData(output_data["data_name"], to_solver)
