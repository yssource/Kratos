from __future__ import print_function, absolute_import, division

# Importing the base class
from co_simulation_base_solver import CoSimulationBaseSolver

# other imports
import co_sim_convergence_accelerator_factory as convergence_accelerator_factory
import co_sim_convergence_criteria_factory as convergence_criteria_factory
import io_factory
import co_simulation_tools as cosim_tools
from co_simulation_tools import csprint, red, green, cyan, bold

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
            self.solvers[solver_name] = solvers_wrapper.CreateSolver(
                self.cosim_solver_settings["solvers"][solver_name], self.lvl)

        self.solver_cosim_details = cosim_tools.GetSolverCoSimulationDetails(
            self.cosim_solver_settings["coupling_loop"])

    def Initialize(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].Initialize()

        self.num_coupling_iterations = self.cosim_solver_settings["num_coupling_iterations"]
        self.convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(
            self.cosim_solver_settings["convergence_accelerator_settings"])
        self.convergence_criteria = convergence_criteria_factory.CreateConvergenceCriteria(
            self.cosim_solver_settings["convergence_criteria_settings"], self.solvers, self.lvl)

    def Finalize(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].Finalize()

    def AdvanceInTime(self, current_time):
        new_time = self.solvers[self.solver_names[0]].AdvanceInTime(current_time)
        for solver_name in self.solver_names[1:]:
            new_time_2 = self.solvers[solver_name].AdvanceInTime(current_time)
            if abs(new_time- new_time_2) > 1e-12:
                raise Exception("Solver time mismatch")

        self.convergence_accelerator.AdvanceTimeStep()

        return new_time

    def Predict(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].Predict()

    def InitializeSolutionStep(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].InitializeSolutionStep()
        self.convergence_criteria.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].FinalizeSolutionStep()

    def SolveSolutionStep(self):
        for k in range(self.num_coupling_iterations):
            csprint(self.lvl, cyan("Coupling iteration: ")+bold(str(k+1)+" / " + str(self.num_coupling_iterations)))
            for solver_name in self.solver_names:
                solver = self.solvers[solver_name]
                # self.__SynchronizeInputData(solver, solver_name)
                solver.SolveSolutionStep()
                # self.__SynchronizeOutputData(solver, solver_name)

            ## TODO print coupling information here => then it is printed all the time! .... possible? => maybe print from convergence criteria ...?
            if self.convergence_criteria.IsConverged():
                csprint(self.lvl, green("##### CONVERGENCE AT INTERFACE WAS ACHIEVED #####"))
                break
            # else:
            #     self.convergence_accelerator.ComputeUpdate(...)
            if k+1 >= self.num_coupling_iterations:
                csprint(self.lvl, red("XXXXX CONVERGENCE AT INTERFACE WAS NOT ACHIEVED XXXXX"))



    def ImportData(self, DataName, FromClient):
        raise NotImplementedError("This needs to be implemented!")
    def ImportMesh(self, MeshName, FromClient):
        raise NotImplementedError("This needs to be implemented!")

    def ExportData(self, DataName, ToClient):
        raise NotImplementedError("This needs to be implemented!")
    def ExportMesh(self, MeshName, ToClient):
        raise NotImplementedError("This needs to be implemented!")

    def MakeDataAvailable(self, DataName, ToClient):
        raise NotImplementedError("This needs to be implemented!")
    def MakeMeshAvailable(self, MeshName, ToClient):
        raise NotImplementedError("This needs to be implemented!")

    def __SynchronizeInputData(self, solver, solver_name):
        input_data_list = self.solver_cosim_details[solver_name]["input_data_list"]
        for input_data in input_data_list:
            from_solver = self.solvers[input_data["from_solver"]]
            solver.ImportData(input_data["data_name"], from_solver)

    def __SynchronizeOutputData(self, solver, solver_name):
        output_data_list = self.solver_cosim_details[solver_name]["output_data_list"]
        for output_data in output_data_list:
            to_solver = self.solvers[output_data["to_solver"]]
            solver.ImportData(output_data["data_name"], to_solver)
