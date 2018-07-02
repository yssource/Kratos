from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

class CoSimulationAnalysis(object):
    """The base class for the CoSimulation-AnalysisStage
    It mimicks the AnalysisStage of Kratos but does NOT derive from it
    The reason is that it also needs to work without Kratos
    """
    def __init__(self, cosim_settings):
        if (type(cosim_settings) != dict):
            raise Exception("Input is expected to be provided as a python dictionary")

        self.cosim_settings = cosim_settings

    def Run(self):
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def RunSolutionLoop(self):
        while self.time < self.end_time:
            self.step += 1
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def Initialize(self):
        self._GetSolver().Initialize()

        ## Stepping and time settings
        self.end_time = self.cosim_settings["problem_data"]["end_time"]
        self.time = self.cosim_settings["problem_data"]["start_time"]
        self.step = 0

    def Finalize(self):
        self._GetSolver().Finalize()

    def InitializeSolutionStep(self):
        print("\x1b[1;34m", "\nCoSimulation:","\x1b[0m\x1b[1;1m", "time={0:.12g}".format(self.time), " | step="+ str(self.step), "\x1b[0m")

        self._GetSolver().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetSolver().FinalizeSolutionStep()

    def OutputSolutionStep(self):
        self._GetSolver().OutputSolutionStep()

    def _GetSolver(self):
        if not hasattr(self, '_solver'):
            self._solver = self._CreateSolver()
        return self._solver

    def _CreateSolver(self):
        """Create the solver
        """
        import python_solvers_wrapper_co_simulation as solvers_wrapper
        return solvers_wrapper.CreateSolver(self.cosim_settings["solver_settings"])

## TODO add the if name==main stuff like in fluid and structure