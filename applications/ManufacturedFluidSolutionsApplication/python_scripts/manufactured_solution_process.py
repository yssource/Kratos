# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ManufacturedFluidSolutionsApplication as MS

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ManufacturedSolutionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ManufacturedSolutionProcess(KM.Process):
    '''
    A wrapper for all the processes related to the manufactured solution,
    such as the body force, the boundary conditions o the error computation processes
    '''
    def __init__(self, model, settings ):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"          : "model_part",
                "manufactured_name"        : "manufactured_solution_name",
                "manufactured_parameters"  : {},
                "processes_list"           : []
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        # Build the manufactured solution
        fluid_properties = model[settings["model_part_name"].GetString()].ElementsArray(0)[0].Properties
        manufactured_class = getattr(MS, settings["manufactured_name"].GetString())
        self.manufactured = manufactured_class(fluid_properties, settings["manufactured_parameters"])

        # Build the sub processes
        from process_factory import KratosProcessFactory
        factory = KratosProcessFactory(model)
        self.processes = factory.ConstructListOfProcesses(settings["processes_list"])

        for process in self.processes:
            process.SetManufacturedSolution(self.manufactured)

    def Check(self):
        for process in self.processes:
            process.Check()

    def ExecuteInitialize(self):
        for process in self.processes:
            process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        for process in self.processes:
            process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        for process in self.processes:
            process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        for process in self.processes:
            process.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        for process in self.processes:
            process.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        for process in self.processes:
            process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        for process in self.processes:
            process.ExecuteFinalize()

    def Clear(self):
        for process in self.processes:
            process.Clear()
