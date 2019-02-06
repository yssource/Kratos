import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return TurbulenceKEpsilonProcess(Model, settings["Parameters"])


class TurbulenceKEpsilonProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)

        # default_settings = KratosMultiphysics.Parameters("""
        # {
        #     "mesh_id"            : 0,
        #     "model_part_name"    : "",
        #     "variable_name"      : "PRESSURE",
        #     "constrained"        : true,
        #     "value"              : 0.0,
        #     "interval"           : [0.0,"End"],
        #     "hydrostatic_outlet" : false,
        #     "h_top"              : 0.0
        # }
        # """)

        # settings.ValidateAndAssignDefaults(default_settings)

    def ExecuteInitializeSolutionStep(self):
        print("ke: ExecuteInitializeSolutionStep")


    def ExecuteInitialize(self):
        print("ke: ExecuteInitialize")


    def Execute(self):
        print("ke: Execute")


    def ExecuteFinalizeSolutionStep(self):
        print("ke: ExecuteFinalizeSolutionStep")