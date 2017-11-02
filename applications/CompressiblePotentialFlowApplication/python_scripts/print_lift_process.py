import KratosMultiphysics
import KratosMultiphysics.AdjointFluidApplication as AdjointFluidApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PrintLiftProcess(Model, settings["Parameters"])

class PrintLiftProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        # default settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
           "response_function_settings" : {
                "response_type"            : "lift",
                "sensitivity_model_part_name"  : "Boundary",
                "nodal_sensitivity_variables" : ["SHAPE_SENSITIVITY"],
                "custom_settings" : {
                    "structure_model_part_name" : "Body2D_Body_Auto1",
                    "lift_direction"            : [0.0, 1.0, 0.0]
                }                
            },
            "lift_model_part_name" : "MainModelPart",
            "lift_file_name"       : "test_ps_sensitivity_2d/one_element_test.dat"
        }""")
        
        settings.ValidateAndAssignDefaults(default_parameters);
        
        self.lift_model_part = Model[settings["lift_model_part_name"].GetString()]
        self.lift_file_name = settings["lift_file_name"].GetString()
        
        #self.response_function = AdjointFluidApplication.PotentialFlowLiftResponseFunction2D(self.model_part, self.settings["response_function_settings"])
        
        self.response_function = AdjointFluidApplication.PotentialFlowLiftResponseFunction2D(self.lift_model_part, default_parameters["response_function_settings"])

    def ExecuteInitialize(self):
        self.lift_file = open(self.lift_file_name, 'w')
        self.lift_file.write("#lift coefficient cl\n")

    def ExecuteFinalizeSolutionStep(self):
        liftc = self.response_function.CalculateValue(self.lift_model_part)
        
        self.lift_file.write('{0:22.15e}\n'.format(liftc))
        self.lift_file.flush()

    def ExecuteFinalize(self):
        self.lift_file.close()
