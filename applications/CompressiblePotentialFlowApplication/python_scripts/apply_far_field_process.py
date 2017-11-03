import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFarFieldProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ApplyFarFieldProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "main_model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "inlet_phi": 1.0,
                "velocity_infinity": [1.0,0.0,0],
                "variable_name":   "POSITIVE_FACE_PRESSURE"
            }  """ );
        
            
        settings.ValidateAndAssignDefaults(default_parameters);
        
        self.main_model_part = Model[settings["main_model_part_name"].GetString()]
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.velocity_infinity = KratosMultiphysics.Vector(3)#array('d', [1.0, 2.0, 3.14])#np.array([0,0,0])#np.zeros(3)#vector(3)
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        #self.density_infinity = settings["density_infinity"].GetDouble() #TODO: must read this from the properties
        self.inlet_phi = settings["inlet_phi"].GetDouble()
        self.variable = getattr(KratosMultiphysics, settings["variable_name"].GetString())
        
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.VELOCITY,self.velocity_infinity)
        
        self.zero_vector = KratosMultiphysics.Vector(3)
        self.zero_vector[0] = 0.0
        self.zero_vector[1] = 0.0
        self.zero_vector[2] = 0.0
        
        for cond in self.main_model_part.Conditions:
            cond.SetValue(KratosMultiphysics.VELOCITY, self.zero_vector)
        
        
    def Execute(self):
        for cond in self.model_part.Conditions:
            cond.SetValue(KratosMultiphysics.VELOCITY, self.velocity_infinity)
        

        #select the first node
        for node in self.model_part.Nodes:
            node1 = node
            break
        
        #find the node with the minimal x
        x0 = node1.X
        y0 = node1.X
        z0 = node1.X
        
        pos = 1e30
        for node in self.model_part.Nodes:
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0
            
            tmp = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]
            
            if(tmp < pos):
                pos = tmp

        for node in self.model_part.Nodes:
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0
            
            tmp = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]
            
            if(tmp < pos+1e-9):
                node.Fix(self.variable)
                node.SetSolutionStepValue(self.variable,0,self.inlet_phi)
                #print(node)
        
    def ExecuteInitializeSolutionStep(self):
        self.Execute()
        
        