import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import math

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFarFieldProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyFarFieldProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "inlet_phi": 1.0,
                "velocity_infinity": [3.4,0.0,0],
                "density_infinity"  : 1.0,
                "mach_infinity": 0.01,
                "gamma": 1.4,
                "pressure_infinity": 101325
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.fluid_model_part = self.model_part.GetRootModelPart()

        self.inlet_phi = settings["inlet_phi"].GetDouble()
        self.velocity_inf = KratosMultiphysics.Vector(3)
        self.velocity_inf[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_inf[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_inf[2] = settings["velocity_infinity"][2].GetDouble()
        self.density_infinity = settings["density_infinity"].GetDouble()
        self.mach_infinity = settings["mach_infinity"].GetDouble()
        self.gamma = settings["gamma"].GetDouble()
        self.pressure_infinity = settings["pressure_infinity"].GetDouble()

        self.u_infinity = math.sqrt(
            self.velocity_inf[0]**2 + self.velocity_inf[1]**2 + self.velocity_inf[2]**2)
        self.a_infinity = self.u_infinity / self.mach_infinity

        # For the model part
        self.model_part.ProcessInfo.SetValue(CPFApp.VELOCITY_INFINITY, self.velocity_inf)

        # For the conditions
        self.fluid_model_part.GetProperties()[0].SetValue(CPFApp.DENSITY_INFINITY, self.density_infinity)

        # For the elements
        self.fluid_model_part.GetProperties()[1].SetValue(CPFApp.VELOCITY_INFINITY, self.velocity_inf)
        self.fluid_model_part.GetProperties()[1].SetValue(CPFApp.DENSITY_INFINITY, self.density_infinity)
        self.fluid_model_part.GetProperties()[1].SetValue(CPFApp.MACH_INFINITY, self.mach_infinity)
        self.fluid_model_part.GetProperties()[1].SetValue(CPFApp.GAMMA, self.gamma)
        self.fluid_model_part.GetProperties()[1].SetValue(KratosMultiphysics.SOUND_VELOCITY, self.a_infinity)
        self.fluid_model_part.GetProperties()[1].SetValue(CPFApp.PRESSURE_INFINITY, self.pressure_infinity)

    def Execute(self):
        #KratosMultiphysics.VariableUtils().SetVectorVar(CPFApp.VELOCITY_INFINITY, self.velocity_inf, self.model_part.Conditions)
        for cond in self.model_part.Conditions:
            cond.SetValue(CPFApp.VELOCITY_INFINITY, self.velocity_inf)

        # Select the first node in the mesh as reference
        for node in self.model_part.Nodes:
            x0 = node.X
            y0 = node.Y
            z0 = node.Z
            break

        # Find smallest distance_to_reference
        pos = 1e30
        for node in self.model_part.Nodes:
            # Computing distance to reference
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0

            distance_to_reference = dx*self.velocity_inf[0] + dy*self.velocity_inf[1] + dz*self.velocity_inf[2]

            if(distance_to_reference < pos):
                pos = distance_to_reference
                self.inlet_node = node

        # Fix nodes in the inlet
        for node in self.model_part.Nodes:
            # Computing distance to reference
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0

            distance_to_reference = dx*self.velocity_inf[0] + dy*self.velocity_inf[1] + dz*self.velocity_inf[2]

            if(distance_to_reference < pos+1e-9):
                node.Fix(CPFApp.VELOCITY_POTENTIAL)
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,0,self.inlet_phi)
                if self.model_part.HasNodalSolutionStepVariable(CPFApp.ADJOINT_VELOCITY_POTENTIAL):
                    node.Fix(CPFApp.ADJOINT_VELOCITY_POTENTIAL)
                    node.SetSolutionStepValue(CPFApp.ADJOINT_VELOCITY_POTENTIAL,0,0.0)

        #TODO: Check best way of initializing the flow field
        '''
        # Compute free stream direction
        self.free_stream_direction = self.velocity_inf
        self.free_stream_direction[0] /= self.u_infinity
        self.free_stream_direction[1] /= self.u_infinity
        self.free_stream_direction[2] /= self.u_infinity
        
        for node in self.fluid_model_part.Nodes:
            # Computing distance to inlet
            dx = node.X - self.inlet_node.X
            dy = node.Y - self.inlet_node.Y
            dz = node.Z - self.inlet_node.Z

            distance_to_inlet = dx*self.free_stream_direction[0] + dy*self.free_stream_direction[1] + dz*self.free_stream_direction[2]

            initial_potential = self.inlet_phi + self.u_infinity * distance_to_inlet
            #node.SetValue(KratosMultiphysics.WATER_PRESSURE,initial_potential)
            #node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,0,initial_potential)
            #node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, 0, initial_potential)
        '''

    def ExecuteInitializeSolutionStep(self):
        self.Execute()
