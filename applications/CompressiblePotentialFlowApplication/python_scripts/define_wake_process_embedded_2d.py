import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import KratosMultiphysics.MeshingApplication as MeshingApplication
from KratosMultiphysics.CompressiblePotentialFlowApplication.define_wake_process_2d import DefineWakeProcess2D
import math

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return DefineWakeProcessEmbedded(Model, settings["Parameters"])

class DefineWakeProcessEmbedded(DefineWakeProcess2D):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"                   : 0,
                "model_part_name"           : "MainModelPart",
                "wake_direction"                 : [1.0,0.0,0.0],
                "epsilon"    : 1e-9
            }
            """)

        settings.ValidateAndAssignDefaults(default_settings)
        # TODO Implement this process in C++ and make it open mp parallel to save time selecting the wake elements

        self.wake_direction = settings["wake_direction"].GetVector()
        if(self.wake_direction.Size() != 3):
            raise Exception('The wake direction should be a vector with 3 components!')

        dnorm = math.sqrt(
            self.wake_direction[0]**2 + self.wake_direction[1]**2 + self.wake_direction[2]**2)
        self.wake_direction[0] /= dnorm
        self.wake_direction[1] /= dnorm
        self.wake_direction[2] /= dnorm

        self.wake_normal = KratosMultiphysics.Vector(3)
        self.wake_normal[0] = -self.wake_direction[1]
        self.wake_normal[1] = self.wake_direction[0]
        self.wake_normal[2] = 0.0

        self.epsilon = settings["epsilon"].GetDouble()

        self.fluid_model_part = Model[settings["model_part_name"].GetString()].GetRootModelPart()
        self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart(
            "trailing_edge_model_part")

    def ExecuteInitialize(self):
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.fluid_model_part,
                                                        self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
        # Neigbour search tool instance
        AvgElemNum = 10
        AvgNodeNum = 10
        nodal_neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(
            self.fluid_model_part, AvgElemNum, AvgNodeNum)
        # Find neighbours
        nodal_neighbour_search.Execute()
        super(DefineWakeProcessEmbedded, self).ExecuteInitialize()
    
    def SaveTrailingEdgeNode(self):
        # This function finds and saves the trailing edge for further computations
        kutta_model_part=self.fluid_model_part.GetSubModelPart('KuttaLS')
        for node in kutta_model_part.Nodes:
            self.te = node

