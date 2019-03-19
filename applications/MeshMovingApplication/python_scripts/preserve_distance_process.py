import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshMovingApplication as MeshMovingApplication

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return PreserveDistanceProcess(model, settings["Parameters"])



class PreserveDistanceProcess(KratosMultiphysics.Process):

    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "constraint_set_name"           : "LinearMasterSlaveConstraint",
            "master_sub_model_part_name"    : "master_connect",
            "slave_sub_model_part_name"     : "slave_connect",
            "variable_names"                : ["MESH_DISPLACEMENT_X","MESH_DISPLACEMENT_Y","MESH_DISPLACEMENT_Z"],
            "debug_info"                    : true,
            "must_find_neighbor"            : true,
            "neighbor_search_radius"        : 0.40,
            "bucket_size"                   : 10
        }
        """)
        #self.settings = settings.Clone()
        #settings.ValidateAndAssignDefaults(default_settings)
        default_settings.ValidateAndAssignDefaults(settings)
        self.settings = settings
        # The computing model part
        self.model = model
        #process_settings = settings["Parameters"]
        #self.model_part = model["FluidModelPart_MeshPart"]
        #model_part_name = self.settings["model_part_name"].GetString()
        #model_part_name = "Parts_Fluid"
        #computing_model_part = self.model.GetModelPart(model_part_name)



    def ExecuteBeforeSolutionLoop(self):

        self.computing_model_part = self.model["FluidModelPart_MeshPart"]
        self.master_model_part = self.model["FluidModelPart"].GetSubModelPart(self.settings["master_sub_model_part_name"].GetString())
        self.preserve_distance_process = MeshMovingApplication.PreserveDistanceProcess(self.model["FluidModelPart"], self.computing_model_part, self.settings)
        self.preserve_distance_process.ExecuteInitialize()
        #for constraint in self.master_model_part.MasterSlaveConstraints:
        #    self.computing_model_part.AddMasterSlaveConstraint(constraint)

    #def ExecuteFinalizeSolutionStep(self):
    #    self.preserve_distance_process.ExecuteFinalizeSolutionStep()
