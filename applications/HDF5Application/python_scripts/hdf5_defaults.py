"""Defaults used in hdf file settings."""
import KratosMultiphysics

hdf5_default_settings = """
            {
                "model_part_name" : "MainModelPart",
                "file_settings" : {},
                "model_part_output_settings" : {},
                "nodal_results_settings" : {},
                "element_results_settings" : {},
                "output_time_settings" : {}
            }
            """

model_part_output_default_settings = """
            {
                "prefix" : "/ModelData"
            }
            """

temporal_default_settings = """
        {
            "prefix" : "/ResultsData",
            "list_of_variables": []
        }
        """

def GetPrefixes(filename):
    prefixes = {}
    parameter_file = open(filename,'r')
    ProjectParameters = KratosMultiphysics.Parameters( parameter_file.read())
    parameter_file.close()

    list_other_processes = ProjectParameters["list_other_processes"]
    for i in range(list_other_processes.size()):
        other_process = list_other_processes[i]
        
        if other_process["python_module"].GetString()!="single_mesh_temporal_output_process" and \
            other_process["python_module"].GetString()!="single_mesh_primal_output_process":
            continue
        
        settings = other_process["Parameters"]
        default_settings = KratosMultiphysics.Parameters(hdf5_default_settings)
        settings = settings.Clone()
        settings.ValidateAndAssignDefaults(default_settings)

        default_settings = KratosMultiphysics.Parameters(model_part_output_default_settings)
        model_part_settings = settings["model_part_output_settings"].Clone()
        model_part_settings.ValidateAndAssignDefaults(default_settings)

        if other_process["python_module"].GetString()=="single_mesh_temporal_output_process":
            nodal_settings = GetOutputPrefixes(settings["nodal_results_settings"])
        elif other_process["python_module"].GetString()=="single_mesh_primal_output_process":
            nodal_settings = GetOutputPrefixes(settings["nodal_results_settings"])
        

        if other_process["python_module"].GetString()=="single_mesh_temporal_output_process":
            element_settings = GetOutputPrefixes(settings["element_results_settings"])
        elif other_process["python_module"].GetString()=="single_mesh_primal_output_process":
            element_settings = GetOutputPrefixes(settings["element_results_settings"])     
        
    prefixes["hdf_file"] = settings["model_part_name"].GetString()
    prefixes["model_data"] = model_part_settings["prefix"].GetString()
    prefixes["nodal_data"] = nodal_settings["prefix"].GetString()
    prefixes["element_data"] = element_settings["prefix"].GetString()

    return prefixes

def GetOutputPrefixes(settings):
    default_settings = KratosMultiphysics.Parameters(temporal_default_settings)
    nodal_settings = settings.Clone()
    nodal_settings.ValidateAndAssignDefaults(default_settings)
    return nodal_settings