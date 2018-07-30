from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

CheckRegisteredApplications("FluidDynamicsApplication")

def Factory(process_settings, model):
    if(type(process_settings) != Parameters):
        raise Exception("Unexpected type for \"process_settings\" argument. Expected input is a Kratos::Parameters object.")

    default_settings = Parameters(r'''{
        "model_part_name" : "",
        "start_time" : 0.05,
        "reset_containers" : true
    }''')

    settings = process_settings["Parameters"]
    settings.ValidateAndAssignDefaults(default_settings)

    model_part_name = settings["model_part_name"].GetString()
    if model_part_name == "":
        raise Exception("Please define the name of the model part (set the \"model_part_name\" (string) argument.")
    model_part = model.GetModelPart(model_part_name)

    start_time = settings["start_time"].GetDouble()
    reset_containers = settings["reset_containers"].GetBool()

    return TurbulenceStatisticsProcess(model_part,start_time,reset_containers)
