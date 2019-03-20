from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

def CreateConvergenceAccelerator(convergence_accelerator_settings, data):
    """This function creates and returns the Convergence Accelerator used for CoSimulation
    The convergence-accelerator-module has to be on the PYTHONPATH
    Naming-Convention: The module-file has to end with "_convergence_accelerator"
    """
    folder_path = "KratosMultiphysics.CoSimulationApplication.co_simulation_convergence_accelerators."
    conv_acc_module_name = convergence_accelerator_settings["type"].GetString() + "_convergence_accelerator"
    total_path = folder_path+conv_acc_module_name
    accelerator_module = __import__(total_path,fromlist=[conv_acc_module_name])

    return accelerator_module.Create(convergence_accelerator_settings, data)
