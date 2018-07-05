from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

available_ios = {
    "dummy"         : "co_sim_base_io",
    "kratos_field"  : "kratos_field_io",
    "kratos_signal" : "kratos_signal_io"
}

def CreateIO(io_settings, solvers, solver_name, cosim_solver_details, level):
    """This function creates and returns the IO used for CoSimulation
    New IOs have to be registered by adding them to "available_ios"
    """
    if (type(io_settings) != dict):
        raise Exception("Input is expected to be provided as a python dictionary")

    io_type = io_settings["io_type"]

    if io_type in available_ios:
        io_module = __import__(available_ios[io_type])
        return io_module.Create(io_settings, solvers, solver_name, cosim_solver_details, level)
    else:
        err_msg  = 'The requested IO "' + io_type + '" is not available!\n'
        err_msg += 'The following IOs are available:\n'
        for avail_io in available_ios:
            err_msg += "\t" + avail_io + "\n"
        raise NameError(err_msg)
