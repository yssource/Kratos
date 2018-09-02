from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

available_ios = {
    "kratos" : "kratos_io",
    "sdof"   : "sdof_io",
    "mdof"   : "mdof_io"
}

def CreateIO(io_name, solvers, solver_name, level):
    """This function creates and returns the IO used for CoSimulation
    New IOs have to be registered by adding them to "available_ios"
    """

    import os
    from os import listdir
    from os.path import isfile, join
    onlyfiles = [print(f) for f in listdir(os.path.dirname(os.path.abspath(__file__))) if f.endswith("_io.py")]

    if io_name in available_ios:
        io_module = __import__(available_ios[io_name])
        return io_module.Create(solvers, solver_name, level)
    else:
        err_msg  = 'The requested IO "' + io_name + '" is not available!\n'
        err_msg += 'The following IOs are available:\n'
        for avail_io in available_ios:
            err_msg += "\t" + avail_io + "\n"
        raise NameError(err_msg)
