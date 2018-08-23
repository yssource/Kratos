def GetSolverCoSimulationDetails(co_simulation_solver_settings):
    num_solvers = len(co_simulation_solver_settings)
    solver_cosim_details = {}
    for solver_settings in co_simulation_solver_settings:
        solver_name = solver_settings["name"]
        solver_cosim_details[solver_name] = solver_settings
    # TODO check if the data is consitently defined! => maybe do at another place though...
    # - input in one is output in another
    # - one IO is defined for each data_name
    # - if the same data is defined multiple times
    # - check if data format has been specified
    return solver_cosim_details

def ImportArrayFromSolver(solver, data_name, data_array, buffer_index=0):
    data_settings = {
        "data_format"  : "numpy_array",
        "data_name"    : data_name,
        "data_array"   : data_array,
        "buffer_index" : buffer_index
    }

    solver.ExportData(data_settings, solver)

def ExportArrayToSolver(solver, data_name, data_array, buffer_index=0):
    data_settings = {
        "data_format"  : "numpy_array",
        "data_name"    : data_name,
        "data_array"   : data_array,
        "buffer_index" : buffer_index
    }

    solver.ImportData(data_settings, solver)



class CoSimulationParameters(object):
    """Object to mimic Kratos::Parameters
    It works with a python-dictionary, thus it does not
    need many of the auxilliary methods
    """
    def __init__(self):
        pass

    def ValidateAndAssignDefaults(self, defaults):
        pass

    def RecursivelyValidateAndAssignDefaults(self, defaults):
        pass

    def Merge(self):
        """TODO: not really sure how to implement this, have to take a look at
        the issue again
        """
        pass


class CoSimulationSpace(object):
    pass

    def IsDistributed(self):
        return False

    def Barrier(self):
        pass

    def Rank(self):
        return 0

    def Size(self):
        return 1

# class CoSimulationMPISpace(object):
#     """This requires some MPI-commands exposed to Python
#     This is currently only available with Kratos,
#     i.e. MPI can only be used with Kratos compiled with MPI
#     """
#     try:
#         import KratosMultiphysics
#         from KratosMultiphysics.mpi import mpi
#     except ImportError:
#         raise Exception("Running in MPI currently requires Kratos-MPI!")

#     def __init__(self):

#         # Precompute rank and size such that they don't have to be recomputed all the time
#         self.comm_rank = self.mpi.rank
#         self.comm_size = self.mpi.size

#         if self.comm_size < 2:
#             raise Exception("Running in MPI requires at least 2 processes!")

#     def IsDistributed(self):
#         return True

#     def Barrier(self):
#         self.mpi.world.barrier()

#     def Rank(self):
#         return self.comm_rank

#     def Size(self):
#         return self.comm_size



# Global Object for wrapping calls that are different in OpenMP/MPI
COSIM_SPACE = CoSimulationSpace() # I think this works, could be overridden by the CoSimAnalysis





PRINT_COLORS = False # Global var to specify if colors should be printed
PRINTING_RANK = True # Global var to specify if this rank is the printing one in MPI

def color_string(string2color, color_code):
    if PRINT_COLORS:
        return "\x1b["+color_code+"m" + str(string2color) + "\x1b[0m"
    else:
        return string2color

def bold(string2color):
    return color_string(string2color, "1;1")
def italic(string2color):
    return color_string(string2color, "1;3")
def darkify(string2color):
    return bold(color_string(string2color, "1;2")) # bold is needed bcs it is removed otherwise
def underline(string2color):
    return color_string(string2color, "1;4")

def blue(string2color):
    return color_string(string2color, "1;34")
def darkblue(string2color):
    return (darkify(blue(string2color)))

def red(string2color):
    return color_string(string2color, "1;31")
def darkred(string2color):
    return (darkify(red(string2color)))

def green(string2color):
    return color_string(string2color, "1;32")
def darkgreen(string2color):
    return (darkify(green(string2color)))

def yellow(string2color):
    return color_string(string2color, "1;33")
def darkyellow(string2color):
    return (darkify(yellow(string2color)))

def cyan(string2color):
    return color_string(string2color, "1;36")
def darkcyan(string2color):
    return (darkify(cyan(string2color)))

def magenta(string2color):
    return color_string(string2color, "1;35")
def darkmagenta(string2color):
    return (darkify(magenta(string2color)))

SPACE = 4 * " "
def csprint(level, *args):
    if PRINTING_RANK:
        print(level*SPACE + blue("<CS-"+str(level)+">"), " ".join(map(str,args)))
    COSIM_SPACE.Barrier()

def solverprint(level, solver_name, *args):
    csprint(level, yellow(solver_name + ":"), *args)

def couplingsolverprint(level, solver_name, *args):
    csprint(level, darkyellow(solver_name + ":"), *args)

def classprint(level, solver_name, *args):
    csprint(level, magenta(solver_name + ":"), *args)

if __name__ == "__main__":
    print("printing all color options:\n")

    str2print = "MyCustomString"

    PRINT_COLORS = True

    print("print:", str2print)

    print("bold:", bold(str2print))
    print("italic:", italic(str2print))
    print("darkify:", darkify(str2print))
    print("underline:", underline(str2print))

    print("blue:", blue(str2print))
    print("darkblue:", darkblue(str2print))

    print("red:", red(str2print))
    print("darkred:", darkred(str2print))

    print("green:", green(str2print))
    print("darkgreen:", darkgreen(str2print))

    print("yellow:", yellow(str2print))
    print("darkyellow:", darkyellow(str2print))

    print("cyan:", cyan(str2print))
    print("darkcyan:", darkcyan(str2print))

    print("magenta:", magenta(str2print))
    print("darkmagenta:", darkmagenta(str2print))