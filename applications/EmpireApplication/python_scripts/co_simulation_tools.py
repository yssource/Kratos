import json

def GetSolverCoSimulationDetails(co_simulation_solver_settings):
    num_solvers = len(co_simulation_solver_settings)
    solver_cosim_details = {}
    for solver_settings in co_simulation_solver_settings:
        solver_name = solver_settings["name"]
        solver_cosim_details[solver_name] = solver_settings
    return solver_cosim_details

def CheckCoSimulationSettings(co_simulation_settings):
    # TODO check if the data is consitently defined! => maybe do at another place though...
    # - input in one is output in another
    # - one IO is defined for each data_name
    # - if the same data is defined multiple times
    # - check if data format has been specified

    if not "problem_data" in co_simulation_settings:
        raise Exception('"problem_data" not existing in settings!')
    if not "solver_settings" in co_simulation_settings:
        raise Exception('"solver_settings" not existing in settings!')

    problem_data = co_simulation_settings["problem_data"]
    if not "start_time" in problem_data:
        raise Exception('"start_time" not defined in "problem_data"!')
    start_time = problem_data["start_time"]
    if not type(start_time) is float:
        raise Exception('"start_time" has to be defined as a float!')
    if start_time < 0.0:
        raise Exception('"start_time" has to be >= 0.0!')

    if not "end_time" in problem_data:
        raise Exception('"end_time" not defined in "problem_data"!')
    end_time = problem_data["end_time"]
    if not type(start_time) is float:
        raise Exception('"end_time" has to be defined as a float!')
    if start_time < 0.0:
        raise Exception('"end_time" has to be >= 0.0!')
    if start_time > end_time:
        raise Exception('"end_time" has to be > "start_time"!')




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


def ValidateAndAssignDefaults(defaults, settings, recursive=False):
    for key, val in settings.items():
        # check if the current entry also exists in the defaults
        if not key in defaults.keys():
            err_msg  = 'The item with name "' + key + '" is present in this '
            err_msg += 'settings\nbut NOT in the defaults!\n'
            err_msg += 'settings are:\n'
            err_msg += json.dumps(settings, indent=4)
            err_msg += '\ndefaults are:\n'
            err_msg += json.dumps(defaults, indent=4)
            raise Exception(err_msg)

        # check if the type is the same in the defaults
        if type(settings[key]) != type(defaults[key]):
            err_msg  = 'The type of the item with name "' + key + '" (type: "'
            err_msg += str(type(settings[key]).__name__)+'") in this '
            err_msg += 'settings\nis NOT the same as in the defaults (type: "'
            err_msg += str(type(defaults[key]).__name__)+'")!\n'
            err_msg += 'settings are:\n'
            err_msg += json.dumps(settings, indent=4)
            err_msg += '\ndefaults are:\n'
            err_msg += json.dumps(defaults, indent=4)
            raise Exception(err_msg)

    # loop the defaults and add the missing entries
    for key_d, val_d in defaults.items():
        if key_d not in settings: # add the default in case the setting is not present
            settings[key_d] = val_d
        elif recursive and type(val_d) is dict:
            RecursivelyValidateAndAssignDefaults(val_d, val)

def RecursivelyValidateAndAssignDefaults(defaults, settings):
    ValidateAndAssignDefaults(defaults, settings, recursive=True)


class CoSimulationSpace(object):

    mpi_op_types = ["min", "max", "sum"]

    def IsDistributed(self):
        return False

    def Barrier(self):
        pass

    def Rank(self):
        return 0

    def Size(self):
        return 1

    def BroadCast(self, val_to_broadcast, rank_to_broadcast_from):
        self._CheckInputFor_BroadCast(val_to_broadcast, rank_to_broadcast_from)
        return val_to_broadcast

    def Gather(self, my_value, rank_to_gather_on):
        self._CheckInputFor_Gather(my_value, rank_to_gather_on)
        return [my_value] # the mpi-call also returns a list

    def AllGather(self, my_value):
        self._CheckInputFor_AllGather(my_value)
        return [my_value] # the mpi-call also returns a list

    def GatherV(self, my_values, rank_to_gather_on):
        self._CheckInputFor_GatherV(my_values, rank_to_gather_on)
        return [my_value] # the mpi-call also returns a list of lists

    def Scatter(self, vals_to_scatter, rank_to_scatter_from):
        self._CheckInputFor_Scatter(vals_to_scatter, rank_to_scatter_from)
        return vals_to_scatter[self.Rank]

    def ScatterV(self, vals_to_scatter, rank_to_scatter_from):
        self._CheckInputFor_ScatterV(vals_to_scatter, rank_to_scatter_from)
        return vals_to_scatter[self.Rank]

    def Reduce(self, my_value, rank_to_reduce_on, mpi_op):
        self._CheckInputFor_Reduce(my_value, rank_to_reduce_on, mpi_op)
        return my_value

    def AllReduce(self, my_value, rank_to_reduce_on, mpi_op):
        self._CheckInputFor_AllReduce(my_value, rank_to_reduce_on, mpi_op)
        return my_value

    def _CheckInputFor_BroadCast(val_to_broadcast, rank_to_broadcast_from):
        if not type(val_to_broadcast) is int or not type(val_to_broadcast) is float:
            raise Exception('Acceptable types are "int" and "float"')

    def _CheckInputFor_Gather(my_value, rank_to_gather_on):
        if not type(val_to_broadcast) is int or not type(val_to_broadcast) is float:
            raise Exception('Acceptable types are "int" and "float"')

    def _CheckInputFor_AllGather(my_value):
        if not type(val_to_broadcast) is int or not type(val_to_broadcast) is float:
            raise Exception('Acceptable types are "int" and "float"')

    def _CheckInputFor_GatherV(my_values, rank_to_gather_on):
        if not type(my_values) is list:
            raise Exception('Acceptable types are "list"')

    def _CheckInputFor_Scatter(vals_to_scatter, rank_to_scatter_from):
        if not type(my_values) is list:
            raise Exception('Acceptable types are "list"')
        if len(vals_to_scatter) != self.Size():
            raise Exception('The number of values to scatter is different from Size!')
        if not type(vals_to_scatter[self.Rank]) is list:
            raise Exception('Acceptable types are "list of lists"')

    def _CheckInputFor_ScatterV(vals_to_scatter, rank_to_scatter_from):
        if not type(my_values) is list:
            raise Exception('Acceptable types are "list"')
        if len(vals_to_scatter) != self.Size():
            raise Exception('The number of values to scatter is different from Size!')

    def _CheckInputFor_Reduce(my_value, rank_to_reduce_on, mpi_op):
        if not type(val_to_broadcast) is int or not type(val_to_broadcast) is float:
            raise Exception('Acceptable types are "int" and "float"')
        if not mpi_op in self.mpi_op_types:
            raise Exception('Acceptable operations are ' + str(self.mpi_op_types))

    def _CheckInputFor_AllReduce(my_value, rank_to_reduce_on, mpi_op):
        if not type(val_to_broadcast) is int or not type(val_to_broadcast) is float:
            raise Exception('Acceptable types are "int" and "float"')
        if not mpi_op in self.mpi_op_types:
            raise Exception('Acceptable operations are ' + str(self.mpi_op_types))


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