def GetSolverCoSimulationDetails(co_simulation_solver_settings):
    num_solvers = len(co_simulation_solver_settings)
    solver_cosim_details = {}
    for solver_settings in co_simulation_solver_settings:
        solver_name = solver_settings["name"]
        solver_cosim_details[solver_name] = solver_settings

    return solver_cosim_details

PRINT_COLORS = False # Global var to specify if colors should be printed

def color_string(string2color, color_code):
    if PRINT_COLORS:
        return "\x1b["+color_code+"m" + str(string2color) + "\x1b[0m"
    else:
        return string2color

def blue(string2color):
    return color_string(string2color, "1;34")
def red(string2color):
    return color_string(string2color, "1;31")
def green(string2color):
    return color_string(string2color, "1;32")
def yellow(string2color):
    return color_string(string2color, "1;33")
def cyan(string2color):
    return color_string(string2color, "1;36")
def magenta(string2color):
    return color_string(string2color, "1;35")
def bold(string2color):
    return color_string(string2color, "1;1")

SPACE = 4 * " "
def csprint(level, *args):
    print(level*SPACE + blue(":CS-"+str(level)+":"), " ".join(map(str,args)))
