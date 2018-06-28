from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.EmpireApplication as KratosEmpire

# Importing the base class
from co_simulation_base_solver import CoSimulationBaseSolver

# Other imports
from io_factory import IOFactory


def CreateSolver(cosim_solver_settings):
    return KratosEmpireSolver(cosim_solver_settings)

class KratosEmpireSolver(CoSimulationBaseSolver):
    """The base class for the Python Solvers in the applications
    Changes to this BaseClass have to be discussed first!
    """
    def __init__(cosim_solver_settings):
        """The constructor of the PythonSolver-Object.

        It is intended to be called from the constructor
        of deriving classes:
        super(DerivedSolver, self).__init__(settings)

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        settings -- The solver settings used
        """
        pass

