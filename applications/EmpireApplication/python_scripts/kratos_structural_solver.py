from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Importing the base class
from co_simulation_base_solver import CoSimulationBaseSolver

# Other imports
from structural_mechanics_analysis import StructuralMechanicsAnalysis

def CreateSolver(model, custom_settings):
    return KratosStructuralSolver(model, custom_settings)

class KratosStructuralSolver(CoSimulationBaseSolver):
    """The base class for the Python Solvers in the applications
    Changes to this BaseClass have to be discussed first!
    """
    def __init__(self, model, settings):
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

