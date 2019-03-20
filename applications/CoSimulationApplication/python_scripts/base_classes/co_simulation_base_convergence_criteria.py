from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
data_structure = cs_tools.cs_data_structure
from decimal import Decimal
import numpy as np

def Create(settings, data_object):
    return CoSimulationConvergenceCriteria(settings, data_object)

class CoSimulationConvergenceCriteria(object):
    ## __init__ : Initializer
    #
    #  @param settings: setting of the convergence criteria
    #  @param data: data to which convergence criteria is to be attached
    def __init__(self, settings, data):
        default_settings = self._GetDefaultSettings()
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.data = data
        self.echo_level = 3
        self.data_name = self.data.name
        self.abs_tolerance = self.settings["abs_tolerance"].GetDouble()
        self.rel_tolerance = self.settings["rel_tolerance"].GetDouble()
        self.iteration = 0
        self.initial_residual_norm = 0.0

        self.data_prev_iter = []
        self.data_current_iter = []

    def _GetDefaultSettings(self):
        default_setting = data_structure.Parameters("""
            {
                "solver"        : "",
                "data_name"     : "",
                "abs_tolerance" : 1e-5,
                "rel_tolerance" : 1e-5,
                "echo_level"    : 3
            }
        """)
        return default_setting

    ## Initialize : Initialize function of the class
    #                   To be called once at the beginning of the simulation.
    #
    def Initialize(self):
        pass

    ## Finalize : Finalize function of the class
    #                   To be called once at the end of the simulation.
    #
    def Finalize(self):
        pass

    ## InitializeSolutionStep : InitializeSolutionStep function of the class.
    #                           To be called once at the beginning of the SolutionStep.
    #
    def InitializeSolutionStep(self):
        self.iteration = 0

    ## FinalizeSolutionStep : FinalizeSolutionStep function of the class.
    #                           To be called once at the end of the SolutionStep.
    #
    def FinalizeSolutionStep(self):
        pass

    ## InitializeCouplingIteration : InitializeCouplingIteration function of the class.
    #                           To be called once at the beginning of the non-linear coupling iteration.
    #
    def InitializeCouplingIteration(self):
        # storing the data for residual calculation
        self.data_prev_iter = self.data.GetNumpyArray()

    ## FinalizeCouplingIteration : FinalizeCouplingIteration function of the class.
    #                           To be called once at the end of the non-linear coupling iteration.
    #
    def FinalizeCouplingIteration(self):
        #self.data_current_iter = self.data.GetNumpyArray()
        self.iteration = self.iteration + 1


    ## IsConverged : IsConverged function of the class.
    #                  To be called called when the convergence is to be enquired for this criteria
    #
    def IsConverged(self):
        self.data_current_iter = self.data.GetNumpyArray()
        norm_current_data = np.linalg.norm(self.data_current_iter)/np.sqrt(self.data_current_iter.size)
        residual = self._CalculateResidual()
        abs_residual_norm = np.linalg.norm(residual) / np.sqrt(residual.size)
        if(abs_residual_norm == 0):
            abs_residual_norm = 1.0
        if(self.iteration == 1):
            self.initial_residual_norm = abs_residual_norm

        rel_residual_norm = abs_residual_norm / self.initial_residual_norm
        #rel_residual_norm = abs_residual_norm / norm_current_data

        is_converged = abs_residual_norm < self.abs_tolerance*10 or rel_residual_norm < self.rel_tolerance*10
        if self.echo_level > 1:
            if is_converged:
                info_msg = cs_tools.bcolors.GREEN+ "ACHIEVED"
            else:
                info_msg = cs_tools.bcolors.FAIL + "NOT ACHIEVED"
            info_msg += cs_tools.bcolors.ENDC
            cs_tools.PrintInfo(cs_tools.bcolors.HEADER + '\tConvergence for "' + cs_tools.bcolors.BOLD + self.data_name +'"', info_msg)
        if self.echo_level > 2:
            info_msg  = cs_tools.bcolors.HEADER + "\tabs_norm" + " = " + cs_tools.bcolors.BOLD + str("{:.3E}".format(Decimal(abs_residual_norm))) + " | "
            info_msg += "abs_tol" + " = " + cs_tools.bcolors.BOLD + str(self.abs_tolerance)
            info_msg += " || " + "rel_norm" + " = " + cs_tools.bcolors.BOLD + str("{:.3E}".format(Decimal(rel_residual_norm))) + " | "
            info_msg += "rel_tol = " + cs_tools.bcolors.BOLD + str(self.rel_tolerance) + cs_tools.bcolors.ENDC
            cs_tools.PrintInfo(info_msg)

        return is_converged

    ## PrintInfo : PrintInfo function of the class.
    #                  Prints the information about the object
    #
    def PrintInfo(self):
        cs_tools.PrintInfo("Convergence criteria for data : ", self.data_name)
        cs_tools.PrintInfo("Absolute Tolerace  : ", self.abs_tolerance)
        cs_tools.PrintInfo("Relative Tolerace  : ", self.rel_tolerance)
        print()

    ## Check : Check function of the class.
    #
    def Check(self):
        cs_tools.PrintInfo("Check from Base Convergence Criteria")


    ## SetEchoLevel : Function to set the echo level of this class.
    #                   Default is 3
    #
    def SetEchoLevel(self, level):
        self.echo_level = level


    ## _Name : Gets the name of the class
    #
    def _Name(self):
        return self.__class__.__name__

    ## _CalculateResidual : Calculates residual of the data specified in the settings
    #                       Numpy can be used in the variants of this class.
    #                       residual = data_in_current_iter - data_in_previous_iter
    #
    #  @param self            The object pointer
    def _CalculateResidual(self):
        return self.data_current_iter - self.data_prev_iter
