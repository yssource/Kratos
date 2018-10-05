from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_solvers.co_simulation_base_solver import CoSimulationBaseSolver
from co_simulation_tools import RecursivelyValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

available_mdof_models = ["sdof",
                         "cantilever_shear_2d",
                         "cantilever_eb_beam_2d",
                         "bridge_2dof"]

available_schemes = ["forward_euler1",
                     "euler12",
                     "generalized_alpha",
                     "bdf2"]

def CreateSolver(cosim_solver_settings, level):
    return MDoFSolver(cosim_solver_settings, level)

class MDoFSolver(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings, level):
        super(MDoFSolver, self).__init__(cosim_solver_settings, level)

        input_file_name = self.cosim_solver_settings["input_file"]
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as ProjectParameters:
            parameters = json.load(ProjectParameters)

        # creating model using a certain module
        model_type = parameters["model_parameters"]["type"]
        if not(model_type in available_mdof_models):
            err_msg  = 'The requested mdof model "' + model_type + '" is not available!\n'
            err_msg += 'The following mdof models are available:\n'
            for avail_model in available_mdof_models:
                err_msg += "\t" + avail_model + "\n"
            raise NameError(err_msg)
        model_module = __import__("mdof_" + model_type + "_model")
        self.model = model_module.CreateModel(parameters["model_parameters"])

        # creating model using a certain module
        scheme_type = parameters["time_integration_scheme_parameters"]["type"]
        if not(scheme_type in available_schemes):
            err_msg  = 'The requested time integration scheme "' + model_type + '" is not available!\n'
            err_msg += 'The following time integration schemes are available:\n'
            for avail_scheme in available_schemes:
                err_msg += "\t" + avail_scheme + "\n"
            raise NameError(err_msg)
        scheme_module = __import__("time_integration_" + scheme_type + "_scheme")
        self.scheme = scheme_module.CreateScheme(parameters["time_integration_scheme_parameters"])

        default_solver_settings = {
                    "buffer_size": 3
                }
        default_output_settings = {
                    "file_name": "results_mdof.dat"
                }

        # solver settings
        RecursivelyValidateAndAssignDefaults(default_solver_settings, parameters["solver_parameters"])
        # output settings
        RecursivelyValidateAndAssignDefaults(default_output_settings, parameters["output_parameters"])

        self.buffer_size = parameters["solver_parameters"]["buffer_size"]
        self.output_file_name = parameters["output_parameters"]["file_name"]

        # PMT this might not be needed
        self.load_vector = None

    def Initialize(self):
        if os.path.isfile(self.output_file_name):
            os.remove(self.output_file_name)

        # 1st dimension: variables: disp, acc, vel
        # 2nd dimension: buffer size
        # 3rd dimension: number of dofs
        self.x = np.zeros((3, self.buffer_size, len(self.model.u0)))

        # initialize scheme parameters
        self.scheme.Initialize(self.model)

        self.dx = np.array([self.scheme.GetDisplacement(),
                            self.scheme.GetVelocity(),
                            self.scheme.GetAcceleration()])

        # PMT is this needed?
        self.load_vector = self.scheme.GetLoad()

    def Predict(self):
        return self.scheme.Predict()

    def OutputSolutionStep(self):
        # PMT: check if this syntax is still ok for MDoF
        with open(self.output_file_name, "a") as results_sdof:
            #outputs displacements
            results_sdof.write(str(self.time) + "\t" + " ".join(str(value) for value in self.dx[0]) + "\n")

    def AdvanceInTime(self, current_time):
        # PMT: check if this syntax is still ok for MDoF

        # similar to the Kratos CloneTimeStep function
        # advances values along the buffer axis (so rolling columns) using numpy's roll

        self.x = np.roll(self.x,1,axis=1)
        # overwriting at the buffer_idx=0 the newest values
        #buffer_idx = 0
        #self.x[:,buffer_idx] = self.dx

        self.dx[0] = self.scheme.GetPreviousDisplacement()
        self.dx[1] = self.scheme.GetPreviousVelocity()
        self.dx[2] = self.scheme.GetPreviousAcceleration()

        self.x[0,0,:] = self.dx[0]
        self.x[1,0,:] = self.dx[1]
        self.x[2,0,:] = self.dx[2]

        # update displacement, velocity and acceleration
        self.scheme.AdvanceScheme()
        self.time = current_time + self.scheme.dt

        return self.time

    def SolveSolutionStep(self):
        # sys of eq reads: LHS * u1 = RHS
        self.scheme.Solve(self.model)
        # from u1 determine v1, a1
        self.scheme.UpdateDerivedValues()

    def GetBufferSize(self):
        return self.buffer_size

    def GetDeltaTime(self):
        return self.scheme.dt

    def GetSolutionStepValue(self, identifier, buffer_idx=0):
        #
        # PMT: check if this syntax is still ok for MDoF
        #

        if identifier == "DISPLACEMENT":
            return self.x[:,buffer_idx][0]
        elif identifier == "VELOCITY":
            return self.x[:,buffer_idx][1]
        elif identifier == "ACCELERATION":
            return self.x[:,buffer_idx][2]
        else:
            raise Exception("Identifier is unknown!")

    def SetSolutionStepValue(self, identifier, value, buffer_idx=0):
        #
        # PMT: check if this syntax is still ok for MDoF
        #

        if identifier == "DISPLACEMENT":
            self.x[:,buffer_idx][0] = value
        elif identifier == "VELOCITY":
            self.x[:,buffer_idx][1] = value
        elif identifier == "ACCELERATION":
            self.x[:,buffer_idx][2] = value
        else:
            raise Exception("Identifier is unknown!")

    def SetData(self, identifier, data):
        #
        # PMT: check if this syntax is still ok for MDoF
        #

        if identifier == "LOAD":
            # last index is the external force
            self.load_vector[-1] = data
        elif identifier == "DISPLACEMENT":
            # first index is displacement
            # maybe use buffer index
            self.SetSolutionStepValue("DISPLACEMENT", data,0)
        else:
            raise Exception("Identifier is unknown!")

    def GetData(self, identifier):
        #
        # PMT: check if this syntax is still ok for MDoF
        #

        if identifier == "LOAD":
            # last index is the external force
            return self.load_vector[-1]
        elif identifier == "DISPLACEMENT":
            # first index is displacement
            return self.GetSolutionStepValue("DISPLACEMENT",0)
        else:
            raise Exception("Identifier is unknown!")

    def _GetIOName(self):
        return "mdof"

    def _Name(self):
        return self.__class__.__name__