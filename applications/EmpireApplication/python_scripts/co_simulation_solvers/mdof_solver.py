from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_solvers.co_simulation_base_solver import CoSimulationBaseSolver
from co_simulation_tools import RecursivelyValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

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
        model_module = __import__("mdof_" + parameters["model_parameters"]["type"] + "_model")
        model = model_module.CreateModel(parameters["model_parameters"])

        default_solver_settings = json.loads("""{
                    "buffer_size": 3
                }""")
        default_output_settings = json.loads("""{
                    "file_name": "results_mdof.dat"
                }""")

        RecursivelyValidateAndAssignDefaults(default_solver_settings, parameters["solver_parameters"])
        RecursivelyValidateAndAssignDefaults(default_output_settings, parameters["output_parameters"])

        ##
        #PMT: paramaters received from parameters JSON string ->

        self.buffer_size = parameters["solver_parameters"]["buffer_size"]
        self.output_file_name = parameters["output_parameters"]["file_name"]

        ##
        #PMT: to do for (dt, mM, mB, mK, pInf, vu0, vv0, va0)
        # time step
        self.dt = parameters["time_integration_scheme_parameters"]["time_step"]

        # mass, damping and spring stiffness
        self.M = model.m
        self.B = model.b
        self.K = model.k

        # what to do with model["nodal_coordinates?"]

        # generalized alpha parameters (to ensure unconditional stability, 2nd order accuracy)
        pInf = parameters["time_integration_scheme_parameters"]["settings"]["p_inf"]
        self.alphaM = (2.0 * pInf - 1.0) / (pInf + 1.0)
        self.alphaF = pInf / (pInf + 1.0)
        self.beta = 0.25 * (1 - self.alphaM + self.alphaF)**2
        self.gamma = 0.5 - self.alphaM + self.alphaF

        # coefficients for LHS
        self.a1h = (1.0 - self.alphaM) / (self.beta * self.dt**2)
        self.a2h = (1.0 - self.alphaF) * self.gamma / (self.beta * self.dt)
        self.a3h = 1.0 - self.alphaF

        # coefficients for mass
        self.a1m = self.a1h
        self.a2m = self.a1h * self.dt
        self.a3m = (1.0 - self.alphaM - 2.0 * self.beta) / (2.0 * self.beta)

        #coefficients for damping
        self.a1b = (1.0 - self.alphaF) * self.gamma / (self.beta * self.dt)
        self.a2b = (1.0 - self.alphaF) * self.gamma / self.beta - 1.0
        self.a3b = (1.0 - self.alphaF) * (0.5 * self.gamma / self.beta - 1.0) * self.dt

        # coefficient for stiffness
        self.a1k = -1.0 * self.alphaF

        # coefficients for velocity update
        self.a1v = self.gamma / (self.beta * self.dt)
        self.a2v = 1.0 - self.gamma / self.beta
        self.a3v = (1.0 - self.gamma / (2 * self.beta)) * self.dt

        # coefficients for acceleration update
        self.a1a = self.a1v / (self.dt * self.gamma)
        self.a2a = -1.0 / (self.beta * self.dt)
        self.a3a = 1.0 - 1.0 / (2.0 * self.beta)

		#structure
        # initial displacement, velocity and acceleration
        self.u0 = model.u0
        self.v0 = model.v0
        self.a0 = model.a0
        # initial force
        self.force = model.f0

        # initial displacement, velocity and acceleration
        self.u1 = self.u0
        self.v1 = self.v0
        self.a1 = self.a0

		# force from a previous time step (initial force)
        self.f0 = np.dot(self.M,self.a0) + np.dot(self.B,self.v0) + np.dot(self.K,self.u0)
        self.f1 = np.dot(self.M,self.a1) + np.dot(self.B,self.v1) + np.dot(self.K,self.u1)

    def Initialize(self):
        # 1st dimension: variables: disp, acc, vel
        # 2nd dimension: buffer size
        # 3rd dimension: number of dofs
        self.x = np.zeros((3, self.buffer_size, len(self.u0)))

        # PMT: check if this syntax is still ok for MDoF
        self.dx = np.array([self.u0,
                            self.v0,
                            self.a0])

        #x and dx contain: [displacement, velocity, acceleration]

        if os.path.isfile(self.output_file_name):
            os.remove(self.output_file_name)

        #apply external load as an initial impulse
        self.load_vector = self.force

    ## PMT this is new for MDoF, maybe do it for SDoF as well
    def Predict(self):
        return 2.0 * self.u1 - self.u0

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

        self.dx[0] = self.u1
        self.dx[1] = self.v1
        self.dx[2] = self.a1

        self.x[0,0,:] = self.dx[0]
        self.x[1,0,:] = self.dx[1]
        self.x[2,0,:] = self.dx[2]
        ##
        # PMT: this might be needed
        # update displacement, velocity and acceleration
        self.u0 = self.u1
        self.v0 = self.v1
        self.a0 = self.a1

        # update the force
        self.f0 = self.f1

        self.time = current_time + self.dt
        return self.time

    def SolveSolutionStep(self):
        # sys of eq reads: LHS * u1 = RHS

        #F = (1.0 - self.alphaF) * f1 + self.alphaF * self.f0
        F = (1.0 - self.alphaF) * self.load_vector + self.alphaF * self.f0

        LHS = self.a1h * self.M + self.a2h * self.B + self.a3h * self.K
        RHS = np.dot(self.M,(self.a1m * self.u0 + self.a2m * self.v0 + self.a3m * self.a0))
        RHS += np.dot(self.B,(self.a1b * self.u0 + self.a2b * self.v0 + self.a3b * self.a0))
        RHS += np.dot(self.a1k * self.K, self.u0) + F

        # update self.f1
        #self.f1 = f1
        self.f1 = self.load_vector

        # updates self.u1,v1,a1
        self.u1 = np.linalg.solve(LHS, RHS)
        self.v1 = self.a1v * (self.u1 - self.u0) + self.a2v * self.v0 + self.a3v * self.a0
        self.a1 = self.a1a * (self.u1 - self.u0) + self.a2a * self.v0 + self.a3a * self.a0

    def GetBufferSize(self):
        return self.buffer_size

    def GetDeltaTime(self):
        return self.dt

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