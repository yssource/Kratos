# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.SwimmingDEMApplication
from importlib import import_module
import numpy as np
from KratosMultiphysics import Vector

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyCustomBodyForceProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyCustomBodyForceProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"          : "please_specify_model_part_name",
                "variable_name"            : "BODY_FORCE",
                "benchmark_name"           : "custom_body_force.vortex",
                "benchmark_parameters"     : {},
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
            }
            """
            )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[self.settings["model_part_name"].GetString()]
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())

        benchmark_module = import_module(self.settings["benchmark_name"].GetString())
        self.benchmark = benchmark_module.CreateManufacturedSolution(self.settings["benchmark_parameters"])

        self.compute_error = self.settings["compute_nodal_error"].GetBool()

        self.print_output = self.settings["print_convergence_output"].GetBool()
        if self.print_output:
            from custom_body_force.hdf5_output_tool import Hdf5OutputTool
            self.output_process = Hdf5OutputTool(model, self.settings["output_parameters"])

    def ExecuteBeforeSolutionLoop(self):
        current_time = 0.0
        x = np.array([node.X for node in self.model_part.Nodes])
        y = np.array([node.Y for node in self.model_part.Nodes])
        z = np.array([node.Z for node in self.model_part.Nodes])

        value = np.array([self.benchmark.Velocity(current_time, x, y, z) for x, y, z in zip(x, y, z)])

        iterator = 0
        for node in self.model_part.Nodes:
            vel_value = Vector(list(value[iterator]))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, vel_value)
            node.SetValue(KratosMultiphysics.Y, vel_value)
            iterator += 1

    def ExecuteInitializeSolutionStep(self):
        self._SetBodyForce()

    def ExecuteFinalizeSolutionStep(self):
        if self.compute_error:
            self._ComputeVelocityError()

    def ExecuteBeforeOutputStep(self):
        self._ComputeVelocityBenchmark()

    def ExecuteFinalize(self):
        if self.compute_error:
            rel_err = self._SumNodalError()
            if self.print_output:
                self.output_process.WriteBodyForceAttributes(self.settings["benchmark_parameters"])
                self.output_process.WriteAverageRelativeError(rel_err)

    def _SetBodyForce(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        x = np.array([node.X for node in self.model_part.Nodes])
        y = np.array([node.Y for node in self.model_part.Nodes])
        z = np.array([node.Z for node in self.model_part.Nodes])

        value = np.array([self.benchmark.BodyForce(current_time, x, y, z) for x, y, z in zip(x, y, z)])
        iterator = 0
        for node in self.model_part.Nodes:
            var_value = Vector(list(value[iterator]))
            node.SetSolutionStepValue(self.variable, var_value)
            iterator += 1

    def _ComputeVelocityError(self):
        epsilon = 1e-16
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        x = np.array([node.X for node in self.model_part.Nodes])
        y = np.array([node.Y for node in self.model_part.Nodes])
        z = np.array([node.Z for node in self.model_part.Nodes])

        value = np.array([self.benchmark.Velocity(current_time, x, y, z) for x, y, z in zip(x, y, z)])
        iterator = 0

        for node in self.model_part.Nodes:
            fem_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            exact_vel = Vector(list(value[iterator]))
            fem_vel_modulus = (fem_vel[0]**2 + fem_vel[1]**2 + fem_vel[2]**2)**0.5
            exact_vel_modulus = (exact_vel[0]**2 + exact_vel[1]**2 + exact_vel[2]**2)**0.5
            error = abs(fem_vel_modulus - exact_vel_modulus)
            error = error / abs(exact_vel_modulus + epsilon)
            node.SetValue(KratosMultiphysics.NODAL_ERROR, error)
            node.SetSolutionStepValue(KratosMultiphysics.SwimmingDEMApplication.EXACT_VELOCITY, exact_vel)
            if node.X == min(x) or node.X == max(x)  or node.Y == min(y) or node.Y == max(y):
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, exact_vel)
            iterator += 1

    def _CopyVelocityAsNonHistorical(self):
        for node in self.model_part.Nodes:
            vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            node.SetValue(KratosMultiphysics.VELOCITY, vel)

    def _ComputeVelocityBenchmark(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        x = np.array([node.X for node in self.model_part.Nodes])
        y = np.array([node.Y for node in self.model_part.Nodes])
        z = np.array([node.Z for node in self.model_part.Nodes])

        value = np.array([self.benchmark.Velocity(current_time, x, y, z) for x, y, z in zip(x, y, z)])
        iterator = 0
        for node in self.model_part.Nodes:
            vel_value = Vector(list(value[iterator]))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, vel_value)
            node.SetValue(KratosMultiphysics.Y, vel_value)
            iterator += 1

    def _SumNodalError(self):
        err_sum = 0.0
        for node in self.model_part.Nodes:
            err_sum = err_sum + node.GetValue(KratosMultiphysics.NODAL_ERROR)
        rel_err = err_sum / self.model_part.Nodes.__len__()
        KratosMultiphysics.Logger.PrintInfo("Benchmark", "The nodal error average is : ", rel_err)
        return rel_err
