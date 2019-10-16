# This class can be taken as an abstract template for derivation. It can be used
# for default passive behaviour.

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.SwimmingDEMApplication as SDEM
from . import recoverer

class ManufacturedSolutionRecoverer(recoverer.DerivativesRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.DerivativesRecoverer.__init__(self, project_parameters, model_part)

class CasasSolutionRecoverer(ManufacturedSolutionRecoverer):
    def __init__(self, project_parameters, model_part):
        super(CasasSolutionRecoverer, self).__init__(project_parameters, model_part)
        self.model_part = model_part
        self.settings = project_parameters["fluid_parameters"]["processes"]["boundary_conditions_process_list"][1]["Parameters"]["benchmark_parameters"]

    def RecoverFluidFractionGradient(self):
        from KratosMultiphysics.SwimmingDEMApplication.custom_body_force.casas_fluid_fraction_solution import CasasFluidFractionSolution
        self.fluid_fraction_field = CasasFluidFractionSolution(self.settings)
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        for node in self.model_part.Nodes:
            alpha1 = self.fluid_fraction_field.alpha1(time, node.X, node.Y, node.Z)
            alpha2 = self.fluid_fraction_field.alpha2(time, node.X, node.Y, node.Z)
            alpha3 = self.fluid_fraction_field.alpha3(time, node.X, node.Y, node.Z)
            fluid_fraction_gradient_field = [alpha1, alpha2, alpha3]
            node.SetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION_GRADIENT, fluid_fraction_gradient_field)

class StationarySolutionRecoverer(ManufacturedSolutionRecoverer):
    def __init__(self, project_parameters, model_part):
        super(StationarySolutionRecoverer, self).__init__(project_parameters, model_part)
        self.model_part = model_part
        self.settings = project_parameters["fluid_parameters"]["processes"]["boundary_conditions_process_list"][1]["Parameters"]["benchmark_parameters"]

    def RecoverFluidFractionGradient(self):
        from KratosMultiphysics.SwimmingDEMApplication.custom_body_force.stationary_fluid_fraction_solution import StationaryFluidFractionSolution
        self.fluid_fraction_field = StationaryFluidFractionSolution(self.settings)
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        for node in self.model_part.Nodes:
            alpha1 = self.fluid_fraction_field.alpha1(time, node.X, node.Y, node.Z)
            alpha2 = self.fluid_fraction_field.alpha2(time, node.X, node.Y, node.Z)
            alpha3 = self.fluid_fraction_field.alpha3(time, node.X, node.Y, node.Z)
            fluid_fraction_gradient_field = [alpha1, alpha2, alpha3]
            node.SetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION_GRADIENT, fluid_fraction_gradient_field)