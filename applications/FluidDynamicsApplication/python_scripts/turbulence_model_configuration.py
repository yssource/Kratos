from __future__ import print_function, absolute_import, division

import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication
from python_solver import PythonSolver

from kratos_utilities import CheckIfApplicationsAvailable
if CheckIfApplicationsAvailable("RANSConstitutiveLawsApplication"):
    import KratosMultiphysics.RANSConstitutiveLawsApplication as RANS

def CreateTurbulenceModel(model, settings):
    if not CheckIfApplicationsAvailable("RANSConstitutiveLawsApplication"):
        raise Exception(
            "Using a turbulence model requires the RANSConstitutiveLawsApplication, Please re-install/re-compile with RANSConstitutiveLawsApplication."
        )

    model_type = settings["model_type"].GetString()
    if model_type == "k_epsilon":
        from evm_k_epsilon_configuration import TurbulenceKEpsilonConfiguration
        return TurbulenceKEpsilonConfiguration(model, settings)
    else:
        raise Exception("Unknown turbulence model_type: " + model_type)


class TurbulenceModelConfiguration(Kratos.Process):
    def __init__(self, model, settings):
        super(TurbulenceModelConfiguration, self).__init__()

        if (type(model) != KratosMultiphysics.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if (type(settings) != KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.model = model
        self.settings = settings

        self.echo_level = self.settings["echo_level"].GetInt()

        self.model_elements = []
        self.model_conditions = []

    def CreateTurbulenceModelParts(self):
        self.domain_size = self.fluid_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]

        tubulence_model_parts_list = []
        connectivity_preserve_modeler = Kratos.ConnectivityPreserveModeler()
        for elem, cond in zip(self.model_elements, self.model_conditions):
            model_part_name = "TurbulenceModelPart_" + elem
            if not self.model.HasModelPart(model_part_name):
                model_part = self.model.CreateModelPart(model_part_name)

                self.AddVariables(model_part)

                element_name = "{0}{1}D{2}N".format(elem, self.domain_size,
                                                    self.domain_size + 1)
                condition_name = "{0}{1}D{2}N".format(cond, self.domain_size,
                                                      self.domain_size)
                connectivity_preserve_modeler.GenerateModelPart(
                    self.fluid_model_part, model_part, element_name,
                    condition_name)
                Kratos.Logger.PrintInfo(self.__class__.__name__, model_part_name + " created successfully.")
            else:
                Kratos.Logger.PrintInfo(self.__class__.__name__, model_part_name + " already exists. Using the existing model part.")

            tubulence_model_parts_list.append(self.model.GetModelPart(model_part_name))

        return tubulence_model_parts_list

    def AddVariables(self, model_part):
        pass

    def AddDofs(self, model_part):
        pass

    def Execute(self):
        print(" Test Trampoline")

