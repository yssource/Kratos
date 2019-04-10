from __future__ import print_function, absolute_import, division

import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication
from python_solver import PythonSolver

from KM.kratos_utilities import IsApplicationAvailable
if IsApplicationAvailable("RANSConstitutiveLawsApplication"):
    import KM.RANSConstitutiveLawsApplication as RANS


def CreateTurbulenceModel(model, settings):
if not IsApplicationAvailable("RANSConstitutiveLawsApplication"):
        raise Exception(
            "Using a turbulence model requires the RANSConstitutiveLawsApplication, Please re-install/re-compile with RANSConstitutiveLawsApplication."
        )

    model_type = settings["model_type"].GetString()
    if model_type == "k_epsilon":
        from RANS.k_epsilon_configuration import KEpsilonConfiguration
        return KEpsilonConfiguration(model, settings)
    else:
        raise Exception("Unknown turbulence model_type: " + model_type)


class TurbulenceModelConfiguration(PythonSolver):
    def __init__(self, model, settings):
        super(TurbulenceModelConfiguration, self).__init__(model, settings)

        self.model_elements = []
        self.model_conditions = []

    def CreateTurbulenceModelParts(self, original_model_part):
        self.domain_size = original_model_part.ProcessInfo[KM.DOMAIN_SIZE]

        tubulence_model_parts_list = []
        connectivity_preserve_modeler = KM.ConnectivityPreserveModeler()
        for elem, cond in zip(self.model_elements, self.model_conditions):
            model_part_name = "TurbulenceModelPart_" + elem
            if not self.model.HasSubModelPart(model_part_name):
                model_part = self.model.CreateSubModelPart(model_part_name)
                element_name = "{0}{1}D{2}N".format(elem, self.domain_size,
                                                    self.domain_size + 1)
                condition_name = "{0}{1}D{2}N".format(cond, self.domain_size,
                                                      self.domain_size)
                connectivity_preserve_modeler.GenerateModelPart(
                    original_model_part, model_part, element_name,
                    condition_name)
                KM.Logger.PrintInfo(self.__name__, model_part_name + " created successfully.")
            else:
                KM.Logger.PrintInfo(self.__name__, model_part_name + " already exists. Using the existing model part.")
            
            tubulence_model_parts_list.append(self.model.GetSubModelPart(model_part_name))
        
        return tubulence_model_parts_list

    def AddVariables(self):
        pass

    def AddDofs(self):
        pass

    def GetTurbulenceModelProcess(self):
        raise Exception(
            "Calling base class TurbulenceModelConfiguration::GetTurbulenceModelProcess"
        )
