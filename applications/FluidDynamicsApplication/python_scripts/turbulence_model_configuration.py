from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

from KratosMultiphysics.kratos_utilities import IsApplicationAvailable
if IsApplicationAvailable("RANSConstitutiveLawsApplication"):
    import KratosMultiphysics.RANSConstitutiveLawsApplication as RANS

def CreateTurbulenceModel(model, settings):

    if not IsApplicationAvailable("RANSConstitutiveLawsApplication"):
        raise Exception("Using a turbulence model requires the RANSConstitutiveLawsApplication, which is unavailable.")

    model_type = settings["model_type"].GetString()
    if model_type == "k_epsilon":
        from RANS.k_epsilon_configuration import KEpsilonConfiguration
        return KEpsilonConfiguration(model,settings)
    else:
        raise Exception("Unknown turbulence model_type: "+model_type)

class TurbulenceModelConfiguration(object):

    def __init__(self, model, parameters):
        self.model = model
        self.default_settings = KratosMultiphysics.Parameters(r'''{
            "model_type" : "",
            "fluid_model_part" : "PLEASE_SPECIFY_FLUID_MODEL_PART",
            "model_settings" : {}
        }''')

        self.fluid_model_part = self.model.GetSubModelPart(self.settings["fluid_model_part"].GetString())

        self.model_elements = []
        self.model_conditions = []

    def ValidateInputSettings(self, parameters, defaults):
        return parameters.ValidateAndAssignDefaults(self.default_settings)

    def CreateTurbulenceModelParts(self):
        domain_size = self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        connectivity_preserve_modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        for elem,cond in zip(self.model_elements, self.model_conditions):
            if not self.model.HasSubModelPart("TurbulenceModelPart_"+elem):
                model_part = self.model.CreateSubModelPart("TurbulenceModelPart_"+elem)
                element_name = "{0}{1}D{2}N".format(elem,domain_size,domain_size+1)
                condition_name = "{0}{1}D{2}N".format(cond,domain_size,domain_size)
                connectivity_preserve_modeler.GenerateModelPart(self.fluid_model_part, model_part, element_name, condition_name)



    def AddVariables(self):
        self.fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)
        self.fluid_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
        self.fluid_model_part.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY)

    def AddDofs(self):
        pass

    def GetTurbulenceModelProcess(self):
        raise Exception("Calling base class TurbulenceModelConfiguration::GetTurbulenceModelProcess")