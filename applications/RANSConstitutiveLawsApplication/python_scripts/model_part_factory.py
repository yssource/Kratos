import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSConstitutiveLawsApplication as KratosRANS


def CreateDuplicateModelPart(origin_modelpart,
                             destination_modelpart_name,
                             element_name,
                             condition_name):
    domain_size = origin_modelpart.ProcessInfo[Kratos.DOMAIN_SIZE]
    model = origin_modelpart.GetModel()
    connectivity_preserve_modeler = Kratos.ConnectivityPreserveModeler()

    if not model.HasModelPart(destination_modelpart_name):
        model_part = model.CreateModelPart(destination_modelpart_name)
        KratosRANS.RansVariableUtils().CopyNodalSolutionStepVariablesList(
            origin_modelpart, model_part)

        connectivity_preserve_modeler.GenerateModelPart(
            origin_modelpart, model_part, element_name, condition_name)

    return model.GetModelPart(destination_modelpart_name)