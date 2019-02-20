from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

def CreateMeshMotionModelPart(main_model_part, element_name):
    mesh_motion_elements_model_part_name = main_model_part.Name + "_" + element_name
    mesh_motion_elements_model_part = main_model_part.CreateSubModelpart(mesh_motion_elements_model_part_name)

    connectivity_preserve_modeler_settings = KM.Parameters("""{
        "element_name"              : \"""" + element_name"""\",
        "duplicate_sub_model_parts" : false
    }""")

    modeler = KM.ConnectivityPreserveModeler()
    modeler.GenerateModelPart(main_model_part,
                              mesh_motion_elements_model_part,
                              connectivity_preserve_modeler_settings)

    return mesh_motion_elements_model_part