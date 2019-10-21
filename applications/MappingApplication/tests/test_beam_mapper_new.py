from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM

import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.KratosUnittest as KratosUnittest

mdpa_file_name_beam    = "mdpa_files/beam_mesh"
mdpa_file_name_surface = "mdpa_files/surface_mesh"

def WriteGiDOutput(model_part):
    from gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(model_part,
        "gid_output_"+model_part.Name,
        KM.Parameters("""
            {
                "result_file_configuration" : {
                    "gidpost_flags": {
                        "GiDPostMode"           : "GiD_PostAscii",
                        "WriteDeformedMeshFlag" : "WriteUndeformed",
                        "WriteConditionsFlag"   : "WriteConditions",
                        "MultiFileFlag"         : "SingleFile"
                    },
                    "nodal_results"       : ["DISPLACEMENT"],
                    "gauss_point_results" : []
                }
            }
            """)
        )

    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

class TestBeamMapper(KratosUnittest.TestCase):
    def setUp(self):
        self.current_model = KM.Model()
        self.model_part_beam = self.current_model.CreateModelPart("beam")
        self.model_part_surface = self.current_model.CreateModelPart("surface")

        # list of variables involved in the Mapper-Tests
        self.model_part_beam.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.model_part_beam.AddNodalSolutionStepVariable(KM.ROTATION)

        self.model_part_surface.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

        KM.ModelPartIO(mdpa_file_name_beam).ReadModelPart(self.model_part_beam)
        KM.ModelPartIO(mdpa_file_name_surface).ReadModelPart(self.model_part_surface)
    
    def addDofs(self):
        for node in self.model_part_beam.Nodes:
            node.AddDof(KM.DISPLACEMENT_X)
            node.AddDof(KM.DISPLACEMENT_Y)
            node.AddDof(KM.DISPLACEMENT_Z)
            node.AddDof(KM.ROTATION_X)
            node.AddDof(KM.ROTATION_Y)
            node.AddDof(KM.ROTATION_Z)
        for node in self.model_part_surface.Nodes:
            node.AddDof(KM.DISPLACEMENT_X)
            node.AddDof(KM.DISPLACEMENT_Y)
            node.AddDof(KM.DISPLACEMENT_Z)

    def test_beam_mapper(self):
        mapper_settings = KM.Parameters("""{
            "mapper_type": "beam_mapper",
            "echo_level" : 3
        }""")

        self.mapper = KratosMapping.MapperFactory.CreateMapper(self.model_part_beam, self.model_part_surface, mapper_settings)
        #WriteGiDOutput(self.model_part_surface)
        for node in self.model_part_beam.Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT_Y, 10*node.X)
            node.SetSolutionStepValue(KM.ROTATION_X, 10*node.X)
            
        self.mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT)
        #self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)
        #self.mapper.InverseMap(KM.DISPLACEMENT, KM.TEMPERATURE)
        #self.mapper.InverseMap(KM.DISPLACEMENT_X, KM.TEMPERATURE)
        #self.mapper.InverseMap(KM.DISPLACEMENT_Y, KM.TEMPERATURE)
        #self.mapper.InverseMap(KM.DISPLACEMENT_Z, KM.TEMPERATURE)

        #WriteGiDOutput(self.model_part_beam)
        WriteGiDOutput(self.model_part_surface)


if __name__ == '__main__':
    KratosUnittest.main()