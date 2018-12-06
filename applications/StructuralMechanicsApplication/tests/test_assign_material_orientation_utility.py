from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt, sin, cos, pi, exp, atan

class TestAssignMaterialOrientation(KratosUnittest.TestCase):
    # muting the output
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    def _apply_orthotropic_shell_material_properties(self,mp):
        #define properties
        # we specify only the properties we need (others are youngs modulus etc)
        num_plies = 3
        orthotropic_props = KratosMultiphysics.Matrix(num_plies,16)
        for row in range(num_plies):
            for col in range(16):
                orthotropic_props[row,col] = 0.0

        # Orthotropic mechanical moduli
        orthotropic_props[0,0] = 0.005 # lamina thickness
        orthotropic_props[0,2] = 2200  # density
        orthotropic_props[1,0] = 0.01  # lamina thickness
        orthotropic_props[1,2] = 1475  # density
        orthotropic_props[2,0] = 0.015 # lamina thickness
        orthotropic_props[2,2] = 520   # density

        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.SHELL_ORTHOTROPIC_LAYERS,orthotropic_props)

        cl = StructuralMechanicsApplication.LinearElasticOrthotropic2DLaw()

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)


    def _create_shell_nodes(self,mp):
       
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 0.3, 0.0, 0.0)
        mp.CreateNewNode(3, 0.3, 0.2, 0.0)
        mp.CreateNewNode(4, 0.0, 0.2, 0.0)
        mp.CreateNewNode(5, 0.8, 0.0, 0.0)
        mp.CreateNewNode(6, 0.8, 0.2, 0.0)
        mp.CreateNewNode(7, 0.8, 0.5, 0.0)
        mp.CreateNewNode(8, 0.3, 0.5, 0.0)
        mp.CreateNewNode(9, 0.0, 0.5, 0.0)

    def _create_shell_elements(self,mp,element_name = "ShellThinElementCorotational3D4N"):
    
        mp.CreateNewElement(element_name, 1, [1,2,3,4], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 2, [2,5,6,3], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 3, [3,6,7,8], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 4, [4,3,8,9], mp.GetProperties()[1])


    def test_orthotropic_shell_cog(self):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part_orthotropic_shells")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        mp.SetBufferSize(2)

        self._apply_orthotropic_shell_material_properties(mp)
        self._create_shell_nodes(mp)
        self._create_shell_elements(mp)

        composite_property_alignment_settings_11 = KratosMultiphysics.Parameters("""
        {
            "method": "simple",
            "echo_level"      : 1,
            "method_specific_settings" : {
                "global_fiber_direction" : [0.7071,0.7071,0]
            }
        }
        """)

        mat_orient_assign_util = StructuralMechanicsApplication.AssignMaterialOrientationUtility(mp)
        mat_orient_assign_util.Execute(composite_property_alignment_settings_11)
        for elem in mp.Elements:
            angle_theta = elem.GetValue(StructuralMechanicsApplication.MATERIAL_ORIENTATION_ANGLE)
            print(angle_theta)
            self.assertAlmostEqual(pi/4, angle_theta)
            

if __name__ == '__main__':
    KratosUnittest.main()


