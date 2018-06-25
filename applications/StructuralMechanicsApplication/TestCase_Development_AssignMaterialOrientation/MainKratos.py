import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

from structural_mechanics_analysis import StructuralMechanicsAnalysis

class AssignMaterialDirectionStage(StructuralMechanicsAnalysis):

    def RunSolutionLoop(self):
        while self.time < self.end_time:
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            # self._GetSolver().SolveSolutionStep() # Do nothing bcs I am only interested in assigning the material direction
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def Initialize(self):
        super(AssignMaterialDirectionStage, self).Initialize()

        composite_property_alignment_settings_11 = KratosMultiphysics.Parameters("""
        {
            "method": "simple",
            "echo_level"      : 1,
            "method_specific_settings" : {
                "global_fiber_direction" : [0,0,1]
            }
        }
        """)
        composite_property_alignment_settings_12 = KratosMultiphysics.Parameters("""
        {
            "method": "simple",
            "echo_level"      : 1,
            "method_specific_settings" : {
                "cs_axis_1" : [0,0,1],
                "cs_axis_2" : [1,0,0],
                "cs_normal_axis" : 3
            }
        }
        """)
        composite_property_alignment_settings_13 = KratosMultiphysics.Parameters("""
        {
            "method": "simple",
            "echo_level"      : 0,
            "method_specific_settings" : {
                "cs_axis_1" : [0,0,1],
                "cs_axis_2" : [1,0,0],
                "cs_rotation_angle" : 15,
                "cs_normal_axis" : 3
            }
        }
        """)

        composite_property_alignment_settings_21 = KratosMultiphysics.Parameters("""
        {
            "method": "advanced",
            "echo_level"      : 2,
            "method_specific_settings" : {
                "global_fiber_direction" : [0,0,1],
                "normal_vector"   : [1,0,0],
                "smoothness_level" : 1
            }
        }
        """)

        composite_property_alignment_settings_22 = KratosMultiphysics.Parameters("""
        {
            "method": "advanced",
            "echo_level"      : 0,
            "method_specific_settings" : {
                "global_fiber_direction" : [0,0,1],
                "normal_vector"   : [1,0,0],
                "smoothness_level" : 2
            }
        }
        """)

        composite_property_alignment_settings_23 = KratosMultiphysics.Parameters("""
        {
            "method": "advanced",
            "echo_level"      : 0,
            "method_specific_settings" : {
                "global_fiber_direction" : [0,0,1],
                "normal_vector"   : [1,0,0],
                "smoothness_level" : 3
            }
        }
        """)
        main_model_part = self.model["Structure"]
        mat_orient_assign_util = StructuralMechanicsApplication.AssignMaterialOrientationUtility(main_model_part)
        mat_orient_assign_util.Execute(composite_property_alignment_settings_11)
        mat_orient_assign_util.WriteFiberAngles("fiber_angles.dat")


# Printing with color: http://ozzmaker.com/add-colour-to-text-in-python/
# https://stackoverflow.com/questions/287871/print-in-terminal-with-colors
def print_with_color(color_code, *args):
    print("\033["+color_code+"m", " ".join(map(str,args)), "\033[0m")

def print_yellow(*args):
    print_with_color("1;33", *args)

def print_red(*args):
    print_with_color("1;31", *args)

def print_bold(*args):
    print_with_color("1;37", *args)


if __name__ == "__main__":
    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = AssignMaterialDirectionStage(model, parameters)
    simulation.Run()
