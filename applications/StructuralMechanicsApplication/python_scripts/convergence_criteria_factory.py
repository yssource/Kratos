from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication


# Convergence criteria class
class convergence_criterion:
    def __init__(self, convergence_criterion_parameters):
        # Note that all the convergence settings are introduced via a Kratos parameters object.

        D_RT = convergence_criterion_parameters["displacement_relative_tolerance"].GetDouble()
        D_AT = convergence_criterion_parameters["displacement_absolute_tolerance"].GetDouble()
        R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
        R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()

        echo_level = convergence_criterion_parameters["echo_level"].GetInt()
        convergence_crit = convergence_criterion_parameters["convergence_criterion"].GetString()

        if(echo_level >= 1):
            KratosMultiphysics.Logger.PrintInfo("::[Mechanical Solver]:: ", "CONVERGENCE CRITERION : " +
                  convergence_criterion_parameters["convergence_criterion"].GetString())

        rotation_dofs = False
        if(convergence_criterion_parameters.Has("rotation_dofs")):
            if(convergence_criterion_parameters["rotation_dofs"].GetBool()):
                rotation_dofs = True


        # Convergence criteria if there are rotation DOFs in the problem
        if(rotation_dofs == True):
            conv_crit_settings_res = KratosMultiphysics.Parameters("""
            {
                "basis_vector_type" : "residual",
                "variables_to_separate" : ["DISPLACEMENT"],
                "relative_convergence_tolerances" : [1e-6],
                "absolut_convergence_tolerances" : [1e-11]
            }
            """)
            conv_crit_settings_res["relative_convergence_tolerances"][0].SetDouble(R_RT)
            conv_crit_settings_res["absolut_convergence_tolerances"][0].SetDouble(R_AT)

            conv_crit_settings_disp = KratosMultiphysics.Parameters("""
            {
                "basis_vector_type" : "solution_update",
                "variables_to_separate" : ["DISPLACEMENT"],
                "relative_convergence_tolerances" : [1e-6],
                "absolut_convergence_tolerances" : [1e-11]
            }
            """)
            conv_crit_settings_disp["relative_convergence_tolerances"][0].SetDouble(D_RT)
            conv_crit_settings_disp["absolut_convergence_tolerances"][0].SetDouble(D_AT)

            other_dofs_name = "ROTATION"

            if(convergence_crit == "displacement_criterion"):
                self.mechanical_convergence_criterion = StructuralMechanicsApplication.GeneralConvergenceCriteria(D_RT, D_AT, conv_crit_settings_res)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)

            elif(convergence_crit == "residual_criterion"):
                self.mechanical_convergence_criterion = StructuralMechanicsApplication.GeneralConvergenceCriteria(D_RT, D_AT, conv_crit_settings_disp)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)

            elif(convergence_crit == "and_criterion"):
                Displacement = StructuralMechanicsApplication.GeneralConvergenceCriteria(D_RT, D_AT, conv_crit_settings_disp)
                Displacement.SetEchoLevel(echo_level)
                Residual = StructuralMechanicsApplication.GeneralConvergenceCriteria(D_RT, D_AT, conv_crit_settings_res)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)

            elif(convergence_crit == "or_criterion"):
                Displacement = StructuralMechanicsApplication.GeneralConvergenceCriteria(D_RT, D_AT, conv_crit_settings_disp)
                Displacement.SetEchoLevel(echo_level)
                Residual = StructuralMechanicsApplication.GeneralConvergenceCriteria(D_RT, D_AT, conv_crit_settings_res)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
            else:
                err_msg =  "The requested convergence criterion \"" + convergence_crit + "\" is not available!\n"
                err_msg += "Available options are: \"displacement_criterion\", \"residual_criterion\", \"and_criterion\", \"or_criterion\""
                raise Exception(err_msg)

        # Convergence criteria without rotation DOFs
        else:
            conv_crit_settings_res = KratosMultiphysics.Parameters("""
            {
                "basis_vector_type" : "residual"
            }
            """)
            conv_crit_settings_disp = KratosMultiphysics.Parameters("""
            {
                "basis_vector_type" : "solution_update"
            }
            """)

            if(convergence_crit == "displacement_criterion"):
                self.mechanical_convergence_criterion = StructuralMechanicsApplication.GeneralConvergenceCriteria(D_RT, D_AT, conv_crit_settings_disp)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)

            elif(convergence_crit == "residual_criterion"):
                self.mechanical_convergence_criterion = StructuralMechanicsApplication.GeneralConvergenceCriteria(R_RT, R_AT, conv_crit_settings_res)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)

            elif(convergence_crit == "and_criterion"):
                Displacement = StructuralMechanicsApplication.GeneralConvergenceCriteria(D_RT, D_AT, conv_crit_settings_disp)
                Displacement.SetEchoLevel(echo_level)
                Residual = StructuralMechanicsApplication.GeneralConvergenceCriteria(R_RT, R_AT, conv_crit_settings_res)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)

            elif(convergence_crit == "or_criterion"):
                Displacement = StructuralMechanicsApplication.GeneralConvergenceCriteria(D_RT, D_AT, conv_crit_settings_disp)
                Displacement.SetEchoLevel(echo_level)
                Residual = StructuralMechanicsApplication.GeneralConvergenceCriteria(R_RT, R_AT, conv_crit_settings_res)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
            else:
                err_msg =  "The requested convergence criterion \"" + convergence_crit + "\" is not available!\n"
                err_msg += "Available options are: \"displacement_criterion\", \"residual_criterion\", \"and_criterion\", \"or_criterion\""
                raise Exception(err_msg)


