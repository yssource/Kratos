from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.MappingApplication as KratosMapping

# Convergence criteria class
class convergence_criterion:
    def __init__(self, convergence_criterion_parameters):
        # Note that all the convergence settings are introduced via a Kratos parameters object.

        D_RT = convergence_criterion_parameters["displacement_relative_tolerance"].GetDouble()
        D_AT = convergence_criterion_parameters["displacement_absolute_tolerance"].GetDouble()
        R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
        R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()

        echo_level = 1#convergence_criterion_parameters["echo_level"].GetInt()
        convergence_crit = convergence_criterion_parameters["convergence_criterion"].GetString()

        if(echo_level >= 1):
            KratosMultiphysics.Logger.PrintInfo("::[Mechanical Solver]:: ", "CONVERGENCE CRITERION : " +
                  convergence_criterion_parameters["convergence_criterion"].GetString())

        convergence_crit_settings = KratosMultiphysics.Parameters("""{
            "basis_vector_type"     : "residual",
            "variables_to_separate" : [],
            "relative_tolerances"   : [],
            "absolute_tolerances"   : [],
            "other_dofs_name"       : "ROTATION",
            "print_colors"          : true
        }""");

        rotation_dofs = False
        if convergence_criterion_parameters.Has("rotation_dofs"):
            if convergence_criterion_parameters["rotation_dofs"].GetBool():
                convergence_crit_settings["variables_to_separate"].Append("DISPLACEMENT")
                rotation_dofs = True

        if convergence_crit == "displacement_criterion":
            convergence_crit_settings["basis_vector_type"].SetString("solution_update")
            convergence_crit_settings["relative_tolerances"].Append(D_RT)
            convergence_crit_settings["absolute_tolerances"].Append(D_AT)
            if rotation_dofs: # for now using the same tolerance for disp and rot
                convergence_crit_settings["relative_tolerances"].Append(D_RT)
                convergence_crit_settings["absolute_tolerances"].Append(D_AT)

        else:
            convergence_crit_settings["relative_tolerances"].Append(R_RT)
            convergence_crit_settings["absolute_tolerances"].Append(R_AT)
            if rotation_dofs: # for now using the same tolerance for disp and rot
                convergence_crit_settings["relative_tolerances"].Append(R_RT)
                convergence_crit_settings["absolute_tolerances"].Append(R_AT)

        print(convergence_crit_settings.PrettyPrintJsonString())
        print("starting....")

        if(convergence_crit == "displacement_criterion"):
            self.mechanical_convergence_criterion = KratosMapping.GeneralConvergenceCriteria(convergence_crit_settings)
            self.mechanical_convergence_criterion.SetEchoLevel(echo_level)

        elif(convergence_crit == "residual_criterion"):
            self.mechanical_convergence_criterion = KratosMapping.GeneralConvergenceCriteria(convergence_crit_settings)
            self.mechanical_convergence_criterion.SetEchoLevel(echo_level)

        elif(convergence_crit == "and_criterion"):
            err
            Displacement = StructuralMechanicsApplication.DisplacementAndOtherDoFCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)
            Residual = StructuralMechanicsApplication.ResidualDisplacementAndOtherDoFCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)

        elif(convergence_crit == "or_criterion"):
            err
            Displacement = StructuralMechanicsApplication.DisplacementAndOtherDoFCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)
            Residual = StructuralMechanicsApplication.ResidualDisplacementAndOtherDoFCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            self.mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
        else:
            err_msg =  "The requested convergence criterion \"" + convergence_crit + "\" is not available!\n"
            err_msg += "Available options are: \"displacement_criterion\", \"residual_criterion\", \"and_criterion\", \"or_criterion\""
            raise Exception(err_msg)
