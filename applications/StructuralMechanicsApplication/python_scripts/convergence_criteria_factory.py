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

        echo_level = 1#convergence_criterion_parameters["echo_level"].GetInt()
        convergence_crit = convergence_criterion_parameters["convergence_criterion"].GetString()
        convergence_criterion_parameters.RemoveValue("convergence_criterion")

        if(echo_level >= 1):
            KratosMultiphysics.Logger.PrintInfo("::[Mechanical Solver]:: ", "CONVERGENCE CRITERION : " +
                  convergence_crit)

        print(convergence_criterion_parameters.PrettyPrintJsonString())
        print("starting....")

        if(convergence_crit == "residual_criterion"):
            convergence_criterion_parameters.AddEmptyValue("basis_vector_type").SetString("residual")
            self.mechanical_convergence_criterion = KratosMapping.GeneralConvergenceCriteria(convergence_criterion_parameters)
            self.mechanical_convergence_criterion.SetEchoLevel(echo_level)

        elif(convergence_crit == "displacement_criterion"):
            convergence_criterion_parameters.AddEmptyValue("basis_vector_type").SetString("solution_update")
            self.mechanical_convergence_criterion = KratosMapping.GeneralConvergenceCriteria(convergence_criterion_parameters)
            self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
        else:
            errrrr

        # elif(convergence_crit == "and_criterion"):
        #     err
        #     Displacement = StructuralMechanicsApplication.DisplacementAndOtherDoFCriteria(D_RT, D_AT)
        #     Displacement.SetEchoLevel(echo_level)
        #     Residual = StructuralMechanicsApplication.ResidualDisplacementAndOtherDoFCriteria(R_RT, R_AT)
        #     Residual.SetEchoLevel(echo_level)
        #     self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)

        # elif(convergence_crit == "or_criterion"):
        #     err
        #     Displacement = StructuralMechanicsApplication.DisplacementAndOtherDoFCriteria(D_RT, D_AT)
        #     Displacement.SetEchoLevel(echo_level)
        #     Residual = StructuralMechanicsApplication.ResidualDisplacementAndOtherDoFCriteria(R_RT, R_AT)
        #     Residual.SetEchoLevel(echo_level)
        #     self.mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
        # else:
        #     err_msg =  "The requested convergence criterion \"" + convergence_crit + "\" is not available!\n"
        #     err_msg += "Available options are: \"displacement_criterion\", \"residual_criterion\", \"and_criterion\", \"or_criterion\""
        #     raise Exception(err_msg)
