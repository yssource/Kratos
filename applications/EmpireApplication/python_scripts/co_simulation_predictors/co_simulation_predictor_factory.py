from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

available_predictors = {
    "linear_derivative_based"    : "linear_derivative_based_predictor",
    "average_value_based"        : "average_value_based_predictor",
    "standard_linear"            : "standard_linear_predictor"
}

def CreatePredictor(predictor_settings, solvers, cosim_solver_details, level):
    """This function creates and returns the Predictor used for CoSimulation
    New Predictors have to be registered by adding them to "available_predictors"
    """
    if (type(predictor_settings) != dict):
        raise Exception("Input is expected to be provided as a python dictionary")

    predictor_type = predictor_settings["predictor_type"]

    if predictor_type in available_predictors:
        predictor_module = __import__(available_predictors[predictor_type])
        return predictor_module.Create(predictor_settings, solvers, cosim_solver_details, level)
    else:
        err_msg  = 'The requested Predictor "' + predictor_type + '" is not available!\n'
        err_msg += 'The following Predictors are available:\n'
        for avail_predictor in available_predictors:
            err_msg += "\t" + avail_predictor + "\n"
        raise NameError(err_msg)
