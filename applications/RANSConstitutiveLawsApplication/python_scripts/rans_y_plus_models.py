import KratosMultiphysics
import KratosMultiphysics.RANSConstitutiveLawsApplication as KratosRANS


def Factory(model_part, settings):
    if settings["model_type"].GetString() == "logarithmic":
        return KratosRANS.RansLogarithmicYPlusModelProcess(model_part, settings["model_settings"])
    else:
        raise Exception("Unknown y_plus model_type: " + settings["model_type"].GetString())