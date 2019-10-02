from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Application dependent names and paths
import KratosMultiphysics
import KratosExaquteSandboxApplication
application = KratosExaquteSandboxApplication.KratosExaquteSandboxApplication()
application_name = "KratosExaquteSandboxApplication"
application_folder = "ExaquteSandboxApplication"

# The following lines are common for all applications
KratosMultiphysics._ImportApplicationAsModule(application, application_name, application_folder, __path__)