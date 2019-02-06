# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
import KratosMultiphysics as KM
from KratosRANSConstitutiveLawsApplication import *
application = KratosRANSConstitutiveLawsApplication()
application_name = "KratosRANSConstitutiveLawsApplication"
application_folder = "RANSConstitutiveLawsApplication"

from .. import application_importer
application_importer.ImportApplication(application, application_name, application_folder, __path__)
