from __future__ import print_function, absolute_import, division
import os
import sys

# Application dependent names and paths
import KratosMultiphysics as KM
from KratosCoSimulationApplication import *
application = KratosCoSimulationApplication()
application_name = 'KratosCoSimulationApplication'
application_folder = 'CoSimulationApplication'

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)
## adding pyKratos to path
application_path = os.path.join("/media/aditya/750Disk/01_Software/01_Kratos/03_co_sim_new/Kratos/applications", application_folder)
python_path = os.path.join(application_path, 'custom_data_structure')
sys.path.append(python_path)
