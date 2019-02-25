kratos_install_path = "/home/mchung/Kratos"
kratos_python_source_path = "/home/mchung/Kratos"

import os
kratos_libs = os.path.abspath(os.path.join(kratos_install_path, "libs"))
kratos_applications = os.path.abspath(os.path.join(kratos_python_source_path, "applications"))
kratos_scripts = os.path.abspath(os.path.join(kratos_python_source_path, "kratos/python_scripts"))
kratos_tests = os.path.abspath(os.path.join(kratos_python_source_path, "kratos/tests"))
