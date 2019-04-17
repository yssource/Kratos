# Import Kratos core and apps
import KratosMultiphysics as km
import KratosMultiphysics.ShapeOptimizationApplication as kso


su2_interface_settings = km.Parameters("""
{
    "su2_related": {
        "config_file"                 : "ubend_coarse.cfg",
        "number_of_cores"             : 12,
        "gradient_method"             : "DISCRETE_ADJOINT",
        "write_frequency_adjoint_run" : 500,
        "residual_min_adjoint_run"    : -10,
        "mesh_file"                   : "csm_model.su2",
        "mesh_motion_file"            : "mesh_motion.dat",
        "use_restart_in_adjoint_run"  : "YES",
        "design_surface_tag"          : "wet_interface"
    },
    "kratos_related": {
        "mdpa_file"      : "csm_model.mdpa",
        "write_elements" : true
    },
    "echo_level" : 2
}""")

from interface_su2 import InterfaceSU2
interface_su2 = InterfaceSU2(su2_interface_settings)
interface_su2.WriteSU2MeshAsMDPA()