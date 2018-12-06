# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Import ANurbs library
import numpy as np

# ==============================================================================
class DisplacementMappingCondition:
    # --------------------------------------------------------------------------
    def __init__(self, fe_node, surface_geometry, nonzero_pole_indices, shape_functions, variabl_to_map):
        self.fe_node = fe_node
        self.surface_geometry = surface_geometry
        self.nonzero_pole_indices = nonzero_pole_indices
        self.shape_functions = shape_functions
        self.variabl_to_map = variabl_to_map

        self.local_system_size = len(nonzero_pole_indices)

    # --------------------------------------------------------------------------
    def CalculateLHS(self):
        local_lhs = np.zeros((self.local_system_size, self.local_system_size))

        for i in range(self.local_system_size):
            for j in range(self.local_system_size):
                local_lhs[i, j] = self.shape_functions[i] * self.shape_functions[j]

        return local_lhs

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        local_rhs = np.zeros((self.local_system_size, 3))
        node_coords = np.array([self.fe_node.X, self.fe_node.Y, self.fe_node.Z])
        nodal_update =  np.array(self.fe_node.GetSolutionStepValue(self.variabl_to_map))

        pole_coords = np.zeros((self.local_system_size, 3))
        for i, (r,s) in enumerate(self.nonzero_pole_indices):
            pole_coords[i,:] = self.surface_geometry.Pole(r,s)

        for i in range(self.local_system_size):

            fem_contribution = self.shape_functions[i] * (node_coords + nodal_update)

            cad_contribution = np.zeros(3)
            for j in range(self.local_system_size):
                cad_contribution += self.shape_functions[i] * self.shape_functions[j] * pole_coords[j]

            local_rhs[i] = -(cad_contribution - fem_contribution)

        return local_rhs

# ==============================================================================