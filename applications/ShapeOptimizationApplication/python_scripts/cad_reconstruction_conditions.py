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
class DistanceMinimizationCondition:
    # --------------------------------------------------------------------------
    def __init__(self, fe_node, surface_geometry, nonzero_pole_indices, shape_functions, variabl_to_map, penalty_fac):
        self.fe_node = fe_node
        self.surface_geometry = surface_geometry
        self.nonzero_pole_indices = nonzero_pole_indices
        self.shape_functions = shape_functions
        self.variabl_to_map = variabl_to_map

        # Penalty factor to overweight special nodes / conditions (e.g. ones on the boundary)
        self.penalty_fac = penalty_fac

        self.local_system_size = 3*len(nonzero_pole_indices)
        self.block_size = len(nonzero_pole_indices)
        node_initial_coords = np.array([self.fe_node.X, self.fe_node.Y, self.fe_node.Z])
        nodal_update = np.array(self.fe_node.GetSolutionStepValue(self.variabl_to_map))
        self.fe_node_coords = node_initial_coords + nodal_update

    # --------------------------------------------------------------------------
    def CalculateLHS(self):
        local_lhs = np.zeros([self.local_system_size,self.local_system_size])

        local_lhs[0:self.block_size,0:self.block_size] = self.penalty_fac * np.outer(self.shape_functions, self.shape_functions)
        local_lhs[self.block_size:2*self.block_size,self.block_size:2*self.block_size] = local_lhs[0:self.block_size,0:self.block_size]
        local_lhs[2*self.block_size:3*self.block_size,2*self.block_size:3*self.block_size] = local_lhs[0:self.block_size,0:self.block_size]

        return local_lhs

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        pole_coords = np.zeros((self.block_size, 3))
        for i, (r,s) in enumerate(self.nonzero_pole_indices):
            pole_coords[i,:] = self.surface_geometry.Pole(r,s)

        local_rhs = -self.penalty_fac * (np.outer(self.shape_functions, self.shape_functions)@pole_coords - np.outer(self.shape_functions, self.fe_node_coords))

        return local_rhs.T.flatten()

# ==============================================================================