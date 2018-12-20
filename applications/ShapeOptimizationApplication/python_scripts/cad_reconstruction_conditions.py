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

        local_rhs = -self.penalty_fac * (np.outer(self.shape_functions, self.shape_functions) @ pole_coords - np.outer(self.shape_functions, self.fe_node_coords))

        return local_rhs.T.flatten()

# ==============================================================================
class TangentEnforcementCondition:
    # --------------------------------------------------------------------------
    def __init__(self, target_normal, surface_geometry, nonzero_pole_indices, shape_function_derivatives_u, shape_function_derivatives_v, penalty_factor):
        self.target_normal = target_normal
        self.surface_geometry = surface_geometry
        self.nonzero_pole_indices = nonzero_pole_indices
        self.shape_function_derivatives_u = shape_function_derivatives_u
        self.shape_function_derivatives_v = shape_function_derivatives_v
        self.penalty_factor = penalty_factor

        self.local_system_size = 3*len(nonzero_pole_indices)
        self.block_size = len(nonzero_pole_indices)

    # --------------------------------------------------------------------------
    def CalculateLHS(self):
        term1 = np.outer(self.shape_function_derivatives_u,self.target_normal)
        term1 = term1.T.flatten()

        term2 = np.outer(self.shape_function_derivatives_v,self.target_normal)
        term2 = term2.T.flatten()

        return self.penalty_factor * np.outer(term1,term1) + self.penalty_factor * np.outer(term2,term2)

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        pole_coords = np.zeros((self.block_size, 3))
        for i, (r,s) in enumerate(self.nonzero_pole_indices):
            pole_coords[i,:] = self.surface_geometry.Pole(r,s)

        a1 = self.shape_function_derivatives_u @ pole_coords
        a2 = self.shape_function_derivatives_v @ pole_coords

        local_rhs_1 = - self.penalty_factor * np.outer(self.shape_function_derivatives_u,self.target_normal) * (a1 @ self.target_normal)
        local_rhs_2 = - self.penalty_factor * np.outer(self.shape_function_derivatives_v,self.target_normal) * (a2 @ self.target_normal)

        return local_rhs_1.T.flatten() + local_rhs_2.T.flatten()

# ==============================================================================
class PositionEnforcementCondition:
    # --------------------------------------------------------------------------
    def __init__(self, target_disp, target_position, surface_geometry, nonzero_pole_indices, shape_functions, penalty_factor):
        self.target_position = target_position
        self.target_disp = target_disp
        self.surface_geometry = surface_geometry
        self.nonzero_pole_indices = nonzero_pole_indices
        self.shape_functions = shape_functions
        self.penalty_fac = penalty_factor

        self.local_system_size = 3*len(nonzero_pole_indices)
        self.block_size = len(nonzero_pole_indices)

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

        local_rhs = -self.penalty_fac * (np.outer(self.shape_functions, self.shape_functions) @ pole_coords - np.outer(self.shape_functions, self.target_position))

        return local_rhs.T.flatten()


# ==============================================================================