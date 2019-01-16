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

# For automatic differentiation
from HyperJet import HyperJet, Jet
from numpy.testing import assert_array_almost_equal

# ==============================================================================
class DistanceMinimizationCondition:
    # --------------------------------------------------------------------------
    def __init__(self, fe_node, surface_geometry, nonzero_pole_indices, shape_functions, variabl_to_map, weight):
        self.fe_node = fe_node
        self.surface_geometry = surface_geometry
        self.nonzero_pole_indices = nonzero_pole_indices
        self.shape_functions = shape_functions
        self.variabl_to_map = variabl_to_map
        self.weight = weight

        self.local_system_size = 3*len(nonzero_pole_indices)
        self.block_size = len(nonzero_pole_indices)
        node_initial_coords = np.array([self.fe_node.X, self.fe_node.Y, self.fe_node.Z])
        nodal_update = np.array(self.fe_node.GetSolutionStepValue(self.variabl_to_map))
        self.fe_node_coords = node_initial_coords + nodal_update

    # --------------------------------------------------------------------------
    def CalculateLHS(self):
        local_lhs = np.zeros([self.local_system_size,self.local_system_size])

        local_lhs[0:self.block_size,0:self.block_size] = self.weight * np.outer(self.shape_functions, self.shape_functions)
        local_lhs[self.block_size:2*self.block_size,self.block_size:2*self.block_size] = local_lhs[0:self.block_size,0:self.block_size]
        local_lhs[2*self.block_size:3*self.block_size,2*self.block_size:3*self.block_size] = local_lhs[0:self.block_size,0:self.block_size]

        return local_lhs

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        pole_coords = np.zeros((self.block_size, 3))
        for i, (r,s) in enumerate(self.nonzero_pole_indices):
            pole_coords[i,:] = self.surface_geometry.Pole(r,s)

        local_rhs = -self.weight * np.outer(self.shape_functions, (self.shape_functions @ pole_coords - self.fe_node_coords))

        return local_rhs.T.flatten()

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

        local_rhs = -self.penalty_fac * np.outer(self.shape_functions, (self.shape_functions @ pole_coords - self.target_position))

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
class ReconstructionConditionWithADBase():
    # --------------------------------------------------------------------------
    def ComputeActual(self, shape):
        value = np.zeros(3)

        for i, (r,s) in enumerate(self.nonzero_pole_indices):
            value += shape[i]*self.surface_geometry.Pole(r,s)

        return value

    # --------------------------------------------------------------------------
    def ComputeActualJet(self, shape):
        value = np.zeros(3)

        for i, (r,s) in enumerate(self.nonzero_pole_indices):
            value += shape[i]*self.surface_geometry.Pole(r,s)

        dx = np.zeros(shape.size * 3)
        dy = np.zeros(shape.size * 3)
        dz = np.zeros(shape.size * 3)

        dx[:shape.size] = shape
        dy[shape.size:2*shape.size] = shape
        dz[2*shape.size:] = shape

        return np.array([Jet(value[0], dx),
                         Jet(value[1], dy),
                         Jet(value[2], dz)])

    # --------------------------------------------------------------------------
    def ComputeActualHyperJet(self, shape):
        value = np.zeros(3)

        for i, (r,s) in enumerate(self.nonzero_pole_indices):
            value += shape[i]*self.surface_geometry.Pole(r,s)

        dx = np.zeros(shape.size * 3)
        dy = np.zeros(shape.size * 3)
        dz = np.zeros(shape.size * 3)

        dx[:shape.size] = shape
        dy[shape.size:2*shape.size] = shape
        dz[2*shape.size:] = shape

        return np.array([HyperJet(value[0], dx),
                         HyperJet(value[1], dy),
                         HyperJet(value[2], dz)])

# ==============================================================================
class TangentEnforcementConditionWithAD( ReconstructionConditionWithADBase ):
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
        a1 = super().ComputeActualHyperJet(self.shape_function_derivatives_u)
        a2 = super().ComputeActualHyperJet(self.shape_function_derivatives_v)

        f_1 = self.penalty_factor * np.dot(a1,self.target_normal)**2 / 2
        f_2 = self.penalty_factor * np.dot(a2,self.target_normal)**2 / 2

        local_lhs_1 = f_1.h
        local_lhs_2 = f_2.h

        lhs_new = local_lhs_1 + local_lhs_2

        # # Compare against old way of computing it
        # term1 = np.outer(self.shape_function_derivatives_u,self.target_normal)
        # term1 = term1.T.flatten()

        # term2 = np.outer(self.shape_function_derivatives_v,self.target_normal)
        # term2 = term2.T.flatten()

        # lhs_old = self.penalty_factor * np.outer(term1,term1) + self.penalty_factor * np.outer(term2,term2)

        # assert_array_almost_equal(lhs_new,lhs_old,10)

        return lhs_new

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        a1 = super().ComputeActualJet(self.shape_function_derivatives_u)
        a2 = super().ComputeActualJet(self.shape_function_derivatives_v)

        f_1 = self.penalty_factor * np.dot(a1,self.target_normal)**2 / 2
        f_2 = self.penalty_factor * np.dot(a2,self.target_normal)**2 / 2

        local_rhs_1 = - f_1.g
        local_rhs_2 = - f_2.g

        rhs_new = local_rhs_1 + local_rhs_2

        # # Compare against old way of computing it
        # pole_coords = np.zeros((self.block_size, 3))
        # for i, (r,s) in enumerate(self.nonzero_pole_indices):
        #     pole_coords[i,:] = self.surface_geometry.Pole(r,s)

        # a1 = self.shape_function_derivatives_u @ pole_coords
        # a2 = self.shape_function_derivatives_v @ pole_coords

        # local_rhs_1 = - self.penalty_factor * np.outer(self.shape_function_derivatives_u,self.target_normal) * (a1 @ self.target_normal)
        # local_rhs_2 = - self.penalty_factor * np.outer(self.shape_function_derivatives_v,self.target_normal) * (a2 @ self.target_normal)

        # rhs_old = local_rhs_1.T.flatten() + local_rhs_2.T.flatten()

        # assert_array_almost_equal(rhs_new,rhs_old,10)

        return rhs_new

# ==============================================================================
class CurvatureMinimizationConditionWithAD( ReconstructionConditionWithADBase ):
    # --------------------------------------------------------------------------
    def __init__(self, surface_geometry, nonzero_pole_indices, shape_function_derivatives_u, shape_function_derivatives_v, shape_function_derivatives_uu, shape_function_derivatives_uv, shape_function_derivatives_vv, penalty_factor):
        self.surface_geometry = surface_geometry
        self.nonzero_pole_indices = nonzero_pole_indices
        self.shape_function_derivatives_u = shape_function_derivatives_u
        self.shape_function_derivatives_v = shape_function_derivatives_v
        self.shape_function_derivatives_uu = shape_function_derivatives_uu
        self.shape_function_derivatives_uv = shape_function_derivatives_uv
        self.shape_function_derivatives_vv = shape_function_derivatives_vv
        self.penalty_factor = penalty_factor

        self.local_system_size = 3*len(nonzero_pole_indices)
        self.block_size = len(nonzero_pole_indices)

        # Reference configuration
        A1 = super().ComputeActual(shape_function_derivatives_u)
        A2 = super().ComputeActual(shape_function_derivatives_v)
        A1_1 = super().ComputeActual(shape_function_derivatives_uu)
        A1_2 = super().ComputeActual(shape_function_derivatives_uv)
        A2_2 = super().ComputeActual(shape_function_derivatives_vv)

        A3 = np.cross(A1,A2)
        self.dA = np.linalg.norm(A3)
        A3 = A3/self.dA

        self.B11 = np.dot(A1_1,A3)
        self.B12 = np.dot(A1_2,A3)
        self.B22 = np.dot(A2_2,A3)

    # --------------------------------------------------------------------------
    def CalculateLHS(self):
        a1 = super().ComputeActualHyperJet(self.shape_function_derivatives_u)
        a2 = super().ComputeActualHyperJet(self.shape_function_derivatives_v)
        a1_1 = super().ComputeActualHyperJet(self.shape_function_derivatives_uu)
        a1_2 = super().ComputeActualHyperJet(self.shape_function_derivatives_uv)
        a2_2 = super().ComputeActualHyperJet(self.shape_function_derivatives_vv)

        a3 = np.cross(a1,a2)
        a3 /= np.linalg.norm(a3)

        b11 = np.dot(a1_1,a3)
        b12 = np.dot(a1_2,a3)
        b22 = np.dot(a2_2,a3)

        f_1 = self.penalty_factor * (self.B11 - b11)**2 / 2
        f_2 = self.penalty_factor * (self.B12 - b12)**2 / 2
        f_3 = self.penalty_factor * (self.B22 - b22)**2 / 2

        local_lhs_1 = f_1.h
        local_lhs_2 = f_2.h
        local_lhs_3 = f_3.h

        return local_lhs_1 + local_lhs_2 + local_lhs_3

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        a1 = super().ComputeActualJet(self.shape_function_derivatives_u)
        a2 = super().ComputeActualJet(self.shape_function_derivatives_v)
        a1_1 = super().ComputeActualJet(self.shape_function_derivatives_uu)
        a1_2 = super().ComputeActualJet(self.shape_function_derivatives_uv)
        a2_2 = super().ComputeActualJet(self.shape_function_derivatives_vv)

        a3 = np.cross(a1,a2)
        a3 /= np.linalg.norm(a3)

        b11 = np.dot(a1_1,a3)
        b12 = np.dot(a1_2,a3)
        b22 = np.dot(a2_2,a3)

        f_1 = self.penalty_factor * (self.B11 - b11)**2 / 2
        f_2 = self.penalty_factor * (self.B12 - b12)**2 / 2
        f_3 = self.penalty_factor * (self.B22 - b22)**2 / 2

        local_rhs_1 = - f_1.g
        local_rhs_2 = - f_2.g
        local_rhs_3 = - f_3.g

        return local_rhs_1 + local_rhs_2 + local_rhs_3

# ==============================================================================
class KLShellConditionWithAD( ReconstructionConditionWithADBase ):
    # --------------------------------------------------------------------------
    def __init__(self, surface_geometry, nonzero_pole_indices, shape_function_derivatives_u, shape_function_derivatives_v, shape_function_derivatives_uu, shape_function_derivatives_uv, shape_function_derivatives_vv, penalty_factor):
        self.surface_geometry = surface_geometry
        self.nonzero_pole_indices = nonzero_pole_indices
        self.shape_function_derivatives_u = shape_function_derivatives_u
        self.shape_function_derivatives_v = shape_function_derivatives_v
        self.shape_function_derivatives_uu = shape_function_derivatives_uu
        self.shape_function_derivatives_uv = shape_function_derivatives_uv
        self.shape_function_derivatives_vv = shape_function_derivatives_vv
        self.penalty_factor = penalty_factor

        self.local_system_size = 3*len(nonzero_pole_indices)
        self.block_size = len(nonzero_pole_indices)

        # Reference configuration
        A1 = super().ComputeActual(shape_function_derivatives_u)
        A2 = super().ComputeActual(shape_function_derivatives_v)
        A1_1 = super().ComputeActual(shape_function_derivatives_uu)
        A1_2 = super().ComputeActual(shape_function_derivatives_uv)
        A2_2 = super().ComputeActual(shape_function_derivatives_vv)

        A11, A22, A12 = np.dot(A1, A1), np.dot(A2, A2), np.dot(A1, A2)
        self.A11, self.A22, self.A12 = A11, A22, A12

        A3 = np.cross(A1,A2)
        self.dA = np.linalg.norm(A3)
        A3 = A3/self.dA

        B11, B12, B22 = np.dot([A1_1, A1_2, A2_2], A3)
        self.B11, self.B12, self.B22 = B11, B12, B22

        # Material properties
        thickness = 1.0
        youngsmodulus = 1.0
        poissonsratio = 0.3

        self.Dm = np.array([
            [1.0, poissonsratio, 0],
            [poissonsratio, 1.0, 0],
            [0, 0, (1.0 - poissonsratio) / 2.0],
        ]) * youngsmodulus * thickness / (1.0 - np.power(poissonsratio,2))

        self.Db = np.array([
            [1.0, poissonsratio, 0],
            [poissonsratio, 1.0, 0],
            [0, 0, (1.0 - poissonsratio) / 2.0],
        ]) * youngsmodulus * np.power(thickness, 3) / (12.0 * (1.0 - np.power(poissonsratio,2)))

        # Transformation
        e1 = A1 / np.linalg.norm(A1)
        e2 = A2 - np.dot(A2, e1) * e1
        e2 /= np.linalg.norm(e2)

        det = A11 * A22 - A12 * A12

        g_ab_con = np.array([A22 / det, A11 / det, -A12 / det])

        g_con1 = g_ab_con[0] * A1 + g_ab_con[2] * A2
        g_con2 = g_ab_con[2] * A1 + g_ab_con[1] * A2

        eg11 = np.dot(e1, g_con1)
        eg12 = np.dot(e1, g_con2)
        eg21 = np.dot(e2, g_con1)
        eg22 = np.dot(e2, g_con2)

        self.Tm = np.array([
            [eg11 * eg11, eg12 * eg12, 2 * eg11 * eg12],
            [eg21 * eg21, eg22 * eg22, 2 * eg21 * eg22],
            [2 * eg11 * eg21, 2 * eg12 * eg22, 2 * (eg11 * eg22 + eg12 * eg21)],
        ])

    # --------------------------------------------------------------------------
    def CalculateLHS(self):
        B11, B12, B22 = self.B11, self.B12, self.B22
        A11, A22, A12 = self.A11, self.A22, self.A12

        a1 =super().ComputeActualHyperJet(self.shape_function_derivatives_u)
        a2 =super().ComputeActualHyperJet(self.shape_function_derivatives_v)
        a1_1 =super().ComputeActualHyperJet(self.shape_function_derivatives_uu)
        a1_2 =super().ComputeActualHyperJet(self.shape_function_derivatives_uv)
        a2_2 =super().ComputeActualHyperJet(self.shape_function_derivatives_vv)

        a3 = np.cross(a1,a2)
        a3 /= np.linalg.norm(a3)

        a11, a22, a12 = np.dot(a1, a1), np.dot(a2, a2), np.dot(a1, a2)
        b11, b12, b22 = np.dot([a1_1, a1_2, a2_2], a3)

        eps = np.dot(self.Tm, [a11 - A11, a22 - A22, a12 - A12]) / 2
        kap = np.dot(self.Tm, [B11 - b11, B12 - b12, B22 - b22])

        n = np.dot(self.Dm, [j.f for j in eps])
        m = np.dot(self.Db, [j.f for j in kap])

        k_act = np.zeros((self.local_system_size, self.local_system_size))
        for k in range(3):
            for l in range(3):
                k_act += np.outer(eps[k].g, eps[l].g) * self.Dm[k, l]
                k_act += np.outer(kap[k].g, kap[l].g) * self.Db[k, l]

            k_act += n[k] * eps[k].h + m[k] * kap[k].h

        return self.penalty_factor * self.dA * k_act

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        B11, B12, B22 = self.B11, self.B12, self.B22
        A11, A22, A12 = self.A11, self.A22, self.A12

        a1 =super().ComputeActualHyperJet(self.shape_function_derivatives_u)
        a2 =super().ComputeActualHyperJet(self.shape_function_derivatives_v)
        a1_1 =super().ComputeActualHyperJet(self.shape_function_derivatives_uu)
        a1_2 =super().ComputeActualHyperJet(self.shape_function_derivatives_uv)
        a2_2 =super().ComputeActualHyperJet(self.shape_function_derivatives_vv)

        a3 = np.cross(a1,a2)
        a3 /= np.linalg.norm(a3)

        a11, a22, a12 = np.dot(a1, a1), np.dot(a2, a2), np.dot(a1, a2)
        b11, b12, b22 = np.dot([a1_1, a1_2, a2_2], a3)

        eps = np.dot(self.Tm, [a11 - A11, a22 - A22, a12 - A12]) / 2
        kap = np.dot(self.Tm, [B11 - b11, B12 - b12, B22 - b22])

        n = np.dot(self.Dm, [j.f for j in eps])
        m = np.dot(self.Db, [j.f for j in kap])

        f_act = np.zeros(self.local_system_size)
        for k in range(3):
            f_act -= n[k] * eps[k].g + m[k] * kap[k].g

        return self.penalty_factor * self.dA * f_act

# ==============================================================================