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
class ReconstructionCondition:
    # --------------------------------------------------------------------------
    @staticmethod
    def ComputeActual(geometry_data, nonzero_pole_ids, shape):
        value = np.zeros(3)

        for i, (r,s) in enumerate(nonzero_pole_ids):
            value += shape[i]*geometry_data.Pole(r,s)

        return value

    # --------------------------------------------------------------------------
    @staticmethod
    def GenerateDofList(surface_geometry_key, nonzero_pole_ids):
        dof_list = []
        for (r,s) in nonzero_pole_ids:
            dof_list.append( (surface_geometry_key,r,s,"x") )
        for (r,s) in nonzero_pole_ids:
            dof_list.append( (surface_geometry_key,r,s,"y") )
        for (r,s) in nonzero_pole_ids:
            dof_list.append( (surface_geometry_key,r,s,"z") )
        return dof_list

# ==============================================================================
class ReconstructionConditionWithAD(ReconstructionCondition):
    # --------------------------------------------------------------------------
    @staticmethod
    def ComputeActualJet(geometry_data, nonzero_pole_ids, shape):
        value = np.zeros(3)

        for i, (r,s) in enumerate(nonzero_pole_ids):
            value += shape[i]*geometry_data.Pole(r,s)

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
    @staticmethod
    def ComputeActualHyperJet(geometry_data, nonzero_pole_ids, shape):
        value = np.zeros(3)

        for i, (r,s) in enumerate(nonzero_pole_ids):
            value += shape[i]*geometry_data.Pole(r,s)

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
class DistanceMinimizationCondition(ReconstructionCondition):
    # --------------------------------------------------------------------------
    def __init__(self, fe_node, surface_geometry, nonzero_pole_ids, shape_functions, variabl_to_map, weight):
        self.fe_node = fe_node
        self.geometry_data = surface_geometry.Data()
        self.nonzero_pole_ids = nonzero_pole_ids
        self.shape_functions = shape_functions
        self.variabl_to_map = variabl_to_map
        self.weight = weight

        node_initial_coords = np.array([self.fe_node.X, self.fe_node.Y, self.fe_node.Z])
        nodal_update = np.array(self.fe_node.GetSolutionStepValue(self.variabl_to_map))
        self.fe_node_coords = node_initial_coords + nodal_update

        self.dof_list = super().GenerateDofList(surface_geometry.Key(), nonzero_pole_ids)
        self.num_local_dofs = len(self.dof_list)
        self.block_size = len(nonzero_pole_ids)

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        pole_coords = np.zeros((self.block_size, 3))
        for i, (r,s) in enumerate(self.nonzero_pole_ids):
            pole_coords[i,:] = self.geometry_data.Pole(r,s)

        local_rhs = -self.weight * np.outer(self.shape_functions, (self.shape_functions @ pole_coords - self.fe_node_coords))
        local_rhs = local_rhs.T.flatten()

        return local_rhs, self.dof_list

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self):
        # LHS
        local_lhs = np.zeros([self.num_local_dofs,self.num_local_dofs])

        local_lhs_block = self.weight * np.outer(self.shape_functions, self.shape_functions)
        local_lhs[0*self.block_size:1*self.block_size,0*self.block_size:1*self.block_size] = local_lhs_block
        local_lhs[1*self.block_size:2*self.block_size,1*self.block_size:2*self.block_size] = local_lhs_block
        local_lhs[2*self.block_size:3*self.block_size,2*self.block_size:3*self.block_size] = local_lhs_block

        # RHS
        local_rhs, _ = self.CalculateRHS()

        return local_lhs, local_rhs, self.dof_list

# ==============================================================================
class PositionEnforcementCondition(ReconstructionCondition):
    # --------------------------------------------------------------------------
    def __init__(self, target_position, surface_geometry, nonzero_pole_ids, shape_functions, penalty_factor):
        self.target_position = target_position
        self.geometry_data = surface_geometry.Data()
        self.nonzero_pole_ids = nonzero_pole_ids
        self.shape_functions = shape_functions
        self.penalty_fac = penalty_factor

        self.dof_list = super().GenerateDofList(surface_geometry.Key(), nonzero_pole_ids)
        self.num_local_dofs = len(self.dof_list)
        self.block_size = len(nonzero_pole_ids)

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        pole_coords = np.zeros((self.block_size, 3))
        for i, (r,s) in enumerate(self.nonzero_pole_ids):
            pole_coords[i,:] = self.geometry_data.Pole(r,s)

        local_rhs = -self.penalty_fac * np.outer(self.shape_functions, (self.shape_functions @ pole_coords - self.target_position))
        local_rhs = local_rhs.T.flatten()

        return local_rhs, self.dof_list

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self):
        # LHS
        local_lhs = np.zeros([self.num_local_dofs,self.num_local_dofs])

        local_lhs_block = self.penalty_fac * np.outer(self.shape_functions, self.shape_functions)
        local_lhs[0*self.block_size:1*self.block_size,0*self.block_size:1*self.block_size] = local_lhs_block
        local_lhs[1*self.block_size:2*self.block_size,1*self.block_size:2*self.block_size] = local_lhs_block
        local_lhs[2*self.block_size:3*self.block_size,2*self.block_size:3*self.block_size] = local_lhs_block

        # RHS
        local_rhs, _ = self.CalculateRHS()

        return local_lhs, local_rhs, self.dof_list

# ==============================================================================
class TangentEnforcementCondition(ReconstructionCondition):
    # --------------------------------------------------------------------------
    def __init__(self, target_normal, surface_geometry, nonzero_pole_ids, shape_function_derivatives_u, shape_function_derivatives_v, penalty_factor):
        self.target_normal = target_normal
        self.geometry_data = surface_geometry.Data()
        self.nonzero_pole_ids = nonzero_pole_ids
        self.shape_function_derivatives_u = shape_function_derivatives_u
        self.shape_function_derivatives_v = shape_function_derivatives_v
        self.penalty_factor = penalty_factor

        self.dof_list = super().GenerateDofList(surface_geometry.Key(), nonzero_pole_ids)
        self.block_size = len(nonzero_pole_ids)

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        pole_coords = np.zeros((self.block_size, 3))
        for i, (r,s) in enumerate(self.nonzero_pole_ids):
            pole_coords[i,:] = self.geometry_data.Pole(r,s)

        a1 = self.shape_function_derivatives_u @ pole_coords
        a2 = self.shape_function_derivatives_v @ pole_coords

        local_rhs_1 = - self.penalty_factor * np.outer(self.shape_function_derivatives_u,self.target_normal) * (a1 @ self.target_normal)
        local_rhs_2 = - self.penalty_factor * np.outer(self.shape_function_derivatives_v,self.target_normal) * (a2 @ self.target_normal)

        local_rhs = local_rhs_1.T.flatten() + local_rhs_2.T.flatten()

        return local_rhs, self.dof_list

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self):
        # LHS
        term1 = np.outer(self.shape_function_derivatives_u,self.target_normal)
        term1 = term1.T.flatten()

        term2 = np.outer(self.shape_function_derivatives_v,self.target_normal)
        term2 = term2.T.flatten()

        local_lhs = self.penalty_factor * np.outer(term1,term1) + self.penalty_factor * np.outer(term2,term2)

        # RHS
        local_rhs, _ = self.CalculateRHS()

        return local_lhs, local_rhs, self.dof_list

# ==============================================================================
class TangentEnforcementConditionWithAD( ReconstructionConditionWithAD ):
    # --------------------------------------------------------------------------
    def __init__(self, target_normal, surface_geometry, nonzero_pole_ids, shape_function_derivatives_u, shape_function_derivatives_v, penalty_factor):
        self.target_normal = target_normal
        self.geometry_data = surface_geometry.Data()
        self.nonzero_pole_ids = nonzero_pole_ids
        self.shape_function_derivatives_u = shape_function_derivatives_u
        self.shape_function_derivatives_v = shape_function_derivatives_v
        self.penalty_factor = penalty_factor

        self.dof_list = super().GenerateDofList(surface_geometry.Key(), nonzero_pole_ids)
        self.block_size = len(nonzero_pole_ids)

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        a1 = super().ComputeActualJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_u)
        a2 = super().ComputeActualJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_v)

        f_1 = self.penalty_factor * np.dot(a1,self.target_normal)**2 * 0.5
        f_2 = self.penalty_factor * np.dot(a2,self.target_normal)**2 * 0.5

        local_rhs_1 = - f_1.g
        local_rhs_2 = - f_2.g

        rhs_new = local_rhs_1 + local_rhs_2

        # # Compare against old way of computing it
        # pole_coords = np.zeros((self.block_size, 3))
        # for i, (r,s) in enumerate(self.nonzero_pole_ids):
        #     pole_coords[i,:] = self.surface_geometry.Pole(r,s)

        # a1 = self.shape_function_derivatives_u @ pole_coords
        # a2 = self.shape_function_derivatives_v @ pole_coords

        # local_rhs_1 = - self.penalty_factor * np.outer(self.shape_function_derivatives_u,self.target_normal) * (a1 @ self.target_normal)
        # local_rhs_2 = - self.penalty_factor * np.outer(self.shape_function_derivatives_v,self.target_normal) * (a2 @ self.target_normal)

        # rhs_old = local_rhs_1.T.flatten() + local_rhs_2.T.flatten()

        # assert_array_almost_equal(rhs_new,rhs_old,10)


        # # Compare against old way of computing it
        # term1 = np.outer(self.shape_function_derivatives_u,self.target_normal)
        # term1 = term1.T.flatten()

        # term2 = np.outer(self.shape_function_derivatives_v,self.target_normal)
        # term2 = term2.T.flatten()

        # lhs_old = self.penalty_factor * np.outer(term1,term1) + self.penalty_factor * np.outer(term2,term2)

        # assert_array_almost_equal(lhs_new,lhs_old,10)

        return rhs_new, self.dof_list

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self):
        a1 = super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_u)
        a2 = super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_v)

        f_1 = self.penalty_factor * np.dot(a1,self.target_normal)**2 * 0.5
        f_2 = self.penalty_factor * np.dot(a2,self.target_normal)**2 * 0.5

        # LHS
        local_lhs_1 = f_1.h
        local_lhs_2 = f_2.h

        lhs_new = local_lhs_1 + local_lhs_2

        # RHS
        local_rhs_1 = - f_1.g
        local_rhs_2 = - f_2.g

        rhs_new = local_rhs_1 + local_rhs_2

        return lhs_new, rhs_new, self.dof_list

# ==============================================================================
class DisplacementCouplingCondition(ReconstructionCondition):
    # --------------------------------------------------------------------------
    def __init__(self, geometry_a, geometry_b, parameters_a, parameters_b, nonzero_pole_ids_a, nonzero_pole_ids_b, shape_functions_a, shape_functions_b, penalty_factor):
        self.geometry_a = geometry_a
        self.geometry_b = geometry_b
        self.parameters_a = parameters_a
        self.parameters_b = parameters_b
        self.nonzero_pole_ids_a = nonzero_pole_ids_a
        self.nonzero_pole_ids_b = nonzero_pole_ids_b
        self.shape_functions_a = shape_functions_a
        self.shape_functions_b = shape_functions_b
        self.penalty_fac = penalty_factor

        dof_list_a = super().GenerateDofList(geometry_a.Key(), nonzero_pole_ids_a)
        dof_list_b = super().GenerateDofList(geometry_b.Key(), nonzero_pole_ids_b)

        self.dof_list = dof_list_a + dof_list_b

        self.num_local_dofs = len(self.dof_list)
        self.block_size_a = len(nonzero_pole_ids_a)
        self.block_size_b = len(nonzero_pole_ids_b)

        self.pos_0_a = super().ComputeActual(self.geometry_a.Data(), nonzero_pole_ids_a, shape_functions_a)
        self.pos_0_b = super().ComputeActual(self.geometry_b.Data(), nonzero_pole_ids_b, shape_functions_b)

    # --------------------------------------------------------------------------
    def CalculateQualityIndicator(self):
        disp_a = super().ComputeActual(self.geometry_a.Data(), self.nonzero_pole_ids_a, self.shape_functions_a) - self.pos_0_a
        disp_b = super().ComputeActual(self.geometry_b.Data(), self.nonzero_pole_ids_b, self.shape_functions_b) - self.pos_0_b
        return disp_a - disp_b

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        N_ab = np.zeros(self.num_local_dofs)

        N_ab[0*self.block_size_a:1*self.block_size_a] = self.shape_functions_a
        N_ab[1*self.block_size_a:2*self.block_size_a] = self.shape_functions_a
        N_ab[2*self.block_size_a:3*self.block_size_a] = self.shape_functions_a

        N_ab[3*self.block_size_a+0*self.block_size_b:3*self.block_size_a+1*self.block_size_b] -= self.shape_functions_b
        N_ab[3*self.block_size_a+1*self.block_size_b:3*self.block_size_a+2*self.block_size_b] -= self.shape_functions_b
        N_ab[3*self.block_size_a+2*self.block_size_b:3*self.block_size_a+3*self.block_size_b] -= self.shape_functions_b

        # RHS
        disp_a = super().ComputeActual(self.geometry_a.Data(), self.nonzero_pole_ids_a, self.shape_functions_a) - self.pos_0_a
        disp_b = super().ComputeActual(self.geometry_b.Data(), self.nonzero_pole_ids_b, self.shape_functions_b) - self.pos_0_b
        delta = disp_a - disp_b

        # a
        delta_ab = np.zeros(self.num_local_dofs)
        delta_ab[0*self.block_size_a:1*self.block_size_a] = delta[0]
        delta_ab[1*self.block_size_a:2*self.block_size_a] = delta[1]
        delta_ab[2*self.block_size_a:3*self.block_size_a] = delta[2]

        # b
        delta_ab[3*self.block_size_a+0*self.block_size_b:3*self.block_size_a+1*self.block_size_b] = delta[0]
        delta_ab[3*self.block_size_a+1*self.block_size_b:3*self.block_size_a+2*self.block_size_b] = delta[1]
        delta_ab[3*self.block_size_a+2*self.block_size_b:3*self.block_size_a+3*self.block_size_b] = delta[2]

        local_rhs = - self.penalty_fac * np.multiply(N_ab, delta_ab)

        return local_rhs, self.dof_list

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self):
        # LHS
        local_lhs = np.zeros([self.num_local_dofs,self.num_local_dofs])

        # aa
        local_lhs_aa = self.penalty_fac * np.outer(self.shape_functions_a, self.shape_functions_a)
        local_lhs[0*self.block_size_a:1*self.block_size_a , 0*self.block_size_a:1*self.block_size_a] = local_lhs_aa
        local_lhs[1*self.block_size_a:2*self.block_size_a , 1*self.block_size_a:2*self.block_size_a] = local_lhs_aa
        local_lhs[2*self.block_size_a:3*self.block_size_a , 2*self.block_size_a:3*self.block_size_a] = local_lhs_aa

        # ab
        local_lhs_ab = - self.penalty_fac * np.outer(self.shape_functions_a, self.shape_functions_b)
        local_lhs[0*self.block_size_a:1*self.block_size_a , 3*self.block_size_a+0*self.block_size_b:3*self.block_size_a+1*self.block_size_b] = local_lhs_ab
        local_lhs[1*self.block_size_a:2*self.block_size_a , 3*self.block_size_a+1*self.block_size_b:3*self.block_size_a+2*self.block_size_b] = local_lhs_ab
        local_lhs[2*self.block_size_a:3*self.block_size_a , 3*self.block_size_a+2*self.block_size_b:3*self.block_size_a+3*self.block_size_b] = local_lhs_ab

        # ba
        local_lsh_ba = - self.penalty_fac * np.outer(self.shape_functions_b, self.shape_functions_a)
        local_lhs[3*self.block_size_a+0*self.block_size_b:3*self.block_size_a+1*self.block_size_b , 0*self.block_size_a:1*self.block_size_a] = local_lsh_ba
        local_lhs[3*self.block_size_a+1*self.block_size_b:3*self.block_size_a+2*self.block_size_b , 1*self.block_size_a:2*self.block_size_a] = local_lsh_ba
        local_lhs[3*self.block_size_a+2*self.block_size_b:3*self.block_size_a+3*self.block_size_b , 2*self.block_size_a:3*self.block_size_a] = local_lsh_ba

        # bb
        local_lhs_bb = self.penalty_fac * np.outer(self.shape_functions_b, self.shape_functions_b)
        local_lhs[3*self.block_size_a+0*self.block_size_b:3*self.block_size_a+1*self.block_size_b , 3*self.block_size_a+0*self.block_size_b:3*self.block_size_a+1*self.block_size_b] = local_lhs_bb
        local_lhs[3*self.block_size_a+1*self.block_size_b:3*self.block_size_a+2*self.block_size_b , 3*self.block_size_a+1*self.block_size_b:3*self.block_size_a+2*self.block_size_b] = local_lhs_bb
        local_lhs[3*self.block_size_a+2*self.block_size_b:3*self.block_size_a+3*self.block_size_b , 3*self.block_size_a+2*self.block_size_b:3*self.block_size_a+3*self.block_size_b] = local_lhs_bb

        # RHS
        local_rhs, _ = self.CalculateRHS()

        return local_lhs, local_rhs, self.dof_list

# ==============================================================================
class DisplacementCouplingConditionWithAD(ReconstructionConditionWithAD):
    # --------------------------------------------------------------------------
    def __init__(self, geometry_a, geometry_b, nonzero_pole_ids_a, nonzero_pole_ids_b, shape_functions_a, shape_functions_b, penalty_factor):
        self.geometry_a = geometry_a
        self.geometry_b = geometry_b
        self.nonzero_pole_ids_a = nonzero_pole_ids_a
        self.nonzero_pole_ids_b = nonzero_pole_ids_b
        self.shape_functions_a = shape_functions_a
        self.shape_functions_b = shape_functions_b
        self.penalty_fac = penalty_factor

        dof_list_a = super().GenerateDofList(geometry_a.Key(), nonzero_pole_ids_a)
        dof_list_b = super().GenerateDofList(geometry_b.Key(), nonzero_pole_ids_b)

        self.dof_list = dof_list_a + dof_list_b

        self.num_local_dofs = len(self.dof_list)
        self.block_size_a = len(nonzero_pole_ids_a)
        self.block_size_b = len(nonzero_pole_ids_b)

        self.pos_0_a = super().ComputeActual(self.geometry_a.Data(), nonzero_pole_ids_a, shape_functions_a)
        self.pos_0_b = super().ComputeActual(self.geometry_b.Data(), nonzero_pole_ids_b, shape_functions_b)

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        pos_a = super().ComputeActualJet(self.geometry_a.Data(), self.nonzero_pole_ids_a, self.shape_functions_a)
        pos_b = super().ComputeActualJet(self.geometry_b.Data(), self.nonzero_pole_ids_b, self.shape_functions_b)

        len_a = len(pos_a[0])
        len_b = len(pos_b[0])

        pos_a = np.array([hj.enlarge(len_b, False) for hj in pos_a])
        pos_b = np.array([hj.enlarge(len_a, True ) for hj in pos_b])

        position_difference = pos_a - pos_b

        function = 0.5  * self.penalty_fac * np.dot(position_difference,position_difference)

        local_rhs = -function.g

        return local_rhs, self.dof_list

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self):
        pos_a = super().ComputeActualHyperJet(self.geometry_a.Data(), self.nonzero_pole_ids_a, self.shape_functions_a)
        pos_b = super().ComputeActualHyperJet(self.geometry_b.Data(), self.nonzero_pole_ids_b, self.shape_functions_b)

        len_a = len(pos_a[0])
        len_b = len(pos_b[0])

        pos_a = np.array([hj.enlarge(len_b, False) for hj in pos_a])
        pos_b = np.array([hj.enlarge(len_a, True ) for hj in pos_b])

        position_difference = pos_a - pos_b

        function = 0.5 * self.penalty_fac * np.dot(position_difference,position_difference)

        local_rhs = -function.g
        local_lhs = function.h

        return local_lhs, local_rhs, self.dof_list

# ==============================================================================
class RotationCouplingConditionWithAD(ReconstructionConditionWithAD):
    # --------------------------------------------------------------------------
    def __init__(self, geometry_a, geometry_b, parameters_a, parameters_b, T2_edge, nonzero_pole_ids_a, nonzero_pole_ids_b, shape_functions_a, shape_functions_b, shape_function_derivatives_u_a, shape_function_derivatives_u_b, shape_function_derivatives_v_a, shape_function_derivatives_v_b, penalty_factor):
        self.geometry_a = geometry_a
        self.geometry_b = geometry_b
        self.parameters_a = parameters_a
        self.parameters_b = parameters_b
        self.T2_edge = T2_edge
        self.nonzero_pole_ids_a = nonzero_pole_ids_a
        self.nonzero_pole_ids_b = nonzero_pole_ids_b
        self.shape_functions_a = shape_functions_a
        self.shape_functions_b = shape_functions_b
        self.shape_function_derivatives_u_a = shape_function_derivatives_u_a
        self.shape_function_derivatives_v_a = shape_function_derivatives_v_a
        self.shape_function_derivatives_u_b = shape_function_derivatives_u_b
        self.shape_function_derivatives_v_b = shape_function_derivatives_v_b
        self.penalty_fac = penalty_factor

        dof_list_a = super().GenerateDofList(geometry_a.Key(), nonzero_pole_ids_a)
        dof_list_b = super().GenerateDofList(geometry_b.Key(), nonzero_pole_ids_b)

        self.dof_list = dof_list_a + dof_list_b

        self.num_local_dofs = len(self.dof_list)
        self.block_size_a = len(nonzero_pole_ids_a)
        self.block_size_b = len(nonzero_pole_ids_b)

        A1_a = super().ComputeActual(self.geometry_a.Data(), self.nonzero_pole_ids_a, self.shape_function_derivatives_u_a)
        A2_a = super().ComputeActual(self.geometry_a.Data(), self.nonzero_pole_ids_a, self.shape_function_derivatives_v_a)
        A3_a = np.cross(A1_a, A2_a)
        self.A3_a = A3_a / np.linalg.norm(A3_a)

        A1_b = super().ComputeActual(self.geometry_b.Data(), self.nonzero_pole_ids_b, self.shape_function_derivatives_u_b)
        A2_b = super().ComputeActual(self.geometry_b.Data(), self.nonzero_pole_ids_b, self.shape_function_derivatives_v_b)
        A3_b = np.cross(A1_b, A2_b)
        self.A3_b = A3_b / np.linalg.norm(A3_b)

    # --------------------------------------------------------------------------
    def CalculateQualityIndicator(self):
        a1_a = super().ComputeActual(self.geometry_a.Data(), self.nonzero_pole_ids_a, self.shape_function_derivatives_u_a)
        a2_a = super().ComputeActual(self.geometry_a.Data(), self.nonzero_pole_ids_a, self.shape_function_derivatives_v_a)
        a3_a = np.cross(a1_a, a2_a)
        a3_a = a3_a / np.linalg.norm(a3_a)

        a1_b = super().ComputeActual(self.geometry_b.Data(), self.nonzero_pole_ids_b, self.shape_function_derivatives_u_b)
        a2_b = super().ComputeActual(self.geometry_b.Data(), self.nonzero_pole_ids_b, self.shape_function_derivatives_v_b)
        a3_b = np.cross(a1_b, a2_b)
        a3_b = a3_b / np.linalg.norm(a3_b)

        w_a = a3_a - self.A3_a
        w_b = a3_b - self.A3_b

        omega_a = np.cross(self.A3_a, w_a)
        omega_b = np.cross(self.A3_b, w_b)

        angle_a = np.arcsin(np.dot(omega_a, self.T2_edge))
        angle_b = np.arcsin(np.dot(omega_b, self.T2_edge))

        return angle_a - angle_b

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        a1_a = super().ComputeActualJet(self.geometry_a.Data(), self.nonzero_pole_ids_a, self.shape_function_derivatives_u_a)
        a2_a = super().ComputeActualJet(self.geometry_a.Data(), self.nonzero_pole_ids_a, self.shape_function_derivatives_v_a)
        a3_a = np.cross(a1_a, a2_a)
        a3_a = a3_a / np.linalg.norm(a3_a)

        a1_b = super().ComputeActualJet(self.geometry_b.Data(), self.nonzero_pole_ids_b, self.shape_function_derivatives_u_b)
        a2_b = super().ComputeActualJet(self.geometry_b.Data(), self.nonzero_pole_ids_b, self.shape_function_derivatives_v_b)
        a3_b = np.cross(a1_b, a2_b)
        a3_b = a3_b / np.linalg.norm(a3_b)

        w_a = a3_a - self.A3_a
        w_b = a3_b - self.A3_b

        omega_a = np.cross(self.A3_a, w_a)
        omega_b = np.cross(self.A3_b, w_b)

        angle_a = np.arcsin(np.dot(omega_a, self.T2_edge))
        angle_b = np.arcsin(np.dot(omega_b, self.T2_edge))

        angle_a = angle_a.enlarge(len(angle_a), False)
        angle_b = angle_b.enlarge(len(angle_b), True)

        angular_difference = angle_a - angle_b

        angle = 0.5 * self.penalty_fac * angular_difference**2

        local_rhs = -angle.g

        return local_rhs, self.dof_list

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self):
        a1_a = super().ComputeActualHyperJet(self.geometry_a.Data(), self.nonzero_pole_ids_a, self.shape_function_derivatives_u_a)
        a2_a = super().ComputeActualHyperJet(self.geometry_a.Data(), self.nonzero_pole_ids_a, self.shape_function_derivatives_v_a)
        a3_a = np.cross(a1_a, a2_a)
        a3_a = a3_a / np.linalg.norm(a3_a)

        a1_b = super().ComputeActualHyperJet(self.geometry_b.Data(), self.nonzero_pole_ids_b, self.shape_function_derivatives_u_b)
        a2_b = super().ComputeActualHyperJet(self.geometry_b.Data(), self.nonzero_pole_ids_b, self.shape_function_derivatives_v_b)
        a3_b = np.cross(a1_b, a2_b)
        a3_b = a3_b / np.linalg.norm(a3_b)

        w_a = a3_a - self.A3_a
        w_b = a3_b - self.A3_b

        omega_a = np.cross(self.A3_a, w_a)
        omega_b = np.cross(self.A3_b, w_b)

        angle_a = np.arcsin(np.dot(omega_a, self.T2_edge))
        angle_b = np.arcsin(np.dot(omega_b, self.T2_edge))

        angle_a = angle_a.enlarge(len(angle_a), False)
        angle_b = angle_b.enlarge(len(angle_b), True)

        angular_difference = angle_a - angle_b

        angle = 0.5 * self.penalty_fac * angular_difference**2

        local_rhs = -angle.g
        local_lhs = angle.h

        return local_lhs, local_rhs, self.dof_list

# ==============================================================================
class CurvatureMinimizationConditionWithAD( ReconstructionConditionWithAD ):
    # --------------------------------------------------------------------------
    def __init__(self, surface_geometry, nonzero_pole_ids, shape_function_derivatives_u, shape_function_derivatives_v, shape_function_derivatives_uu, shape_function_derivatives_uv, shape_function_derivatives_vv, penalty_factor):
        self.geometry_data = surface_geometry.Data()
        self.nonzero_pole_ids = nonzero_pole_ids
        self.shape_function_derivatives_u = shape_function_derivatives_u
        self.shape_function_derivatives_v = shape_function_derivatives_v
        self.shape_function_derivatives_uu = shape_function_derivatives_uu
        self.shape_function_derivatives_uv = shape_function_derivatives_uv
        self.shape_function_derivatives_vv = shape_function_derivatives_vv
        self.penalty_factor = penalty_factor

        self.dof_list = super().GenerateDofList(surface_geometry.Key(), nonzero_pole_ids)
        self.block_size = len(nonzero_pole_ids)

        # Reference configuration
        A1 = super().ComputeActual(self.geometry_data, nonzero_pole_ids, shape_function_derivatives_u)
        A2 = super().ComputeActual(self.geometry_data, nonzero_pole_ids, shape_function_derivatives_v)
        A1_1 = super().ComputeActual(self.geometry_data, nonzero_pole_ids, shape_function_derivatives_uu)
        A1_2 = super().ComputeActual(self.geometry_data, nonzero_pole_ids, shape_function_derivatives_uv)
        A2_2 = super().ComputeActual(self.geometry_data, nonzero_pole_ids, shape_function_derivatives_vv)

        A3 = np.cross(A1,A2)
        self.dA = np.linalg.norm(A3)
        A3 = A3/self.dA

        self.B11 = np.dot(A1_1,A3)
        self.B12 = np.dot(A1_2,A3)
        self.B22 = np.dot(A2_2,A3)

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        a1 = super().ComputeActualJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_u)
        a2 = super().ComputeActualJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_v)
        a1_1 = super().ComputeActualJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_uu)
        a1_2 = super().ComputeActualJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_uv)
        a2_2 = super().ComputeActualJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_vv)

        a3 = np.cross(a1,a2)
        a3 /= np.linalg.norm(a3)

        b11 = np.dot(a1_1,a3)
        b12 = np.dot(a1_2,a3)
        b22 = np.dot(a2_2,a3)

        f_1 = self.penalty_factor * (self.B11 - b11)**2 * 0.5
        f_2 = self.penalty_factor * (self.B12 - b12)**2 * 0.5
        f_3 = self.penalty_factor * (self.B22 - b22)**2 * 0.5

        local_rhs_1 = - f_1.g
        local_rhs_2 = - f_2.g
        local_rhs_3 = - f_3.g

        local_rhs = local_rhs_1 + local_rhs_2 + local_rhs_3

        return local_rhs, self.dof_list

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self):
        a1 = super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_u)
        a2 = super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_v)
        a1_1 = super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_uu)
        a1_2 = super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_uv)
        a2_2 = super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_vv)

        a3 = np.cross(a1,a2)
        a3 /= np.linalg.norm(a3)

        b11 = np.dot(a1_1,a3)
        b12 = np.dot(a1_2,a3)
        b22 = np.dot(a2_2,a3)

        f_1 = self.penalty_factor * (self.B11 - b11)**2 * 0.5
        f_2 = self.penalty_factor * (self.B12 - b12)**2 * 0.5
        f_3 = self.penalty_factor * (self.B22 - b22)**2 * 0.5

        # LHS
        local_lhs_1 = f_1.h
        local_lhs_2 = f_2.h
        local_lhs_3 = f_3.h

        local_lhs = local_lhs_1 + local_lhs_2 + local_lhs_3

        # RHS
        local_rhs_1 = - f_1.g
        local_rhs_2 = - f_2.g
        local_rhs_3 = - f_3.g

        local_rhs = local_rhs_1 + local_rhs_2 + local_rhs_3

        return local_lhs, local_rhs, self.dof_list

# ==============================================================================
class KLShellConditionWithAD( ReconstructionConditionWithAD ):
    # --------------------------------------------------------------------------
    def __init__(self, surface_geometry, nonzero_pole_ids, shape_function_derivatives_u, shape_function_derivatives_v, shape_function_derivatives_uu, shape_function_derivatives_uv, shape_function_derivatives_vv, penalty_factor):
        self.geometry_data = surface_geometry.Data()
        self.nonzero_pole_ids = nonzero_pole_ids
        self.shape_function_derivatives_u = shape_function_derivatives_u
        self.shape_function_derivatives_v = shape_function_derivatives_v
        self.shape_function_derivatives_uu = shape_function_derivatives_uu
        self.shape_function_derivatives_uv = shape_function_derivatives_uv
        self.shape_function_derivatives_vv = shape_function_derivatives_vv
        self.penalty_factor = penalty_factor

        self.dof_list = super().GenerateDofList(surface_geometry.Key(), nonzero_pole_ids)
        self.num_local_dofs = len(self.dof_list)
        self.block_size = len(nonzero_pole_ids)

        # Reference configuration
        A1 = super().ComputeActual(self.geometry_data, nonzero_pole_ids, shape_function_derivatives_u)
        A2 = super().ComputeActual(self.geometry_data, nonzero_pole_ids, shape_function_derivatives_v)
        A1_1 = super().ComputeActual(self.geometry_data, nonzero_pole_ids, shape_function_derivatives_uu)
        A1_2 = super().ComputeActual(self.geometry_data, nonzero_pole_ids, shape_function_derivatives_uv)
        A2_2 = super().ComputeActual(self.geometry_data, nonzero_pole_ids, shape_function_derivatives_vv)

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
            [0, 0, (1.0 - poissonsratio) * 0.5],
        ]) * youngsmodulus * thickness / (1.0 - np.power(poissonsratio,2))

        self.Db = np.array([
            [1.0, poissonsratio, 0],
            [poissonsratio, 1.0, 0],
            [0, 0, (1.0 - poissonsratio) * 0.5],
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
    def CalculateRHS(self):
        B11, B12, B22 = self.B11, self.B12, self.B22
        A11, A22, A12 = self.A11, self.A22, self.A12

        a1 =super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_u)
        a2 =super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_v)
        a1_1 =super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_uu)
        a1_2 =super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_uv)
        a2_2 =super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_vv)

        a3 = np.cross(a1,a2)
        a3 /= np.linalg.norm(a3)

        a11, a22, a12 = np.dot(a1, a1), np.dot(a2, a2), np.dot(a1, a2)
        b11, b12, b22 = np.dot([a1_1, a1_2, a2_2], a3)

        eps = np.dot(self.Tm, [a11 - A11, a22 - A22, a12 - A12]) * 0.5
        kap = np.dot(self.Tm, [B11 - b11, B12 - b12, B22 - b22])

        n = np.dot(self.Dm, [j.f for j in eps])
        m = np.dot(self.Db, [j.f for j in kap])

        f_act = np.zeros(self.num_local_dofs)
        for k in range(3):
            f_act -= n[k] * eps[k].g + m[k] * kap[k].g

        local_rhs = self.penalty_factor * self.dA * f_act

        return local_rhs, self.dof_list

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self):
        B11, B12, B22 = self.B11, self.B12, self.B22
        A11, A22, A12 = self.A11, self.A22, self.A12

        a1 =super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_u)
        a2 =super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_v)
        a1_1 =super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_uu)
        a1_2 =super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_uv)
        a2_2 =super().ComputeActualHyperJet(self.geometry_data, self.nonzero_pole_ids, self.shape_function_derivatives_vv)

        a3 = np.cross(a1,a2)
        a3 /= np.linalg.norm(a3)

        a11, a22, a12 = np.dot(a1, a1), np.dot(a2, a2), np.dot(a1, a2)
        b11, b12, b22 = np.dot([a1_1, a1_2, a2_2], a3)

        eps = np.dot(self.Tm, [a11 - A11, a22 - A22, a12 - A12]) * 0.5
        kap = np.dot(self.Tm, [B11 - b11, B12 - b12, B22 - b22])

        n = np.dot(self.Dm, [j.f for j in eps])
        m = np.dot(self.Db, [j.f for j in kap])

        f_act = np.zeros(self.num_local_dofs)
        k_act = np.zeros((self.num_local_dofs, self.num_local_dofs))
        for k in range(3):
            f_act -= n[k] * eps[k].g + m[k] * kap[k].g
            for l in range(3):
                k_act += np.outer(eps[k].g, eps[l].g) * self.Dm[k, l]
                k_act += np.outer(kap[k].g, kap[l].g) * self.Db[k, l]

            k_act += n[k] * eps[k].h + m[k] * kap[k].h

        local_lhs = self.penalty_factor * self.dA * k_act
        local_rhs = self.penalty_factor * self.dA * f_act

        return local_lhs, local_rhs, self.dof_list

# ==============================================================================
class AlphaRegularizationCondition(ReconstructionCondition):
    # --------------------------------------------------------------------------
    def __init__(self, target_pole_indices, greville_params, surface_geometry, nonzero_pole_ids_of_greville_point, shape_functions, alpha):
        self.target_pole_indices = target_pole_indices
        self.greville_params = greville_params
        self.geometry_data = surface_geometry.Data()
        self.nonzero_pole_ids = nonzero_pole_ids_of_greville_point
        self.alpha = alpha

        self.dof_list = super().GenerateDofList(surface_geometry.Key(), nonzero_pole_ids_of_greville_point)
        self.num_local_dofs = len(self.dof_list)
        self.block_size = len(self.nonzero_pole_ids)

        I = np.zeros(self.block_size)
        tpole_local_system_id = nonzero_pole_ids_of_greville_point.index(target_pole_indices)
        I[tpole_local_system_id] = 1
        self.IminN =  (I - shape_functions)

    # --------------------------------------------------------------------------
    def CalculateRHS(self):
        pole_coords = self.geometry_data.Pole(self.target_pole_indices[0],self.target_pole_indices[1])
        greville_point_coords = self.geometry_data.PointAt(self.greville_params[0], self.greville_params[1])

        distance_vector = pole_coords-greville_point_coords
        distance_vector = np.array(distance_vector)

        # RHS
        local_rhs = np.zeros(self.num_local_dofs)
        local_rhs[0*self.block_size:1*self.block_size] = -self.alpha * self.IminN * distance_vector[0]
        local_rhs[1*self.block_size:2*self.block_size] = -self.alpha * self.IminN * distance_vector[1]
        local_rhs[2*self.block_size:3*self.block_size] = -self.alpha * self.IminN * distance_vector[2]

        return local_rhs, self.dof_list

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self):
        # LHS
        local_lhs = np.zeros([self.num_local_dofs,self.num_local_dofs])
        local_lhs_block = self.alpha * np.outer(self.IminN, self.IminN)
        local_lhs[0*self.block_size:1*self.block_size,0*self.block_size:1*self.block_size] = local_lhs_block
        local_lhs[1*self.block_size:2*self.block_size,1*self.block_size:2*self.block_size] = local_lhs_block
        local_lhs[2*self.block_size:3*self.block_size,2*self.block_size:3*self.block_size] = local_lhs_block

        # RHS
        local_rhs, _ = self.CalculateRHS()

        return local_lhs, local_rhs, self.dof_list

# ==============================================================================